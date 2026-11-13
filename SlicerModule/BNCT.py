# BNCTModule.py
# Modulo unificado BNCT para 3D Slicer
# Combina: BNCT Dosimetry + BNCT DVH Analysis
# Carga manual del modelo STL fijo: SLT/conos_pared_hueco.stl (boton en la interfaz)

import os, json, logging
import vtk, qt, ctk, slicer
from vtk.util.numpy_support import numpy_to_vtk
from slicer.ScriptedLoadableModule import *
import numpy as np
import re
import sys
import unicodedata

#
# Utilidades compartidas
#

_MODULE_DIR = os.path.dirname(__file__)
_PROJECT_DIR = os.path.dirname(_MODULE_DIR)


def _first_existing(*paths):
    """Retorna el primer path que exista, o el primero de la lista si ninguno existe."""
    for p in paths:
        if os.path.exists(p):
            return p
    return paths[0]


def _normName(s):
    """Normaliza un nombre para emparejar segmentos con ROIs de la config.
    Mayusculas, sin acentos, separadores unificados a un espacio."""
    if s is None:
        return ""
    t = str(s).strip().upper()
    t = unicodedata.normalize("NFKD", t)
    t = "".join(c for c in t if not unicodedata.combining(c))
    for ch in "_-./":
        t = t.replace(ch, " ")
    t = re.sub(r"[^A-Z0-9 ]", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


# El STL vive dentro de SlicerModule/SLT
_STL_PATH = _first_existing(
    os.path.join(_MODULE_DIR, "SLT", "conos_pared_hueco.stl"),
)
# La config json vive dentro de SlicerModule/json
_ORGAN_CONFIG_PATH = _first_existing(
    os.path.join(_MODULE_DIR, "json", "organ_config.jsonc"),
)
_MODELS_PATH = _first_existing(
    os.path.join(_MODULE_DIR, "json", "models.jsonc"),
)
# El MAT vive dentro de SlicerModule/mat
_MAT_PATH = _first_existing(
    os.path.join(_MODULE_DIR, "mat", "info_rois.mat"),
)

# Prefijos de nombre usados para identificar los labelmaps (mascaras) creados por el modulo
_PREFIX_INDIVIDUAL = "Mascara_"
_NAME_COMBINED = "Mascaras_Combinadas"
_NAME_COMBINED_GRAY = "Mascaras_Combinadas_Gris"
_NAME_MASKS_FOLDER = "Mascaras_Organos"

# Nodo transform (traslacion pura) que centra el BODY del paciente en el origen RAS
_NAME_CENTER_TRANSFORM = "BNCT_Centrar_BODY"

# Coeficiente del eje X usado para el centrado basado en las dimensiones de la
# matriz de la imagen. Valor solicitado literalmente por el usuario
# (shift_x = N_x * 0.075); es configurable desde aqui si se requiere cambiarlo.
_CENTER_X_COEF = 0.075

# ------------------------------------------------------------------
# Motor radiobiologico UTCP (capa pura, sin Qt/slicer; testeable)
# ------------------------------------------------------------------
# El motor vive en SlicerModule/utcp/ (subcarpeta) para que el escaneo de
# modulos scripted de Slicer no intente instanciarlo como modulo de la GUI.
_UTCP_LIB_DIR = os.path.join(_MODULE_DIR, "utcp")
if _UTCP_LIB_DIR not in sys.path:
    sys.path.insert(0, _UTCP_LIB_DIR)
try:
    from utcp_models import (
        builtin_ntcp_models,
        builtin_tcp_models,
        default_isoe_model_for_roi,
    )
    from utcp_engine import compute_utcp_point
    from utcp_time_optimizer import sweep_utcp

    _UTCP_IMPORT_ERROR = None
except Exception as _e:  # pragma: no cover - dentro de Slicer
    _UTCP_IMPORT_ERROR = str(e)
    logging.warning("Motor UTCP no disponible: %s" % _e)


def load_jsonc(path):
    """Carga un archivo JSONC (JSON con comentarios // y /* */)."""
    with open(path) as f:
        text = f.read()
    text = re.sub(r"//.*?\n|/\*.*?\*/", "", text, flags=re.S)
    return json.loads(text)


def _log(msg):
    """Escribe a un archivo de debug local."""
    log_path = os.path.join(os.path.dirname(__file__), "_bnct_debug.txt")
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


_stl_cargado = False  # Bandera de modulo: garantiza que el STL se cargue UNA sola vez


def cargar_stl_fijo():
    """Carga el STL fijo del cono y pared del reactor en la escena de Slicer.
    Retorna el nodo modelo o None si falla. Idempotente: nunca carga dos veces."""
    global _stl_cargado
    # Doble proteccion: bandera de modulo + verificacion en la escena
    if _stl_cargado:
        return slicer.mrmlScene.GetFirstNodeByName("Conos_Pared_Reactor")
    existing = slicer.mrmlScene.GetFirstNodeByName("Conos_Pared_Reactor")
    if existing:
        _stl_cargado = True
        return existing
    if not os.path.exists(_STL_PATH):
        _log(f"STL no encontrado: {_STL_PATH}")
        return None
    try:
        # En Slicer 5.x loadModel devuelve el nodo directamente (returnNode deprecado)
        model_node = slicer.util.loadModel(_STL_PATH)
        if model_node is None:
            _log("loadModel devolvio None para el STL del reactor")
            return None
        model_node.SetName("Conos_Pared_Reactor")
        # Color gris claro para el reactor
        display = model_node.GetDisplayNode()
        if display:
            display.SetColor(0.7, 0.7, 0.7)
            display.SetOpacity(0.4)
            display.BackfaceCullingOff()
        _stl_cargado = True
        _log(f"STL cargado: {_STL_PATH}")
        # Encadenar slices y posicionar en slice 33
        _encadenarSlicesSlice33()
        return model_node
    except Exception as e:
        _log(f"Error cargando STL: {e}")
        return None


def _encadenarSlicesSlice33():
    """Vincula los 3 slice viewers (Red, Yellow, Green) y posiciona el corte en el slice 33 (mm).
    Los slices quedan encadenados: mover uno mueve los otros."""
    try:
        layoutManager = slicer.app.layoutManager()
        sliceNames = ["Red", "Yellow", "Green"]
        sliceNodes = []
        for name in sliceNames:
            sw = layoutManager.sliceWidget(name)
            if sw:
                sn = sw.mrmlSliceNode()
                if sn:
                    sliceNodes.append(sn)
        if not sliceNodes:
            _log("No se encontraron slice nodes para encadenar")
            return
        # Vincular todos los slices entre sí (el primero es el master)
        master = sliceNodes[0]
        master.SetLink(True)
        for sn in sliceNodes[1:]:
            sn.SetLink(True)
            sn.SetMasterSliceNodeID(master.GetID())
        # Posicionar en slice offset 33 mm en todos
        for sn in sliceNodes:
            sn.SetSliceOffset(33.0)
        _log("Slices encadenados y posicionados en slice offset 33 mm")
    except Exception as e:
        _log(f"Error encadenando slices: {e}")


# ------------------------------------------------------------------
# Fondo celeste de las vistas 2D (Red / Yellow / Green)
# ------------------------------------------------------------------

# Color celeste (r, g, b) en rango [0, 1] para el fondo de las vistas 2D.
_SLICE_BG_CELESTE = (0.7, 0.9, 1.0)

# Nombres de las vistas 2D de Slicer.
_SLICE_VIEW_NAMES = ("Red", "Yellow", "Green")


def _logCeleste(msg):
    """Log del modulo: archivo de debug local (_bnct_debug.txt) + consola de Slicer."""
    try:
        _log(msg)
    except Exception:
        pass
    try:
        logging.info(msg)
    except Exception:
        pass
    try:
        print(msg)
    except Exception:
        pass


def _callLater(callback, delayMs=0):
    """Ejecuta callback en el proximo ciclo del event loop de Qt.

    slicer.app.callLater no tiene la misma firma en todas las versiones de
    Slicer, por eso se prueban ambas variantes de argumentos y, si ninguna
    funciona, se usa QTimer.singleShot como respaldo."""
    try:
        app = slicer.app
    except Exception:
        app = None
    if app is not None:
        try:
            app.callLater(callback, delayMs)
            return True
        except TypeError:
            try:
                app.callLater(delayMs, callback)
                return True
            except Exception:
                pass
        except Exception:
            pass
    try:
        qt.QTimer.singleShot(delayMs, callback)
        return True
    except Exception as e:
        _logCeleste(f"BNCT: no se pudo programar la tarea diferida: {e}")
        return False


def _setSliceBackgroundCeleste(r=None, g=None, b=None):
    """Fija el fondo celeste de las 3 vistas 2D (Red, Yellow, Green).

    Usa el mismo patron que _encadenarSlicesSlice33(): layoutManager ->
    sliceWidget(name) -> mrmlSliceNode(). Como garantia visual tambien se
    fuerza el color en el renderer del slice view.

    Es idempotente: se puede llamar varias veces sin efectos secundarios.
    Devuelve True si se coloreo al menos una vista y False si todavia no hay
    layout manager o slice widgets disponibles (en ese caso NO lanza excepcion).
    """
    if r is None:
        r, g, b = _SLICE_BG_CELESTE
    colored = []
    errors = 0
    try:
        layoutManager = slicer.app.layoutManager()
        if layoutManager is None:
            _logCeleste(
                "BNCT: fondo celeste no aplicado todavia: layoutManager no disponible."
            )
            return False
        for name in _SLICE_VIEW_NAMES:
            try:
                sw = layoutManager.sliceWidget(name)
                if sw is None:
                    continue
                sn = sw.mrmlSliceNode()
                if sn is not None:
                    sn.SetBackgroundColor(r, g, b)
                # Garantia visual: fijar el color tambien en el renderer.
                sliceView = sw.sliceView()
                rw = None
                renderer = None
                if sliceView is not None:
                    rw = sliceView.renderWindow()
                    if rw is not None:
                        renderers = rw.GetRenderers()
                        if renderers is not None:
                            renderer = renderers.GetFirstRenderer()
                    if renderer is None:
                        # Respaldo: el slice view expone su renderer directamente.
                        getRenderer = getattr(sliceView, "renderer", None)
                        if callable(getRenderer):
                            renderer = getRenderer()
                if renderer is not None:
                    renderer.SetBackground(r, g, b)
                    forceRender = getattr(sliceView, "forceRender", None)
                    if callable(forceRender):
                        forceRender()
                    elif rw is not None:
                        rw.ForceRender()
                colored.append(name)
            except Exception as e:
                errors += 1
                _logCeleste(
                    f"BNCT: error aplicando fondo celeste a la vista {name}: {e}"
                )
        if colored:
            _logCeleste(
                "BNCT: fondo celeste aplicado a las vistas 2D: "
                + ", ".join(colored)
                + f" (errores: {errors})"
            )
            return True
        _logCeleste(
            "BNCT: fondo celeste no aplicado: no se encontraron slice widgets "
            "Red/Yellow/Green."
        )
        return False
    except Exception as e:
        _logCeleste(f"BNCT: error global aplicando fondo celeste: {e}")
        return False


def _scheduleSliceBackgroundCeleste(delayMs=0, reintentos=2):
    """Programa la aplicacion del fondo celeste en el proximo ciclo del event loop.

    En setup() el layout manager puede todavia no existir, por eso la rutina
    se difiere y se reintenta un par de veces si el layout aun no esta listo."""
    def _aplicar():
        ok = _setSliceBackgroundCeleste()
        if not ok and reintentos > 0:
            _scheduleSliceBackgroundCeleste(delayMs=300, reintentos=reintentos - 1)

    _callLater(_aplicar, delayMs)


class BNCT(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "BNCT"
        self.parent.categories = ["BNCT"]
        self.parent.dependencies = []
        self.parent.contributors = ["Antigravity AI / BNCTAr Team"]
        self.parent.helpText = """
        <b>Modulo unificado BNCT para 3D Slicer.</b><br>
        Combina dos funcionalidades:<br>
        1. <b>BNCT Dosimetry:</b> carga archivos PHITS, interpola dosis y calcula dosis fisica/RBE/isoefectiva.<br>
        2. <b>BNCT DVH Analysis:</b> calcula histogramas DVH por órgano.<br>
        3. <b>Mascaras:</b> genera labelmaps por ROI desde un volumen de estructura (gray_level) y verifica superposicion.<br>
        <br>
        Se carga manualmente el modelo STL fijo del cono y pared del reactor
        (SLT/conos_pared_hueco.stl) mediante el botón de la interfaz.
        """
        self.parent.acknowledgementText = "BNCTAr en 3D Slicer."


class BNCTWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # Botón manual para cargar STL del reactor (respaldo)
        stlBtnLayout = qt.QHBoxLayout()
        self.btnLoadSTL = qt.QPushButton("Cargar STL Reactor (conos_pared_hueco.stl)")
        self.btnLoadSTL.setToolTip(
            "Carga el modelo 3D fijo de la pared y cono de salida del reactor en la escena de Slicer."
        )
        self.btnLoadSTL.connect("clicked(bool)", self.onLoadSTL)
        stlBtnLayout.addWidget(self.btnLoadSTL)
        stlBtnLayout.addStretch()
        self.layout.addLayout(stlBtnLayout)

        # Estado del STL
        self.stlStatusLabel = qt.QLabel("STL: carga manual con el boton de abajo.")
        self.stlStatusLabel.setWordWrap(True)
        self.stlStatusLabel.setStyleSheet("color: #666; font-size: 11px;")
        self.layout.addWidget(self.stlStatusLabel)

        # ===== PESTAÑAS: Dosimetry | DVH =====
        self.tabs = qt.QTabWidget()
        self.layout.addWidget(self.tabs)

        # --- Pestaña 1: Dosimetry ---
        self.dosimetryTab = qt.QWidget()
        self.dosimetryLayout = qt.QVBoxLayout(self.dosimetryTab)
        self.tabs.addTab(self.dosimetryTab, "Dosimetry")
        self._setupDosimetryUI()

        # --- Pestaña 2: DVH Analysis ---
        self.dvhTab = qt.QWidget()
        self.dvhLayout = qt.QVBoxLayout(self.dvhTab)
        self.tabs.addTab(self.dvhTab, "DVH Analysis")
        self._setupDVHUI()

        # --- Pestaña 3: Mascaras ---
        self.masksTab = qt.QWidget()
        self.masksLayout = qt.QVBoxLayout(self.masksTab)
        self.tabs.addTab(self.masksTab, "Mascaras")
        self._setupMasksUI()

        # --- Pestana 4: Colisiones ---
        self.collisionsTab = qt.QWidget()
        self.collisionsLayout = qt.QVBoxLayout(self.collisionsTab)
        self.tabs.addTab(self.collisionsTab, "Colisiones")
        self._setupCollisionsUI()

        # --- Pestana 5: Voxelizado ---
        self.voxelTab = qt.QWidget()
        self.voxelLayout = qt.QVBoxLayout(self.voxelTab)
        self.tabs.addTab(self.voxelTab, "Voxelizado")
        self._setupVoxelUI()

        _log(
            "BNCT module loaded: v2026-11-06 json-only gray_levels + autocarga voxel + labels wrap"
        )
        self.layout.addStretch(1)

        # Fondo celeste de las vistas 2D (Red/Yellow/Green).
        # En setup() el layout manager puede todavia no existir, por eso la
        # aplicacion se difiere al proximo ciclo del event loop.
        _scheduleSliceBackgroundCeleste()

    def enter(self):
        """Se llama cada vez que el usuario abre el modulo: autocarga los niveles."""
        try:
            # Fondo celeste de las vistas 2D (Red/Yellow/Green).
            _setSliceBackgroundCeleste()
            if getattr(self, "voxelMappingTable", None) is None:
                self._removeLegacyMaskNodes()
                return
            removed = self._removeLegacyMaskNodes()
            if removed:
                self.voxelStatusLabel.setText(
                    f"Se eliminaron {len(removed)} mascaras de una version anterior "
                    "(niveles >=256).\nPresione 'Generar mascaras' para regenerarlas."
                )
                self.voxelStatusLabel.setStyleSheet("color: #800; font-size: 11px;")
            if self.voxelMappingTable.rowCount == 0:
                if self.voxelMaskSelector.currentNode() is not None:
                    self.onLoadVoxelLevels(silent=True)
        except Exception as e:
            _log(f"enter: autocarga de niveles fallo: {e}")

    def onLoadSTL(self):
        """Carga manual del STL del reactor (conos_pared_hueco.stl).
        Se ejecuta cuando el usuario presiona el botón, así la escena
        de Slicer ya está lista y puede aceptar nodos."""
        model_node = cargar_stl_fijo()
        if model_node:
            self.stlStatusLabel.text = f"STL cargado: {_STL_PATH} (modelo en escena)"
            self.stlStatusLabel.setStyleSheet("color: #080; font-size: 11px;")
            slicer.util.infoDisplay(
                f"Modelo 3D del reactor cargado: {os.path.basename(_STL_PATH)}",
                windowTitle="STL Reactor",
            )
        else:
            self.stlStatusLabel.text = "STL: error al cargar. Verifica que el archivo exista en SLT/conos_pared_hueco.stl"
            self.stlStatusLabel.setStyleSheet("color: #800; font-size: 11px;")
            slicer.util.errorDisplay(
                f"No se pudo cargar el STL: {_STL_PATH}", windowTitle="STL Reactor"
            )

    def onDetectCollisions(self):
        model1 = self.patientModelSelector.currentNode()
        model2 = self.coneModelSelector.currentNode()
        if not model1 or not model2:
            slicer.util.errorDisplay(
                "Seleccione el modelo del paciente y el modelo del cono/pared."
            )
            return

        logic = BNCTLogic()
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            result = logic.detectCollisions(model1, model2)
        finally:
            qt.QApplication.restoreOverrideCursor()

        if result is None:
            self.collisionStatusLabel.setText(
                "No se pudo analizar la colision. Ver _bnct_debug.txt."
            )
            self.collisionStatusLabel.setStyleSheet("color: #800; font-size: 11px;")
            self.collisionReportText.setPlainText(
                "Error al ejecutar la deteccion. Ver log."
            )
            slicer.util.errorDisplay(
                "No se pudo ejecutar la deteccion de colisiones. Ver el log."
            )
            return

        contacts = int(result.get("contacts", 0))
        collision = bool(result.get("collision", False))
        if collision:
            self.collisionStatusLabel.setText(
                f"COLISION: {contacts} contactos detectados."
            )
            self.collisionStatusLabel.setStyleSheet("color: #800; font-size: 11px;")
            self.collisionReportText.setPlainText(
                f"Modelo 1: {model1.GetName()}\n"
                f"Modelo 2: {model2.GetName()}\n"
                f"Contactos: {contacts}\n"
                f"Resultado: COLISION"
            )
            slicer.util.warningDisplay(
                f"Se detectaron {contacts} contactos entre los modelos."
            )
        else:
            self.collisionStatusLabel.setText("Sin interseccion detectada.")
            self.collisionStatusLabel.setStyleSheet("color: #080; font-size: 11px;")
            self.collisionReportText.setPlainText(
                f"Modelo 1: {model1.GetName()}\n"
                f"Modelo 2: {model2.GetName()}\n"
                f"Contactos: 0\n"
                f"Resultado: Sin interseccion"
            )
            slicer.util.infoDisplay("No se detectaron colisiones entre los modelos.")

    def onLoadVoxelLevels(self, silent=False):
        """Lee los niveles de gris unicos de la mascara y llena la tabla de mapeo.

        `silent=True` (autocarga automatica al abrir el tab o cambiar de mascara)
        evita el dialogo de error cuando todavia no hay mascara seleccionada.
        Tambien se usa como slot de clicked(bool): en ese caso `silent` llega
        como False y se mantiene el dialogo de error original.
        """
        maskNode = self.voxelMaskSelector.currentNode()
        if not maskNode:
            if silent is True:
                self.voxelMappingTable.clearContents()
                self.voxelMappingTable.setRowCount(0)
                self.voxelStatusLabel.setText(
                    "Sin mascara gris: seleccione Mascaras_Combinadas_Gris o genere las mascaras."
                )
                self.voxelStatusLabel.setStyleSheet("color: #666; font-size: 11px;")
                _log(
                    "onLoadVoxelLevels: sin mascara seleccionada (autocarga silenciosa)"
                )
                return
            slicer.util.errorDisplay(
                "Seleccione la mascara gris (Mascaras_Combinadas_Gris)."
            )
            return
        try:
            arr = slicer.util.arrayFromVolume(maskNode)
            if arr is None:
                if silent is not True:
                    slicer.util.errorDisplay(
                        "No se pudo leer el volumen de la mascara."
                    )
                else:
                    _log("onLoadVoxelLevels: volumen no leible (autocarga silenciosa)")
                return
            arr_int = np.rint(np.asarray(arr, dtype=np.float64)).astype(np.int64)
            levels = [int(v) for v in np.unique(arr_int)]

            defaults = self._voxelUniverseDefaults()
            sin_config = [gl for gl in levels if gl not in defaults]

            # Si la mascara seleccionada tiene niveles de una version anterior
            # (>=256) y existe otra mas nueva cuyos niveles estan todos en la
            # config, cambiar a esa automaticamente (una sola vez, sin recursion).
            if (
                sin_config
                and any(gl >= 256 for gl in sin_config)
                and not getattr(self, "_voxelAutoSwitching", False)
            ):
                clean = self._voxelFindCleanMask()
                if clean is not None and clean.GetID() != maskNode.GetID():
                    _log(
                        f"onLoadVoxelLevels: '{maskNode.GetName()}' tiene niveles >=256; "
                        f"cambiando a la mascara mas reciente '{clean.GetName()}'"
                    )
                    self._voxelAutoSwitching = True
                    try:
                        self.voxelMaskSelector.setCurrentNode(clean)
                    finally:
                        self._voxelAutoSwitching = False
                    self.onLoadVoxelLevels(silent=silent)
                    return

            self.voxelMappingTable.clearContents()
            self.voxelMappingTable.setRowCount(len(levels))
            for row, gl in enumerate(levels):
                if gl in defaults:
                    uni = defaults[gl]
                else:
                    uni = "SIN CONFIG"

                item_gl = qt.QTableWidgetItem(str(gl))
                item_gl.setFlags(qt.Qt.ItemIsEnabled | qt.Qt.ItemIsSelectable)
                self.voxelMappingTable.setItem(row, 0, item_gl)

                item_uni = qt.QTableWidgetItem(str(uni))
                self.voxelMappingTable.setItem(row, 1, item_uni)

            if sin_config:
                sin_config_sorted = sorted(sin_config)
                aviso = (
                    f"AVISO: {len(sin_config)} niveles sin universo:\n"
                    f"{sin_config_sorted}\n"
                    "Corrija la columna 'universo (U=)' o regenere las mascaras."
                )
                if any(gl >= 256 for gl in sin_config):
                    aviso += (
                        "\nNiveles >=256: la mascara seleccionada fue creada por una "
                        "version anterior del modulo.\n"
                        "Presione 'Generar mascaras': el modulo borra esas mascaras "
                        "viejas y crea una nueva con los gray_level de la config."
                    )
                self.voxelStatusLabel.setText(aviso)
                self.voxelStatusLabel.setStyleSheet("color: #800; font-size: 11px;")
                try:
                    _log(
                        f"onLoadVoxelLevels: nodo='{maskNode.GetName()}' id={maskNode.GetID()} "
                        f"max={int(np.max(arr_int))} niveles={sorted(levels)} "
                        f"sin_config={sin_config_sorted}"
                    )
                    for _i in range(
                        slicer.mrmlScene.GetNumberOfNodesByClass(
                            "vtkMRMLLabelMapVolumeNode"
                        )
                    ):
                        _n = slicer.mrmlScene.GetNthNodeByClass(
                            _i, "vtkMRMLLabelMapVolumeNode"
                        )
                        _nm = _n.GetName() or ""
                        if not (
                            _nm.startswith(_PREFIX_INDIVIDUAL)
                            or _nm == _NAME_COMBINED
                            or _nm.startswith(_NAME_COMBINED + "_")
                            or _nm == _NAME_COMBINED_GRAY
                            or _nm.startswith(_NAME_COMBINED_GRAY + "_")
                        ):
                            continue
                        _a = slicer.util.arrayFromVolume(_n)
                        _mx = int(np.max(np.asarray(_a))) if _a is not None else -1
                        _log(f"  mascara en escena: '{_nm}' id={_n.GetID()} max={_mx}")
                except Exception as _e:
                    _log(f"onLoadVoxelLevels: diagnostico fallo: {_e}")
                _log(f"onLoadVoxelLevels: niveles SIN CONFIG: {sin_config}")
            else:
                from_config = sum(1 for gl in levels if gl in defaults)
                self.voxelStatusLabel.setText(
                    f"{len(levels)} gray_levels cargados desde '{maskNode.GetName()}' "
                    f"({from_config} con universo desde organ_config.jsonc)."
                )
                self.voxelStatusLabel.setStyleSheet("color: #060; font-size: 11px;")
                _log(f"onLoadVoxelLevels: {len(levels)} niveles: {levels}")
        except Exception as e:
            _log(f"onLoadVoxelLevels: error: {e}")
            if silent is not True:
                slicer.util.errorDisplay(
                    "Error leyendo los niveles de gris. Ver _bnct_debug.txt."
                )

    def _onVoxelMaskChanged(self, node=None):
        """Autocarga la tabla de universos cada vez que cambia la mascara gris."""
        try:
            if node is None:
                self.voxelMappingTable.clearContents()
                self.voxelMappingTable.setRowCount(0)
                self.voxelStatusLabel.setText(
                    "Seleccione la mascara gris (Mascaras_Combinadas_Gris)."
                )
                self.voxelStatusLabel.setStyleSheet("color: #666; font-size: 11px;")
                return
            self.onLoadVoxelLevels(silent=True)
        except Exception as e:
            _log(f"_onVoxelMaskChanged: error: {e}")

    def _voxelFindCleanMask(self):
        """Devuelve la Mascaras_Combinadas_Gris mas nueva cuyos niveles estan
        todos en organ_config.jsonc (sin valores >=256 de versiones anteriores)."""
        defaults = self._voxelUniverseDefaults()
        candidates = []
        for i in range(
            slicer.mrmlScene.GetNumberOfNodesByClass("vtkMRMLLabelMapVolumeNode")
        ):
            n = slicer.mrmlScene.GetNthNodeByClass(i, "vtkMRMLLabelMapVolumeNode")
            nm = n.GetName() or ""
            if nm == _NAME_COMBINED_GRAY or nm.startswith(_NAME_COMBINED_GRAY + "_"):
                candidates.append(n)
        for node in reversed(candidates):  # el mas nuevo primero
            try:
                arr = slicer.util.arrayFromVolume(node)
                if arr is None:
                    continue
                vals = set(int(v) for v in np.unique(np.asarray(arr)))
                if all(v in defaults for v in vals):
                    return node
            except Exception as e:
                _log(f"_voxelFindCleanMask: error con '{node.GetName()}': {e}")
        return None

    def _removeLegacyMaskNodes(self):
        """Elimina labelmaps del modulo cuyos valores lleguen a 256 o mas.

        Esos niveles los inventaba una version anterior del modulo para los
        segmentos sin ROI en organ_config.jsonc; quedaron grabados dentro de la
        mascara combinada y no se pueden mapear a ningun universo (U=).
        Devuelve la lista de descripciones de los nodos eliminados.
        """
        removed = []
        try:
            n_nodes = slicer.mrmlScene.GetNumberOfNodesByClass(
                "vtkMRMLLabelMapVolumeNode"
            )
            for i in range(
                n_nodes - 1, -1, -1
            ):  # descendente: borrar no invalida indices menores
                node = slicer.mrmlScene.GetNthNodeByClass(
                    i, "vtkMRMLLabelMapVolumeNode"
                )
                if node is None:
                    continue
                nm = node.GetName() or ""
                # Solo los labelmaps COMBINADOS pueden tener valores >=256: las
                # mascaras individuales (Mascara_*) son binarias (0/1). Se evita
                # borrar un labelmap del usuario que se llame "Mascara_algo".
                is_module_mask = (
                    nm == _NAME_COMBINED
                    or nm.startswith(_NAME_COMBINED + "_")
                    or nm == _NAME_COMBINED_GRAY
                    or nm.startswith(_NAME_COMBINED_GRAY + "_")
                )
                if not is_module_mask:
                    continue
                try:
                    arr = slicer.util.arrayFromVolume(node)
                    if arr is None:
                        continue
                    vals = np.unique(np.asarray(arr))
                    if vals.size and int(vals.max()) >= 256:
                        removed.append(f"{nm} (max={int(vals.max())})")
                        slicer.mrmlScene.RemoveNode(node)
                except Exception as e:
                    _log(f"_removeLegacyMaskNodes: no se pudo evaluar '{nm}': {e}")
        except Exception as e:
            _log(f"_removeLegacyMaskNodes: error general: {e}")
        if removed:
            _log(
                f"_removeLegacyMaskNodes: eliminadas {len(removed)} mascaras con niveles >=256: {removed}"
            )
        else:
            _log("_removeLegacyMaskNodes: ninguna mascara con niveles >=256")
        return removed

    def onGenerateVoxel(self):
        """Genera voxel.dat a partir de la mascara gris y la tabla de mapeo."""
        maskNode = self.voxelMaskSelector.currentNode()
        if not maskNode:
            slicer.util.errorDisplay(
                "Seleccione la mascara gris (Mascaras_Combinadas_Gris)."
            )
            return
        if self.voxelMappingTable.rowCount == 0:
            self.onLoadVoxelLevels()
        if self.voxelMappingTable.rowCount == 0:
            slicer.util.errorDisplay(
                "La tabla de universos esta vacia: la mascara seleccionada no tiene "
                "niveles utiles. Verifique la 'Mascara gris' y regenere las mascaras."
            )
            return

        universo_by_gray = {}
        try:
            for row in range(self.voxelMappingTable.rowCount):
                gl_item = self.voxelMappingTable.item(row, 0)
                uni_item = self.voxelMappingTable.item(row, 1)
                if gl_item is None or uni_item is None:
                    continue
                gl = int(str(gl_item.text()).strip())
                uni = int(str(uni_item.text()).strip())
                universo_by_gray[gl] = uni
        except Exception:
            legacy = []
            for row in range(self.voxelMappingTable.rowCount):
                gl_item = self.voxelMappingTable.item(row, 0)
                uni_item = self.voxelMappingTable.item(row, 1)
                if gl_item is None or uni_item is None:
                    continue
                txt = str(uni_item.text()).strip()
                try:
                    int(txt)
                except Exception:
                    try:
                        gl = int(str(gl_item.text()).strip())
                    except Exception:
                        gl = -1
                    legacy.append(gl)
            legacy.sort()
            msg = (
                "Niveles sin universo asignado: "
                + ", ".join(str(g) for g in legacy)
                + ".\n"
            )
            if any(g >= 256 for g in legacy):
                msg += (
                    "Los niveles >= 256 no existen en organ_config.jsonc: la mascara gris "
                    "seleccionada es de una version anterior.\n"
                    "Presione 'Generar mascaras' para regenerarla con los gray_level de la config."
                )
            else:
                msg += "Asigne un universo (U=) valido en la tabla o regenere las mascaras."
            _log("onGenerateVoxel: " + msg.replace("\n", " "))
            slicer.util.errorDisplay(msg)
            return

        orientation = str(self.voxelOrientationCombo.currentText)
        outputPath = str(self.voxelOutputPath.text).strip()
        if not outputPath:
            slicer.util.errorDisplay("Indique la ruta de salida del voxel.dat.")
            return

        resampleX = float(self.voxelResampleX.value)
        resampleY = float(self.voxelResampleY.value)
        resampleZ = float(self.voxelResampleZ.value)

        logic = BNCTLogic()
        voxelNode = maskNode
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            if (
                abs(resampleX - 1.0) > 1e-9
                or abs(resampleY - 1.0) > 1e-9
                or abs(resampleZ - 1.0) > 1e-9
            ):
                voxelNode = logic.resampleMaskVolume(
                    maskNode, resampleX, resampleY, resampleZ
                )
            if voxelNode is None:
                result = None
            else:
                result = logic.generateVoxelDat(
                    voxelNode, universo_by_gray, orientation, outputPath
                )
        finally:
            qt.QApplication.restoreOverrideCursor()

        if result is None:
            self.voxelStatusLabel.setText(
                "No se pudo generar el voxel.dat. Ver _bnct_debug.txt."
            )
            self.voxelStatusLabel.setStyleSheet("color: #800; font-size: 11px;")
            self.voxelReportText.setPlainText(
                "Error al generar el voxel.dat. Ver el log."
            )
            slicer.util.errorDisplay("No se pudo generar el voxel.dat. Ver el log.")
            return

        shape = result.get("shape", ())
        universos = result.get("universos", [])
        self.voxelStatusLabel.setText(f"voxel.dat generado: {result.get('path')}")
        self.voxelStatusLabel.setStyleSheet("color: #060; font-size: 11px;")
        resample_txt = "sin cambios (1,1,1)"
        if voxelNode is not maskNode:
            resample_txt = f"{resampleX:g} x {resampleY:g} x {resampleZ:g}"
        self.voxelReportText.setPlainText(
            f"Mascara origen: {maskNode.GetName()}\n"
            f"Mascara usada: {voxelNode.GetName()}\n"
            f"Resample (X,Y,Z): {resample_txt}\n"
            f"Orientacion: {result.get('orientation')}\n"
            f"Shape (Y x X x Z): {shape}\n"
            f"Voxeles: {result.get('voxels')}\n"
            f"Universos usados (U=): {', '.join(str(u) for u in universos)}\n"
            f"Valores RLE escritos: {result.get('values')}\n"
            f"Archivo: {result.get('path')}"
        )
        slicer.util.infoDisplay(f"voxel.dat generado en:\n{result.get('path')}")

    def onBrowseVoxelOutput(self):
        """Permite elegir la carpeta de salida del voxel.dat."""
        current = os.path.dirname(str(self.voxelOutputPath.text).strip())
        if not current or not os.path.isdir(current):
            current = _MODULE_DIR
        folder = qt.QFileDialog.getExistingDirectory(
            slicer.util.mainWindow(), "Carpeta de salida del voxel.dat", current
        )
        if folder:
            self.voxelOutputPath.setText(os.path.join(folder, "voxel.dat"))

    def _voxelUniverseDefaults(self):
        """Devuelve {gray_level: universo} leido de organ_config.jsonc (incluye AIRE)."""
        result = {}
        try:
            cfg = BNCTLogic.load_jsonc(_ORGAN_CONFIG_PATH)
            aire = cfg.get("aire") or {}
            if "gray_level" in aire and "universo" in aire:
                result[int(aire["gray_level"])] = int(aire["universo"])
            for roi in cfg.get("rois", []):
                if "gray_level" in roi and "universo" in roi:
                    result[int(roi["gray_level"])] = int(roi["universo"])
        except Exception as e:
            _log(f"_voxelUniverseDefaults: no se pudo leer organ_config.jsonc: {e}")
        if 0 not in result:
            result[0] = 1  # AIRE por defecto
        return result

    # ------------------------------------------------------------------
    # UI: Dosimetry
    # ------------------------------------------------------------------
    def _setupDosimetryUI(self):
        fl = qt.QFormLayout()  # sin padre: el tab ya tiene su QVBoxLayout
        self.dosimetryLayout.addLayout(fl)

        # Carpeta PHITS
        self.phitsDirSelector = ctk.ctkPathLineEdit()
        self.phitsDirSelector.filters = ctk.ctkPathLineEdit.Dirs
        fl.addRow("Carpeta PHITS Output:", self.phitsDirSelector)

        # Imagen de referencia
        self.referenceSelector = slicer.qMRMLNodeComboBox()
        self.referenceSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.referenceSelector.selectNodeUponCreation = True
        self.referenceSelector.addEnabled = False
        self.referenceSelector.removeEnabled = False
        self.referenceSelector.noneEnabled = True
        self.referenceSelector.showHidden = False
        self.referenceSelector.setMRMLScene(slicer.mrmlScene)
        self.referenceSelector.connect(
            "currentNodeChanged(vtkMRMLNode*)", self._onReferenceNodeChanged
        )
        fl.addRow("Imagen de Referencia (CT/MRI):", self.referenceSelector)

        # Conc. Boro
        self.boronConcentration = qt.QDoubleSpinBox()
        self.boronConcentration.setRange(0, 1000)
        self.boronConcentration.value = 15.0
        fl.addRow("Conc. Boro en Sangre (ppm):", self.boronConcentration)

        # Factor TIR
        self.factorTIR = qt.QDoubleSpinBox()
        self.factorTIR.value = 1.0
        fl.addRow("Factor TIR:", self.factorTIR)

        # Factor de escala
        self.scaleFactor = qt.QDoubleSpinBox()
        self.scaleFactor.setDecimals(15)
        self.scaleFactor.setRange(0, 1e20)
        self.scaleFactor.value = 1.0
        fl.addRow("Factor de Escala (S):", self.scaleFactor)

        # Offsets
        self.offsetX = qt.QDoubleSpinBox()
        self.offsetX.setRange(-5000, 5000)
        self.offsetX.value = 0.0
        self.offsetY = qt.QDoubleSpinBox()
        self.offsetY.setRange(-5000, 5000)
        self.offsetY.value = 0.0
        self.offsetZ = qt.QDoubleSpinBox()
        self.offsetZ.setRange(-5000, 5000)
        self.offsetZ.value = 0.0
        fl.addRow("Offset X (mm):", self.offsetX)
        fl.addRow("Offset Y (mm):", self.offsetY)
        fl.addRow("Offset Z (mm):", self.offsetZ)

        # Campos
        self.fieldCountCombo = qt.QComboBox()
        self.fieldCountCombo.addItems(["1 Campo", "2 Campos"])
        fl.addRow("Campos:", self.fieldCountCombo)
        self.peso1Spin = qt.QDoubleSpinBox()
        self.peso1Spin.setRange(0, 100)
        self.peso1Spin.setDecimals(4)
        self.peso1Spin.value = 1.0
        fl.addRow("Peso Campo 1:", self.peso1Spin)
        self.peso2Spin = qt.QDoubleSpinBox()
        self.peso2Spin.setRange(0, 100)
        self.peso2Spin.setDecimals(4)
        self.peso2Spin.value = 1.0
        self.peso2Spin.visible = False
        fl.addRow("Peso Campo 2:", self.peso2Spin)
        self.fieldCountCombo.connect(
            "currentIndexChanged(int)", self._onFieldCountChanged
        )

        # Botones de accion
        self.applyButton = qt.QPushButton("Interpolar Dosis")
        self.applyButton.toolTip = (
            "Alinea centros, carga archivos y genera volumen de dosis."
        )
        self.applyButton.connect("clicked(bool)", self.onApplyButton)
        fl.addRow(self.applyButton)

        self.physDoseButton = qt.QPushButton("Calcular Dosis Fisica")
        self.physDoseButton.toolTip = (
            "Componentes KERMA -> dosis fisica con Tissue/Blood Ratio."
        )
        self.physDoseButton.connect("clicked(bool)", self.onPhysDose)
        fl.addRow(self.physDoseButton)

        self.rbeDoseButton = qt.QPushButton("Calcular Dosis Pesada por RBE")
        self.rbeDoseButton.toolTip = "Componentes, RBE, dosis pesada."
        self.rbeDoseButton.connect("clicked(bool)", self.onRBEDose)
        fl.addRow(self.rbeDoseButton)

        self.isoDoseButton = qt.QPushButton("Calcular Dosis Isoefectiva")
        self.isoDoseButton.toolTip = (
            "Dosis isoefectiva con modelo radiobiologico models.jsonc."
        )
        self.isoDoseButton.connect("clicked(bool)", self.onIsoDose)
        fl.addRow(self.isoDoseButton)

        self.btnCenterBody = qt.QPushButton("Centrar BODY en origen")
        self.btnCenterBody.toolTip = (
            "Aplica un transform lineal que lleva el centro del BODY al origen (0,0,0) RAS."
        )
        self.btnCenterBody.connect("clicked(bool)", self.onCenterBody)
        fl.addRow(self.btnCenterBody)

        self.btnUncenterBody = qt.QPushButton("Quitar centrado BODY")
        self.btnUncenterBody.toolTip = (
            "Elimina el transform de centrado y restaura la posicion original."
        )
        self.btnUncenterBody.connect("clicked(bool)", self.onUncenterBody)
        fl.addRow(self.btnUncenterBody)

        # --- Seccion KERMA por ROI ---
        configBtn = ctk.ctkCollapsibleButton()
        configBtn.text = "Configuracion KERMA por ROI"
        configBtn.collapsed = False
        self.dosimetryLayout.addWidget(configBtn)
        cfl = qt.QFormLayout(configBtn)

        btnLay = qt.QHBoxLayout()
        self.btnSaveMat = qt.QPushButton("Guardar MAT")
        self.btnSaveJson = qt.QPushButton("Guardar JSON")
        self.btnReload = qt.QPushButton("Recargar Default")
        self.btnAddRow = qt.QPushButton("+ Fila")
        self.btnDelRow = qt.QPushButton("- Fila")
        for w in (
            self.btnSaveMat,
            self.btnSaveJson,
            self.btnReload,
            self.btnAddRow,
            self.btnDelRow,
        ):
            btnLay.addWidget(w)
        btnLay.addStretch()
        cfl.addRow(btnLay)

        self.configTable = qt.QTableWidget()
        self.configTable.setColumnCount(10)
        self.configTable.setHorizontalHeaderLabels(
            [
                "nroi",
                "nmesh",
                "Nombre",
                "RBE_B",
                "RBE_F",
                "RBE_T",
                "RBE_G",
                "Ratio",
                "Tumor",
                "GrayLvl",
            ]
        )
        self.configTable.horizontalHeader().setStretchLastSection(True)
        self.configTable.setAlternatingRowColors(True)
        self.configTable.setMinimumHeight(250)
        cfl.addRow(self.configTable)

        self.btnSaveMat.connect("clicked(bool)", self.onSaveMat)
        self.btnSaveJson.connect("clicked(bool)", self.onSaveJson)
        self.btnReload.connect("clicked(bool)", self.onReloadDefault)
        self.btnAddRow.connect("clicked(bool)", self.onAddRow)
        self.btnDelRow.connect("clicked(bool)", self.onRemoveRow)
        self.configTable.cellChanged.connect(self._onConfigCellChanged)

        self._organ_config = self._loadOrganConfig()
        if self._organ_config:
            self._populateConfigTable(self._organ_config)

        # Cargar modelos radiobiologicos
        self._models_data = None
        try:
            logic = BNCTLogic()
            self._models_data = logic.loadModels()
        except Exception as e:
            logging.warning(f"No se pudieron cargar modelos: {e}")

    # ------------------------------------------------------------------
    # UI: DVH Analysis
    # ------------------------------------------------------------------
    def _setupDVHUI(self):
        fl = qt.QFormLayout()  # sin padre: el tab ya tiene su QVBoxLayout
        self.dvhLayout.addLayout(fl)

        # Limpiar nodos DVH viejos al iniciar
        try:
            for cls, prefix in [
                ("vtkMRMLPlotSeriesNode", "DVH_"),
                ("vtkMRMLTableNode", "Table_DVH_"),
            ]:
                coll = slicer.mrmlScene.GetNodesByClass(cls)
                for i in range(coll.GetNumberOfItems()):
                    node = coll.GetItemAsObject(i)
                    if node and node.GetName().startswith(prefix):
                        slicer.mrmlScene.RemoveNode(node)
            old = slicer.mrmlScene.GetFirstNodeByName("Chart_DVH")
            if old:
                slicer.mrmlScene.RemoveNode(old)
        except:
            pass

        dvhColl = ctk.ctkCollapsibleButton()
        dvhColl.text = "Analisis de DVH"
        dvhColl.collapsed = True
        fl.addRow(dvhColl)
        dvhFl = qt.QFormLayout(dvhColl)

        organW = qt.QWidget()
        organWl = qt.QVBoxLayout(organW)
        organWl.setContentsMargins(0, 0, 0, 0)
        self.dvhTodosCheck = qt.QCheckBox("Todos")
        self.dvhTodosCheck.checked = True
        organWl.addWidget(self.dvhTodosCheck)
        scroll = qt.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(180)
        sw = qt.QWidget()
        self.dvhScrollLayout = qt.QVBoxLayout(sw)
        self.dvhScrollLayout.setSpacing(2)
        scroll.setWidget(sw)
        organWl.addWidget(scroll)
        dvhFl.addRow("Organos:", organW)
        self.dvhTodosCheck.connect("toggled(bool)", self._onDvhTodosToggled)

        self.dvhDoseTypeCombo = qt.QComboBox()
        self.dvhDoseTypeCombo.addItems(
            ["Dosis Fisica", "Dosis RBE", "Dosis Isoefectiva"]
        )
        dvhFl.addRow("Tipo de Dosis:", self.dvhDoseTypeCombo)

        self.dvhModeloCombo = qt.QComboBox()
        self.dvhModeloCombo.setMinimumWidth(300)
        self.dvhModeloCombo.visible = False
        dvhFl.addRow("Modelo:", self.dvhModeloCombo)

        self.dvhB10Spin = qt.QDoubleSpinBox()
        self.dvhB10Spin.setRange(0, 1000)
        self.dvhB10Spin.value = 15.0
        dvhFl.addRow("Conc. Boro (ppm):", self.dvhB10Spin)

        self.dvhTirSpin = qt.QDoubleSpinBox()
        self.dvhTirSpin.setRange(0, 1000)
        self.dvhTirSpin.value = 1.0
        dvhFl.addRow("TIR (min):", self.dvhTirSpin)

        self.maskVolumeSelector = slicer.qMRMLNodeComboBox()
        self.maskVolumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.maskVolumeSelector.selectNodeUponCreation = True
        self.maskVolumeSelector.addEnabled = False
        self.maskVolumeSelector.removeEnabled = False
        self.maskVolumeSelector.noneEnabled = False
        self.maskVolumeSelector.showHidden = False
        self.maskVolumeSelector.setMRMLScene(slicer.mrmlScene)
        dvhFl.addRow("Mascara ROIs:", self.maskVolumeSelector)

        self.dvhButton = qt.QPushButton("Calcular DVH")
        self.dvhButton.enabled = True
        dvhFl.addRow(self.dvhButton)
        self.dvhButton.connect("clicked(bool)", self.onDVHButton)
        self.dvhDoseTypeCombo.connect(
            "currentIndexChanged(int)", self._onDvhDoseTypeChanged
        )

        self._dvh_models = self._loadModels()
        self.dvhModeloCombo.clear()
        for m in self._dvh_models:
            self.dvhModeloCombo.addItem(f"{m['id']}. {m['name']}")

        self._dvh_config = self._loadConfigForDvh()
        self._setupUTCPUI()

    # ------------------------------------------------------------------
    # UI: Mascaras (labelmaps por ROI)
    # ------------------------------------------------------------------
    def _setupMasksUI(self):
        """UI de la pestana Mascaras: genera labelmaps por ROI desde una segmentacion o un volumen."""
        fl = qt.QFormLayout()  # sin padre: el tab ya tiene su QVBoxLayout
        self.masksLayout.addLayout(fl)

        # Selector de estructura: acepta una Segmentacion (recomendado) o un volumen escalar
        self.structureSelector = slicer.qMRMLNodeComboBox()
        self.structureSelector.nodeTypes = [
            "vtkMRMLSegmentationNode",
            "vtkMRMLScalarVolumeNode",
        ]
        self.structureSelector.selectNodeUponCreation = True
        self.structureSelector.addEnabled = False
        self.structureSelector.removeEnabled = False
        self.structureSelector.noneEnabled = False
        self.structureSelector.showHidden = False
        self.structureSelector.setMRMLScene(slicer.mrmlScene)
        self.structureSelector.setToolTip(
            "Nodo de Segmentacion del Segment Editor: cada segmento se empareja "
            "con su ROI por NOMBRE. Tambien acepta un volumen escalar donde cada ROI "
            "son los voxeles cuyo valor coincide con su gray_level."
        )
        fl.addRow(
            "Estructura (Segmentacion o Volumen con gray_level):",
            self.structureSelector,
        )

        # Volumen de referencia: define la geometria (origin/spacing) de los labelmaps de salida
        self.segReferenceSelector = slicer.qMRMLNodeComboBox()
        self.segReferenceSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.segReferenceSelector.selectNodeUponCreation = True
        self.segReferenceSelector.addEnabled = False
        self.segReferenceSelector.removeEnabled = False
        self.segReferenceSelector.noneEnabled = True
        self.segReferenceSelector.showHidden = False
        self.segReferenceSelector.setMRMLScene(slicer.mrmlScene)
        self.segReferenceSelector.setToolTip(
            "Volumen de referencia que define la geometria (origin/spacing) de los labelmaps de salida."
        )
        fl.addRow(
            "Volumen de Referencia (geometria salida):", self.segReferenceSelector
        )

        # Info de la config
        self.masksConfigLabel = qt.QLabel(self._masksConfigInfo())
        self.masksConfigLabel.setWordWrap(True)
        self.masksConfigLabel.setStyleSheet("color: #666; font-size: 11px;")
        fl.addRow("Config ROIs:", self.masksConfigLabel)

        # Opciones
        self.masksOptionsCollapsible = ctk.ctkCollapsibleButton()
        self.masksOptionsCollapsible.text = "Opciones de generacion"
        self.masksOptionsCollapsible.collapsed = False
        self.masksLayout.addWidget(self.masksOptionsCollapsible)
        optLay = qt.QVBoxLayout(self.masksOptionsCollapsible)  # padre = QWidget -> OK
        self.chkMasksIndividual = qt.QCheckBox("Generar un labelmap binario por ROI")
        self.chkMasksIndividual.checked = True
        self.chkMasksCombined = qt.QCheckBox(
            "Generar labelmap combinado (cada ROI = su gray_level)"
        )
        self.chkMasksCombined.checked = True
        self.chkMasksSkipEmpty = qt.QCheckBox("Omitir ROIs sin voxeles")
        self.chkMasksSkipEmpty.checked = True
        optLay.addWidget(self.chkMasksIndividual)
        optLay.addWidget(self.chkMasksCombined)
        optLay.addWidget(self.chkMasksSkipEmpty)

        # Botones
        btnRow = qt.QHBoxLayout()
        self.btnGenerateMasks = qt.QPushButton("Generar Mascaras")
        self.btnGenerateMasks.toolTip = (
            "Genera los labelmaps y verifica que no haya superposicion."
        )
        self.btnGenerateMasks.connect("clicked(bool)", self.onGenerateMasks)
        self.btnVerifyOverlap = qt.QPushButton("Verificar Superposicion")
        self.btnVerifyOverlap.toolTip = "Solo verifica, sin crear volumenes."
        self.btnVerifyOverlap.connect("clicked(bool)", self.onVerifyOverlap)
        self.btnCleanupMasks = qt.QPushButton("Limpiar Mascaras")
        self.btnCleanupMasks.toolTip = "Borra los labelmaps generados por este modulo."
        self.btnCleanupMasks.connect("clicked(bool)", self.onCleanupMasks)
        btnRow.addWidget(self.btnGenerateMasks)
        btnRow.addWidget(self.btnVerifyOverlap)
        btnRow.addWidget(self.btnCleanupMasks)
        btnRow.addStretch(1)
        self.masksLayout.addLayout(btnRow)

        # Tabla de resultados
        self.masksTable = qt.QTableWidget()
        self.masksTable.setColumnCount(5)
        self.masksTable.setHorizontalHeaderLabels(
            ["ROI", "gray_level", "nmesh", "Voxeles", "Estado"]
        )
        self.masksTable.horizontalHeader().setStretchLastSection(True)
        self.masksTable.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        self.masksTable.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.masksTable.setMinimumHeight(260)
        self.masksLayout.addWidget(self.masksTable)

        # Informe de superposicion
        self.masksReportText = qt.QTextEdit()
        self.masksReportText.readOnly = True
        self.masksReportText.setMinimumHeight(140)
        self.masksReportText.plainText = (
            "Selecciona el volumen de estructura y pulsa 'Generar Mascaras'."
        )
        self.masksLayout.addWidget(self.masksReportText)

    def _setupCollisionsUI(self):
        fl = qt.QFormLayout()
        self.collisionsLayout.addLayout(fl)

        self.patientModelSelector = slicer.qMRMLNodeComboBox()
        self.patientModelSelector.nodeTypes = ["vtkMRMLModelNode"]
        self.patientModelSelector.selectNodeUponCreation = True
        self.patientModelSelector.addEnabled = False
        self.patientModelSelector.removeEnabled = False
        self.patientModelSelector.noneEnabled = True
        self.patientModelSelector.showHidden = False
        self.patientModelSelector.setMRMLScene(slicer.mrmlScene)
        self.patientModelSelector.setToolTip("Modelo 3D del paciente (superficie).")
        fl.addRow("Modelo paciente:", self.patientModelSelector)

        self.coneModelSelector = slicer.qMRMLNodeComboBox()
        self.coneModelSelector.nodeTypes = ["vtkMRMLModelNode"]
        self.coneModelSelector.selectNodeUponCreation = True
        self.coneModelSelector.addEnabled = False
        self.coneModelSelector.removeEnabled = False
        self.coneModelSelector.noneEnabled = True
        self.coneModelSelector.showHidden = False
        self.coneModelSelector.setMRMLScene(slicer.mrmlScene)
        self.coneModelSelector.setToolTip(
            "Modelo 3D del cono/pared del reactor (Conos_Pared_Reactor)."
        )
        fl.addRow("Modelo cono/pared:", self.coneModelSelector)

        self.btnDetectCollisions = qt.QPushButton("Detectar Colisiones")
        self.btnDetectCollisions.setToolTip(
            "Ejecuta vtkCollisionDetectionFilter entre los dos modelos seleccionados."
        )
        self.btnDetectCollisions.connect("clicked(bool)", self.onDetectCollisions)
        fl.addRow(self.btnDetectCollisions)

        self.collisionStatusLabel = qt.QLabel(
            "Seleccione los dos modelos y presione Detectar Colisiones."
        )
        self.collisionStatusLabel.setWordWrap(True)
        self.collisionStatusLabel.setStyleSheet("color: #666; font-size: 11px;")
        fl.addRow(self.collisionStatusLabel)

        self.collisionReportText = qt.QTextEdit()
        self.collisionReportText.setReadOnly(True)
        self.collisionReportText.setMaximumHeight(140)
        self.collisionReportText.setPlaceholderText(
            "Resultado de la deteccion de colisiones..."
        )
        fl.addRow(self.collisionReportText)

        patient = slicer.mrmlScene.GetFirstNodeByName("Paciente")
        if patient:
            self.patientModelSelector.setCurrentNode(patient)
        cone = slicer.mrmlScene.GetFirstNodeByName("Cono")
        if not cone:
            cone = slicer.mrmlScene.GetFirstNodeByName("Conos_Pared_Reactor")
        if cone:
            self.coneModelSelector.setCurrentNode(cone)

    def _setupVoxelUI(self):
        fl = qt.QFormLayout()
        self.voxelLayout.addLayout(fl)

        self.voxelMaskSelector = slicer.qMRMLNodeComboBox()
        self.voxelMaskSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.voxelMaskSelector.selectNodeUponCreation = True
        self.voxelMaskSelector.addEnabled = False
        self.voxelMaskSelector.removeEnabled = False
        self.voxelMaskSelector.noneEnabled = True
        self.voxelMaskSelector.showHidden = False
        self.voxelMaskSelector.setMRMLScene(slicer.mrmlScene)
        self.voxelMaskSelector.setToolTip(
            "Labelmap gris combinado (Mascaras_Combinadas_Gris). Sus valores son los gray_level "
            "que se traducen a universo (U=) en el voxel.dat."
        )
        fl.addRow("Mascara gris:", self.voxelMaskSelector)

        # Purgar mascaras creadas por la version anterior (niveles 256+).
        removed_legacy = self._removeLegacyMaskNodes()
        if removed_legacy:
            _log(f"_setupVoxelUI: purgadas mascaras viejas: {removed_legacy}")

        grayNode = slicer.mrmlScene.GetFirstNodeByName(_NAME_COMBINED_GRAY)
        for i in range(
            slicer.mrmlScene.GetNumberOfNodesByClass("vtkMRMLLabelMapVolumeNode")
        ):
            n = slicer.mrmlScene.GetNthNodeByClass(i, "vtkMRMLLabelMapVolumeNode")
            nm = n.GetName() or ""
            if nm == _NAME_COMBINED_GRAY or nm.startswith(_NAME_COMBINED_GRAY + "_"):
                grayNode = n  # el ultimo encontrado = el mas nuevo
        if grayNode:
            self.voxelMaskSelector.setCurrentNode(grayNode)

        self.voxelLoadLevelsButton = qt.QPushButton("Cargar gray_levels")
        self.voxelLoadLevelsButton.setToolTip(
            "Lee los valores unicos de la mascara y llena la tabla de mapeo."
        )
        self.voxelLoadLevelsButton.connect("clicked(bool)", self.onLoadVoxelLevels)
        fl.addRow(self.voxelLoadLevelsButton)

        self.voxelMappingTable = qt.QTableWidget()
        self.voxelMappingTable.setColumnCount(2)
        self.voxelMappingTable.setHorizontalHeaderLabels(
            ["gray_level", "universo (U=)"]
        )
        self.voxelMappingTable.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.voxelMappingTable.setMinimumHeight(180)
        self.voxelMappingTable.horizontalHeader().setStretchLastSection(True)
        fl.addRow("Mapeo de universos:", self.voxelMappingTable)

        self.voxelOrientationCombo = qt.QComboBox()
        self.voxelOrientationCombo.addItems(["Axial", "Sagital", "Coronal"])
        self.voxelOrientationCombo.setToolTip(
            "Orientacion de los cortes, igual que en genera_voxel_PHITS_v9.m."
        )
        fl.addRow("Orientacion:", self.voxelOrientationCombo)

        # --- Resample opcional (0 < f <= 1; 1 = sin cambios) ---
        resampleLayout = qt.QHBoxLayout()
        self.voxelResampleX = qt.QDoubleSpinBox()
        self.voxelResampleY = qt.QDoubleSpinBox()
        self.voxelResampleZ = qt.QDoubleSpinBox()
        for _box, _axis in (
            (self.voxelResampleX, "X"),
            (self.voxelResampleY, "Y"),
            (self.voxelResampleZ, "Z"),
        ):
            _box.setDecimals(2)
            _box.setMinimum(0.01)
            _box.setMaximum(1.00)
            _box.setSingleStep(0.05)
            _box.setValue(1.00)
            _box.setMaximumWidth(70)
            _box.setToolTip(
                f"Factor de resample en {_axis}. 1,00 = sin cambios; "
                f"0,50 = mitad de resolucion en ese eje."
            )
            resampleLayout.addWidget(qt.QLabel(_axis + ":"))
            resampleLayout.addWidget(_box)
        resampleLayout.addStretch(1)
        fl.addRow("Resample (0 < f <= 1):", resampleLayout)

        pathLayout = qt.QHBoxLayout()
        self.voxelOutputPath = qt.QLineEdit(
            os.path.join(_MODULE_DIR, "voxel", "voxel.dat")
        )
        self.voxelOutputPath.setToolTip("Ruta del archivo voxel.dat a generar.")
        pathLayout.addWidget(self.voxelOutputPath)
        self.voxelBrowseButton = qt.QPushButton("...")
        self.voxelBrowseButton.setMaximumWidth(30)
        self.voxelBrowseButton.connect("clicked(bool)", self.onBrowseVoxelOutput)
        pathLayout.addWidget(self.voxelBrowseButton)
        fl.addRow("Salida:", pathLayout)

        self.voxelGenerateButton = qt.QPushButton("Generar voxel.dat")
        self.voxelGenerateButton.setToolTip(
            "Aplica el mapeo, orienta, comprime con RLE y escribe el voxel.dat para PHITS."
        )
        self.voxelGenerateButton.connect("clicked(bool)", self.onGenerateVoxel)
        fl.addRow(self.voxelGenerateButton)

        self.voxelStatusLabel = qt.QLabel(
            "Seleccione la mascara gris, cargue los gray_levels y presione Generar voxel.dat."
        )
        self.voxelStatusLabel.setWordWrap(True)
        self.voxelStatusLabel.setStyleSheet("color: #666; font-size: 11px;")
        fl.addRow(self.voxelStatusLabel)

        self.voxelReportText = qt.QTextEdit()
        self.voxelReportText.setReadOnly(True)
        self.voxelReportText.setMaximumHeight(120)
        self.voxelReportText.setPlaceholderText(
            "Resultado de la generacion del voxel.dat..."
        )
        fl.addRow(self.voxelReportText)

        # Autocarga: la config de universos (json) se aplica sola al elegir
        # o cambiar la mascara gris, sin presionar 'Cargar gray_levels'.
        self.voxelMaskSelector.connect(
            "currentNodeChanged(vtkMRMLNode*)", self._onVoxelMaskChanged
        )
        try:
            if getattr(self, "_organ_config", None) is None:
                self._organ_config = self._loadOrganConfig()
            if self.voxelMaskSelector.currentNode() is not None:
                self.onLoadVoxelLevels(silent=True)
            else:
                self._onVoxelMaskChanged(None)
        except Exception as e:
            _log(f"_setupVoxelUI: autocarga de niveles fallo: {e}")

    def _masksConfigInfo(self):
        """Devuelve una descripcion corta de la config de ROIs disponible."""
        n = len(getattr(self, "_organ_config", []) or [])
        if n:
            return f"{n} ROIs desde config cargada ({os.path.basename(_ORGAN_CONFIG_PATH)})"
        return f"Se usara {os.path.basename(_MAT_PATH)}"

    # ------------------------------------------------------------------
    # Handlers compartidos / Dosimetry
    # ------------------------------------------------------------------
    def _computeCenterOffset(self):
        # GetBounds devuelve los bounds LOCALES del volumen (sin transforms):
        # es coherente con la interpolacion, que tambien usa la matriz
        # IJK->RAS local de la CT.
        refNode = self.referenceSelector.currentNode()
        if not refNode:
            return None
        phitsDir = self.phitsDirSelector.currentPath
        if not phitsDir or not os.path.isdir(phitsDir):
            return None
        n_files, _ = self._discoverPhitsFiles(phitsDir)
        if not n_files:
            return None
        import numpy as np

        bounds = [0] * 6
        refNode.GetBounds(bounds)
        center_ras = np.array(
            [
                (bounds[0] + bounds[1]) / 2.0,
                (bounds[2] + bounds[3]) / 2.0,
                (bounds[4] + bounds[5]) / 2.0,
            ]
        )
        logic = BNCTLogic()
        phits = logic.parsePHITS(n_files[0])
        origin = np.array(phits["origin"])
        shape = np.array(phits["shape"])
        spacing = np.array(phits["spacing"])
        center_phits = origin + (shape[::-1] - 1) * spacing / 2.0
        offset = center_ras - center_phits
        return (offset[0], offset[1], offset[2])

    def _onReferenceNodeChanged(self, node):
        if not node:
            return
        offset = self._computeCenterOffset()
        if offset is not None:
            self.offsetX.value = offset[0]
            self.offsetY.value = offset[1]
            self.offsetZ.value = offset[2]
            logging.info(f"Auto-centrado: offset={offset}")

    def _onFieldCountChanged(self, idx):
        self.peso2Spin.visible = idx == 1

    def _discoverPhitsFiles(self, directory):
        neutron_files = []
        gamma_path = None
        out_files = sorted(
            [
                f
                for f in os.listdir(directory)
                if f.endswith(".out")
                and "_err" not in f.lower()
                and "_error" not in f.lower()
            ]
        )
        for f in out_files:
            low = f.lower()
            if "photon" in low or "gamma" in low:
                gamma_path = os.path.join(directory, f)
            elif re.match(r"mesh-\d+\.out", low):
                neutron_files.append(os.path.join(directory, f))

        def mesh_number(p):
            m = re.search(r"mesh-(\d+)", os.path.basename(p).lower())
            return int(m.group(1)) if m else 9999

        neutron_files.sort(key=mesh_number)
        return neutron_files, gamma_path

    def _defaultMatPath(self):
        return _MAT_PATH

    def _populateConfigTable(self, config):
        self.configTable.blockSignals(True)
        self.configTable.setRowCount(len(config))
        for i, row in enumerate(config):
            self.configTable.setItem(i, 0, qt.QTableWidgetItem(str(int(row["nroi"]))))
            self.configTable.setItem(i, 1, qt.QTableWidgetItem(str(int(row["nmesh"]))))
            self.configTable.setItem(i, 2, qt.QTableWidgetItem(str(row["name"])))
            self.configTable.setItem(
                i, 3, qt.QTableWidgetItem(f"{row['rbe_boron']:.2f}")
            )
            self.configTable.setItem(
                i, 4, qt.QTableWidgetItem(f"{row['rbe_fast']:.2f}")
            )
            self.configTable.setItem(
                i, 5, qt.QTableWidgetItem(f"{row['rbe_thermal']:.2f}")
            )
            self.configTable.setItem(
                i, 6, qt.QTableWidgetItem(f"{row['rbe_gamma']:.2f}")
            )
            self.configTable.setItem(i, 7, qt.QTableWidgetItem(f"{row['ratio']:.1f}"))
            chk = qt.QTableWidgetItem()
            chk.setCheckState(
                qt.Qt.Checked if row.get("tumor", False) else qt.Qt.Unchecked
            )
            self.configTable.setItem(i, 8, chk)
            self.configTable.setItem(
                i, 9, qt.QTableWidgetItem(str(int(row.get("gray_level", 0))))
            )
        self.configTable.resizeColumnsToContents()
        self.configTable.blockSignals(False)

    def _readConfigFromTable(self):
        # Conserva claves extra de la config persistida ('universo', 'tir',
        # 'utcp_*'): la tabla solo edita 10 columnas y reconstruir las filas
        # sin merge BORRARIA esas claves del jsonc en cada edicion de celda.
        prev = {}
        try:
            if os.path.exists(_ORGAN_CONFIG_PATH):
                for r in load_jsonc(_ORGAN_CONFIG_PATH).get("rois", []):
                    if isinstance(r, dict) and r.get("name"):
                        prev[r["name"]] = r
        except Exception:
            pass
        extra_keys = (
            "universo",
            "tir",
            "utcp_isoe_model",
            "utcp_ntcp_model",
            "utcp_include",
        )
        config = []
        for i in range(self.configTable.rowCount):

            def val(col):
                item = self.configTable.item(i, col)
                return item.text() if item else ""

            row = {
                "nroi": int(float(val(0))) if val(0) else 0,
                "nmesh": int(float(val(1))) if val(1) else 0,
                "name": val(2),
                "rbe_boron": float(val(3)) if val(3) else 0.0,
                "rbe_fast": float(val(4)) if val(4) else 0.0,
                "rbe_thermal": float(val(5)) if val(5) else 0.0,
                "rbe_gamma": float(val(6)) if val(6) else 0.0,
                "ratio": float(val(7)) if val(7) else 0.0,
                "tumor": self.configTable.item(i, 8).checkState() == qt.Qt.Checked
                if self.configTable.item(i, 8)
                else False,
                "gray_level": int(float(val(9))) if val(9) else 0,
            }
            for k in extra_keys:
                if k in prev.get(row["name"], {}) and k not in row:
                    row[k] = prev[row["name"]][k]
            config.append(row)
        return config

    def _syncOrganConfigJson(self):
        config = self._readConfigFromTable()
        try:
            data = {}
            if os.path.exists(_ORGAN_CONFIG_PATH):
                try:
                    data = load_jsonc(_ORGAN_CONFIG_PATH)
                except Exception:
                    data = {}
            data["rois"] = config
            d = os.path.dirname(_ORGAN_CONFIG_PATH)
            if not os.path.exists(d):
                os.makedirs(d)
            with open(_ORGAN_CONFIG_PATH, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            logging.warning(f"No se pudo sincronizar organ_config.jsonc: {e}")

    def _loadOrganConfig(self):
        if os.path.exists(_ORGAN_CONFIG_PATH):
            try:
                logic = BNCTLogic()
                data = logic.load_jsonc(_ORGAN_CONFIG_PATH)
                config = data.get("rois", [])
                if config:
                    return config
            except Exception as e:
                logging.warning(f"Fallback a MAT: {e}")
        if os.path.exists(_MAT_PATH):
            try:
                logic = BNCTLogic()
                config = logic.loadConfigFromMat(_MAT_PATH)
                return config
            except Exception as e:
                logging.warning(f"Fallback MAT fallo: {e}")
        return None

    def onSaveMat(self):
        config = self._readConfigFromTable()
        path = self._saveDialogPath("Guardar MAT", "MAT files (*.mat)")
        if not path:
            return
        logic = BNCTLogic()
        try:
            logic.saveConfigToMat(path, config)
            print(f"Config guardada en {path}")
        except Exception as e:
            slicer.util.errorDisplay(f"Error al guardar MAT: {e}")

    def onSaveJson(self):
        config = self._readConfigFromTable()
        path = self._saveDialogPath("Guardar JSON", "JSON files (*.json)")
        if not path:
            return
        logic = BNCTLogic()
        try:
            logic.saveConfigToJson(path, config)
            print(f"Config guardada en {path}")
        except Exception as e:
            slicer.util.errorDisplay(f"Error al guardar JSON: {e}")

    def onReloadDefault(self):
        if os.path.exists(_MAT_PATH):
            logic = BNCTLogic()
            try:
                config = logic.loadConfigFromMat(_MAT_PATH)
                self._populateConfigTable(config)
                self._syncOrganConfigJson()
                print(f"Config recargada: {len(config)} ROIs.")
                return
            except ImportError:
                slicer.util.infoDisplay("h5py no encontrado. Instalando...")
                slicer.util.pip_install("h5py")
                config = logic.loadConfigFromMat(_MAT_PATH)
                self._populateConfigTable(config)
                self._syncOrganConfigJson()
                print(f"Config recargada: {len(config)} ROIs.")
                return
            except Exception as e:
                slicer.util.errorDisplay(f"Error al recargar MAT: {e}")
        else:
            slicer.util.errorDisplay("No se encuentra info_rois.mat en mat/")

    def onAddRow(self):
        table = self.configTable
        row = table.rowCount
        table.blockSignals(True)
        table.insertRow(row)
        nroi = row + 1
        nroi_prev = table.item(row - 1, 0) if row > 0 else None
        if nroi_prev:
            nroi = int(float(nroi_prev.text())) + 1
        table.setItem(row, 0, qt.QTableWidgetItem(str(nroi)))
        table.blockSignals(False)
        self._syncOrganConfigJson()

    def onRemoveRow(self):
        table = self.configTable
        selected = table.currentRow()
        if selected < 0:
            slicer.util.warningDisplay("Seleccione una fila para eliminar.")
            return
        if table.rowCount <= 1:
            slicer.util.warningDisplay("Debe haber al menos una fila.")
            return
        table.blockSignals(True)
        table.removeRow(selected)
        table.blockSignals(False)
        self._syncOrganConfigJson()

    def _onConfigCellChanged(self, row, col):
        self._syncOrganConfigJson()

    def _saveDialogPath(self, title, filter_str):
        result = qt.QFileDialog.getSaveFileName(self.parent, title, "", filter_str)
        if isinstance(result, tuple):
            return result[0]
        return result

    # ------------------------------------------------------------------
    # Acciones Dosimetry
    # ------------------------------------------------------------------
    def onApplyButton(self):
        logic = BNCTLogic()
        phitsDir = self.phitsDirSelector.currentPath
        refNode = self.referenceSelector.currentNode()
        if not phitsDir or not os.path.isdir(phitsDir):
            slicer.util.errorDisplay("Debe seleccionar una carpeta PHITS valida.")
            return
        n_files, g_path = self._discoverPhitsFiles(phitsDir)
        if not n_files:
            slicer.util.errorDisplay(
                "No se encontraron archivos mesh-N.out (neutrones) en la carpeta."
            )
            return
        if not g_path:
            logging.warning("No se encontro archivo gamma. Solo se cargaran neutrones.")
        if refNode:
            offset = self._computeCenterOffset()
            if offset is not None:
                self.offsetX.value = offset[0]
                self.offsetY.value = offset[1]
                self.offsetZ.value = offset[2]
                logging.info(f"Auto-centrado en Interpolar Dosis: offset={offset}")
        print("--- Archivos PHITS encontrados ---")
        print(f"  Neutrones ({len(n_files)}):")
        for f in n_files:
            print(f"    - {os.path.basename(f)}")
        print(f"  Gamma: {os.path.basename(g_path) if g_path else 'No encontrado'}")
        rbe_factors = {
            "offset": (self.offsetX.value, self.offsetY.value, self.offsetZ.value)
        }
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic._cleanupOldVolumes()
            logic.run(n_files, g_path, refNode, rbe_factors, suffix="")
            if self.fieldCountCombo.currentIndex == 1:
                dir_c2 = os.path.join(phitsDir, "campo_2")
                if os.path.isdir(dir_c2):
                    n_files_2, g_path_2 = self._discoverPhitsFiles(dir_c2)
                    if n_files_2:
                        logic.run(
                            n_files_2, g_path_2, refNode, rbe_factors, suffix="_C2"
                        )
                    else:
                        print("  campo_2 vacio o sin archivos validos")
                else:
                    print("  No existe campo_2/ en la carpeta PHITS")
        except Exception as e:
            slicer.util.errorDisplay(f"Error durante el calculo: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def _getKermaVolumesForMesh(self, mesh_number):
        import re

        exclude = [
            "Dosis_Fisica_BNCT",
            "Dosis_RBE_BNCT",
            "Dosis_Final_BNCT",
            "Dosis_IsoE_BNCT",
        ]
        volumes = []
        for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
            name = node.GetName()
            if name in exclude:
                continue
            if name == "Gamma_Photon" or name == "Gamma_Photon_C2":
                volumes.append((node, name))
                continue
            if "Kerma" in name:
                match = re.search(r"Kerma(\d+)(_C2)?$", name)
                if match and int(match.group(1)) == mesh_number:
                    volumes.append((node, name))
        return volumes

    def _buildOrganDialog(self, title, ask_rbe=False):
        config = self._organ_config
        if not config:
            slicer.util.errorDisplay(
                "No hay configuracion de organos disponible. "
                "Configure los ROIs en la pestaña Dosimetry."
            )
            return None
        dialog = qt.QDialog(self.parent)
        dialog.setWindowTitle(title)
        dialog.setMinimumWidth(440)
        layout = qt.QVBoxLayout(dialog)
        layout.setSpacing(6)
        organLay = qt.QHBoxLayout()
        organLabel = qt.QLabel("Organo:")
        organCombo = qt.QComboBox()
        organCombo.setMinimumWidth(200)
        for row in config:
            organCombo.addItem(f"{row['name']}  [mesh {int(row['nmesh'])}]")
        organLay.addWidget(organLabel)
        organLay.addWidget(organCombo, 1)
        layout.addLayout(organLay)
        meshInfo = qt.QLabel()
        meshInfo.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(meshInfo)
        valuesGroup = qt.QGroupBox("Parametros")
        valuesLay = qt.QFormLayout(valuesGroup)
        ratioDisplay = qt.QLabel()
        valuesLay.addRow("Ratio T/S:", ratioDisplay)
        rbeDisplays = {}
        if ask_rbe:
            for key, label in [
                ("rbe_boron", "RBE Boro"),
                ("rbe_fast", "RBE Fast"),
                ("rbe_thermal", "RBE Thermal"),
                ("rbe_gamma", "RBE Gamma"),
            ]:
                lbl = qt.QLabel()
                valuesLay.addRow(f"{label}:", lbl)
                rbeDisplays[key] = lbl
        layout.addWidget(valuesGroup)

        def update_values():
            idx = organCombo.currentIndex
            if idx < 0 or idx >= len(config):
                return
            row = config[idx]
            ratioDisplay.text = f"{row['ratio']:.2f}"
            nmesh = int(row["nmesh"])
            meshInfo.text = f"Mesh: {nmesh}"
            if ask_rbe:
                rbeDisplays["rbe_boron"].text = f"{row['rbe_boron']:.2f}"
                rbeDisplays["rbe_fast"].text = f"{row['rbe_fast']:.2f}"
                rbeDisplays["rbe_thermal"].text = f"{row['rbe_thermal']:.2f}"
                rbeDisplays["rbe_gamma"].text = f"{row['rbe_gamma']:.2f}"

        organCombo.connect("currentIndexChanged(int)", lambda idx: update_values())
        update_values()
        btnLay = qt.QHBoxLayout()
        okBtn = qt.QPushButton("OK")
        cancelBtn = qt.QPushButton("Cancel")
        btnLay.addStretch()
        btnLay.addWidget(okBtn)
        btnLay.addWidget(cancelBtn)
        layout.addLayout(btnLay)
        okBtn.connect("clicked()", dialog.accept)
        cancelBtn.connect("clicked()", dialog.reject)
        if dialog.exec() != qt.QDialog.Accepted:
            return None
        idx = organCombo.currentIndex
        organRow = config[idx]
        nmesh = int(organRow["nmesh"])
        volumes = self._getKermaVolumesForMesh(nmesh)
        if not volumes:
            slicer.util.errorDisplay(
                f"No se encontraron componentes KERMA para mesh {nmesh}."
            )
            return None
        rbeMapped = []
        for node, vol_name in volumes:
            if "Boron" in vol_name or "boro" in vol_name.lower():
                rbeMapped.append(organRow["rbe_boron"])
            elif "Fast" in vol_name:
                rbeMapped.append(organRow["rbe_fast"])
            elif "Thermal" in vol_name or "thn" in vol_name.lower():
                rbeMapped.append(organRow["rbe_thermal"])
            elif "Gamma" in vol_name or "photon" in vol_name.lower():
                rbeMapped.append(organRow["rbe_gamma"])
            else:
                rbeMapped.append(1.0)
        return (organRow, rbeMapped, volumes)

    def onCenterBody(self):
        """Centra la imagen de referencia en el origen RAS.

        Aplica un transform lineal de traslacion a la CT de referencia, a la
        segmentacion, a las mascaras del modulo y al modelo 'Paciente'. El
        offset sale de las dimensiones y el spacing de la imagen de
        referencia (no del bbox del BODY).
        """
        import traceback

        try:
            refNode = self.referenceSelector.currentNode()
            structNode = self.structureSelector.currentNode()
            logic = BNCTLogic()
            result = logic.centerBody(
                referenceNode=refNode, structureNode=structNode
            )
            if result is None:
                slicer.util.errorDisplay(
                    "No se pudo centrar la imagen de referencia.\n\n"
                    "Selecciona la imagen de referencia (pestana Dosimetry).",
                )
                return
            offset = result["offset"]
            dims = result["dims"]
            spacing = result["spacing"]
            origin = result.get("origin")
            img_bounds = result["image_bounds"]
            bounds = result["bounds"]
            if result["body_found"]:
                body_line = (
                    f"BODY '{result['segment_name']}' ({result['source']})"
                )
            else:
                body_line = "BODY: no encontrado en la escena (solo la imagen)."
            if bounds is None:
                body_bounds = "Bounds finales BODY: no disponibles.\n"
            else:
                body_bounds = (
                    f"Bounds finales BODY (mm):\n"
                    f"  X: [{bounds[0]:.3f}, {bounds[1]:.3f}]\n"
                    f"  Y: [{bounds[2]:.3f}, {bounds[3]:.3f}]\n"
                    f"  Z: [{bounds[4]:.3f}, {bounds[5]:.3f}]\n"
                )
            if origin is None:
                origin_line = "Origin de la imagen: no disponible\n"
            else:
                origin_line = (
                    f"Origin de la imagen (ox,oy,oz): {origin[0]:.4f}, "
                    f"{origin[1]:.4f}, {origin[2]:.4f}\n"
                )
            msg = (
                f"{body_line} | Imagen de referencia: "
                f"{result['reference_name']}\n\n"
                f"Dimensiones (Nx,Ny,Nz): {dims[0]}, {dims[1]}, {dims[2]}\n"
                f"Spacing (dx,dy,dz): {spacing[0]:.4f}, {spacing[1]:.4f}, "
                f"{spacing[2]:.4f}\n"
                f"{origin_line}"
                f"Traslacion aplicada (mm): X={offset[0]:.2f}, "
                f"Y={offset[1]:.2f}, Z={offset[2]:.2f}\n"
                f"{body_bounds}"
                f"Bounds finales imagen (mm):\n"
                f"  X: [{img_bounds[0]:.3f}, {img_bounds[1]:.3f}]\n"
                f"  Y: [{img_bounds[2]:.3f}, {img_bounds[3]:.3f}]\n"
                f"  Z: [{img_bounds[4]:.3f}, {img_bounds[5]:.3f}]\n"
                f"Nodos movidos: {result['n_nodes']}\n"
                f"Transform: {result['transform_name']}"
            )
            slicer.util.infoDisplay(msg, windowTitle="Centrar BODY")
        except Exception:
            _log("Error en onCenterBody:\n" + traceback.format_exc())
            slicer.util.errorDisplay("Error centrando el BODY. Ver _bnct_debug.txt.")

    def onUncenterBody(self):
        """Quita el transform de centrado y restaura la posicion original."""
        import traceback

        try:
            released = BNCTLogic().uncenterBody()
            slicer.util.infoDisplay(
                f"Se libero el transform de centrado en {released} nodo(s).",
                windowTitle="Quitar centrado BODY",
            )
        except Exception:
            _log("Error en onUncenterBody:\n" + traceback.format_exc())
            slicer.util.errorDisplay(
                "Error quitando el centrado del BODY. Ver _bnct_debug.txt."
            )

    def onPhysDose(self):
        result = self._buildOrganDialog("Dosis Fisica", ask_rbe=False)
        if result is None:
            return
        organRow, _, volumes = result
        nodes = [node for node, _ in volumes]
        logic = BNCTLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculatePhysicalDose(
                nodes=nodes,
                b10_conc=self.boronConcentration.value,
                ratio_tb=organRow["ratio"],
                tir=self.factorTIR.value,
                scale=self.scaleFactor.value,
                peso1=self.peso1Spin.value,
                peso2=self.peso2Spin.value,
            )
            slicer.util.infoDisplay("Dosis Fisica calculada correctamente.")
        except Exception as e:
            slicer.util.errorDisplay(f"Error: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def onRBEDose(self):
        result = self._buildOrganDialog("Dosis Pesada por RBE", ask_rbe=True)
        if result is None:
            return
        organRow, rbeMapped, volumes = result
        nodes = [node for node, _ in volumes]
        logic = BNCTLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculateRBEDose(
                nodes=nodes,
                rbe_values=rbeMapped,
                b10_conc=self.boronConcentration.value,
                ratio_tb=organRow["ratio"],
                tir=self.factorTIR.value,
                scale=self.scaleFactor.value,
                peso1=self.peso1Spin.value,
                peso2=self.peso2Spin.value,
            )
            slicer.util.infoDisplay("Dosis RBE calculada correctamente.")
        except Exception as e:
            slicer.util.errorDisplay(f"Error: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def _buildIsoDialog(self):
        config = self._organ_config
        if not config:
            slicer.util.errorDisplay("No hay configuracion de organos disponible.")
            return None
        if not self._models_data:
            slicer.util.errorDisplay(
                "No se pudieron cargar los modelos radiobiologicos."
            )
            return None
        models = self._models_data.get("isoeffective", [])
        if not models:
            slicer.util.errorDisplay("No hay modelos isoefectivos en models.jsonc.")
            return None
        dialog = qt.QDialog(self.parent)
        dialog.setWindowTitle("Dosis Isoefectiva")
        dialog.setMinimumWidth(520)
        layout = qt.QVBoxLayout(dialog)
        layout.setSpacing(6)
        organLay = qt.QHBoxLayout()
        organLabel = qt.QLabel("Organo:")
        organCombo = qt.QComboBox()
        organCombo.setMinimumWidth(200)
        for row in config:
            organCombo.addItem(f"{row['name']}  [mesh {int(row['nmesh'])}]")
        organLay.addWidget(organLabel)
        organLay.addWidget(organCombo, 1)
        layout.addLayout(organLay)
        meshInfo = qt.QLabel()
        meshInfo.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(meshInfo)
        modelLay = qt.QHBoxLayout()
        modelLabel = qt.QLabel("Modelo:")
        modelCombo = qt.QComboBox()
        modelCombo.setMinimumWidth(300)
        for m in models:
            modelCombo.addItem(f"{m['id']}. {m['name']}")
        modelLay.addWidget(modelLabel)
        modelLay.addWidget(modelCombo, 1)
        layout.addLayout(modelLay)
        infoLabel = qt.QLabel()
        infoLabel.setWordWrap(True)
        infoLabel.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(infoLabel)
        paramsGroup = qt.QGroupBox("Parametros del modelo")
        paramsLay = qt.QFormLayout(paramsGroup)
        paramDisplays = {}
        paramKeys = [
            ("alpha_r", "α_r"),
            ("beta_r", "β_r"),
            ("alpha_boro", "α_boro"),
            ("beta_boro", "β_boro"),
            ("alpha_thn", "α_thn"),
            ("beta_thn", "β_thn"),
            ("alpha_fast", "α_fast"),
            ("beta_fast", "β_fast"),
            ("tof", "tof"),
            ("tos", "tos"),
            ("pf_g", "pf_g"),
            ("ps_g", "ps_g"),
            ("pf_bnct", "pf_bnct"),
            ("ps_bnct", "ps_bnct"),
        ]
        for key, display in paramKeys:
            lbl = qt.QLabel()
            paramsLay.addRow(f"{display}:", lbl)
            paramDisplays[key] = lbl
        layout.addWidget(paramsGroup)
        bratioLay = qt.QHBoxLayout()
        bratioLabel = qt.QLabel("Tumor/blood B10 ratio:")
        bratioSpin = qt.QDoubleSpinBox()
        bratioSpin.setRange(0, 100)
        bratioSpin.setDecimals(2)
        bratioSpin.value = 3.5
        bratioSpin.setFixedWidth(120)
        bratioLay.addWidget(bratioLabel)
        bratioLay.addWidget(bratioSpin)
        bratioLay.addStretch()
        layout.addLayout(bratioLay)

        def update_model():
            idx = modelCombo.currentIndex
            if idx < 0 or idx >= len(models):
                return
            m = models[idx]
            infoLabel.text = m.get("ref", "")
            for key in paramKeys:
                k = key[0]
                val = m.get(k, 0)
                paramDisplays[k].text = "0" if abs(val) < 1e-10 else f"{val:.6g}"

        def update_organ():
            idx = organCombo.currentIndex
            if idx < 0 or idx >= len(config):
                return
            row = config[idx]
            meshInfo.text = f"Mesh: {int(row['nmesh'])}"

        organCombo.connect("currentIndexChanged(int)", lambda i: update_organ())
        modelCombo.connect("currentIndexChanged(int)", lambda i: update_model())
        update_organ()
        update_model()
        btnLay = qt.QHBoxLayout()
        okBtn = qt.QPushButton("OK")
        cancelBtn = qt.QPushButton("Cancel")
        btnLay.addStretch()
        btnLay.addWidget(okBtn)
        btnLay.addWidget(cancelBtn)
        layout.addLayout(btnLay)
        okBtn.connect("clicked()", dialog.accept)
        cancelBtn.connect("clicked()", dialog.reject)
        if dialog.exec() != qt.QDialog.Accepted:
            return None
        orgIdx = organCombo.currentIndex
        mdlIdx = modelCombo.currentIndex
        organRow = config[orgIdx]
        nmesh = int(organRow["nmesh"])
        modelParams = models[mdlIdx]
        bratio = bratioSpin.value
        volumes = self._getKermaVolumesForMesh(nmesh)
        if not volumes:
            slicer.util.errorDisplay(
                f"No se encontraron componentes KERMA para mesh {nmesh}."
            )
            return None
        paramsArr = [
            modelParams["alpha_r"],
            modelParams["beta_r"],
            modelParams["alpha_boro"],
            modelParams["beta_boro"],
            modelParams["alpha_thn"],
            modelParams["beta_thn"],
            modelParams["alpha_fast"],
            modelParams["beta_fast"],
            modelParams["tof"],
            modelParams["tos"],
            modelParams["pf_g"],
            modelParams["ps_g"],
            modelParams["pf_bnct"],
            modelParams["ps_bnct"],
        ]
        return (paramsArr, bratio, volumes, organRow)

    def onIsoDose(self):
        result = self._buildIsoDialog()
        if result is None:
            return
        paramsArr, bratio, volumes, organRow = result
        nodes = [node for node, _ in volumes]
        compArrays = {"boro": None, "thn": None, "fast": None, "gamma": None}
        for node, vol_name in volumes:
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if "Boron" in vol_name or "boro" in vol_name.lower():
                compArrays["boro"] = (
                    arr if compArrays["boro"] is None else compArrays["boro"] + arr
                )
            elif "Thermal" in vol_name or "thn" in vol_name.lower():
                compArrays["thn"] = (
                    arr if compArrays["thn"] is None else compArrays["thn"] + arr
                )
            elif "Fast" in vol_name:
                compArrays["fast"] = (
                    arr if compArrays["fast"] is None else compArrays["fast"] + arr
                )
            elif "Gamma" in vol_name or "photon" in vol_name.lower():
                compArrays["gamma"] = (
                    arr if compArrays["gamma"] is None else compArrays["gamma"] + arr
                )
        logic = BNCTLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculateIsoEffectveDose(
                comp_arrays=compArrays,
                params=paramsArr,
                bratio=bratio,
                tir=self.factorTIR.value,
                scale=self.scaleFactor.value,
                origin=nodes[0].GetOrigin(),
                spacing=nodes[0].GetSpacing(),
                transformNodeID=nodes[0].GetTransformNodeID(),
            )
            slicer.util.infoDisplay("Dosis Isoefectiva calculada correctamente.")
        except Exception as e:
            slicer.util.errorDisplay(f"Error: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    # ------------------------------------------------------------------
    # Handlers: DVH Analysis
    # ------------------------------------------------------------------
    def _onDvhTodosToggled(self, checked):
        for i in range(self.dvhScrollLayout.count()):
            w = self.dvhScrollLayout.itemAt(i).widget()
            if isinstance(w, qt.QCheckBox):
                w.blockSignals(True)
                w.checked = checked
                w.blockSignals(False)

    def _populateDvhCheckboxes(self, config):
        while self.dvhScrollLayout.count():
            item = self.dvhScrollLayout.takeAt(0)
            w = item.widget()
            if w:
                w.delete()
        for row in config:
            cb = qt.QCheckBox(f"{row['name']}  [mesh {int(row['nmesh'])}]")
            cb.checked = True
            self.dvhScrollLayout.addWidget(cb)
        self.dvhTodosCheck.blockSignals(True)
        self.dvhTodosCheck.checked = True
        self.dvhTodosCheck.blockSignals(False)

    def _onDvhDoseTypeChanged(self, idx):
        self.dvhModeloCombo.visible = idx == 2

    def _loadModels(self):
        if os.path.exists(_MODELS_PATH):
            try:
                return load_jsonc(_MODELS_PATH).get("isoeffective", [])
            except Exception as e:
                logging.warning(f"No se pudo cargar models.jsonc: {e}")
        return []

    def _loadConfigForDvh(self):
        if os.path.exists(_ORGAN_CONFIG_PATH):
            try:
                data = load_jsonc(_ORGAN_CONFIG_PATH)
                config = data.get("rois", [])
                if config:
                    self._populateDvhCheckboxes(config)
                    return config
            except Exception as e:
                logging.warning(f"No se pudo cargar organ_config.jsonc: {e}")
        if os.path.exists(_MAT_PATH):
            try:
                logic = BNCTLogic()
                config = logic.loadConfigFromMat(_MAT_PATH)
                d = os.path.dirname(_ORGAN_CONFIG_PATH)
                if not os.path.exists(d):
                    os.makedirs(d)
                with open(_ORGAN_CONFIG_PATH, "w", encoding="utf-8") as f:
                    json.dump({"rois": config}, f, indent=2, ensure_ascii=False)
                self._populateDvhCheckboxes(config)
                return config
            except Exception as e:
                logging.warning(f"Fallback MAT fallo: {e}")
        return []

    def onDVHButton(self):
        config = self._dvh_config
        selected = []
        for i in range(self.dvhScrollLayout.count()):
            w = self.dvhScrollLayout.itemAt(i).widget()
            if isinstance(w, qt.QCheckBox) and w.checked and i < len(config):
                selected.append(config[i])
        if not selected:
            slicer.util.errorDisplay("Seleccione al menos un organo.")
            return
        dose_type = self.dvhDoseTypeCombo.currentText
        b10 = self.dvhB10Spin.value
        tir = self.dvhTirSpin.value
        model_params = None
        if dose_type == "Dosis Isoefectiva":
            midx = self.dvhModeloCombo.currentIndex
            if midx < 0 or midx >= len(self._dvh_models):
                slicer.util.errorDisplay("Seleccione un modelo.")
                return
            m = self._dvh_models[midx]
            model_params = [
                m["alpha_r"],
                m["beta_r"],
                m["alpha_boro"],
                m["beta_boro"],
                m["alpha_thn"],
                m["beta_thn"],
                m["alpha_fast"],
                m["beta_fast"],
                m["tof"],
                m["tos"],
                m["pf_g"],
                m["ps_g"],
                m["pf_bnct"],
                m["ps_bnct"],
            ]
        mask_node = self.maskVolumeSelector.currentNode()
        if not mask_node:
            slicer.util.errorDisplay("Seleccione un volumen de mascara.")
            return
        logic = BNCTLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculateDVH(
                organs=selected,
                dose_type=dose_type,
                b10_conc=b10,
                tir=tir,
                model_params=model_params,
                mask_volume=mask_node,
            )
        except Exception as e:
            slicer.util.errorDisplay(f"Error en DVH: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
            _log(f"ERROR: {e}\n{traceback.format_exc()}")
        finally:
            qt.QApplication.restoreOverrideCursor()

    # ------------------------------------------------------------------
    # UTCP: TCP / NTCP / UTCP multi-organo (pestaña DVH Analysis)
    # ------------------------------------------------------------------
    def _loadModelsAll(self):
        if os.path.exists(_MODELS_PATH):
            try:
                return load_jsonc(_MODELS_PATH)
            except Exception as e:
                logging.warning(f"No se pudo cargar models.jsonc: {e}")
        return {}

    def _utcpNtcpCatalog(self):
        """Catalogo NTCP: builtin (mucosa/brain/configurable/riback)."""
        try:
            return builtin_ntcp_models()
        except Exception as e:
            logging.warning(f"Catalogo NTCP fallback: {e}")
            return []

    def _setupUTCPUI(self):
        if _UTCP_IMPORT_ERROR:
            lbl = qt.QLabel(f"UTCP no disponible: {_UTCP_IMPORT_ERROR}")
            self.dvhLayout.addWidget(lbl)
            return
        self._utcp_models_all = self._loadModelsAll()
        self.dvhLayout.addWidget(
            qt.QLabel("UTCP — TCP / NTCP / UTCP multi-órgano (Gy IsoE):")
        )
        self.utcpTable = qt.QTableWidget(0, 5)
        self.utcpTable.setHorizontalHeaderLabels(
            ["ROI", "Tipo", "Modelo IsoE", "Modelo NTCP", "Incluir UTCP"]
        )
        try:
            self.utcpTable.horizontalHeader().setSectionResizeMode(
                qt.QHeaderView.Stretch
            )
        except Exception:
            pass
        self.utcpTable.setMinimumHeight(180)
        self.dvhLayout.addWidget(self.utcpTable)
        rowBtns = qt.QHBoxLayout()
        self.utcpCalcBtn = qt.QPushButton("Calcular UTCP")
        self.utcpOptBtn = qt.QPushButton("Optimizar tiempo (barrido)")
        rowBtns.addWidget(self.utcpCalcBtn)
        rowBtns.addWidget(self.utcpOptBtn)
        self.dvhLayout.addLayout(rowBtns)
        self.utcpCalcBtn.connect("clicked(bool)", self.onUTCPCompute)
        self.utcpOptBtn.connect("clicked(bool)", self.onUTCPOptimize)
        self._utcpPopulateTable()

    def _utcpPopulateTable(self):
        config = self._dvh_config or []
        iso_models = self._utcp_models_all.get("isoeffective", [])
        ntcp_models = self._utcpNtcpCatalog()
        self.utcpTable.blockSignals(True)
        self.utcpTable.setRowCount(len(config))
        for i, row in enumerate(config):
            is_tumor = bool(row.get("tumor"))
            name_item = qt.QTableWidgetItem(str(row.get("name", "")))
            name_item.setFlags(name_item.flags() & ~qt.Qt.ItemIsEditable)
            self.utcpTable.setItem(i, 0, name_item)
            tipo = "Tumor (TCP Provenzano)" if is_tumor else "OAR"
            tipo_item = qt.QTableWidgetItem(tipo)
            tipo_item.setFlags(tipo_item.flags() & ~qt.Qt.ItemIsEditable)
            self.utcpTable.setItem(i, 1, tipo_item)
            # Modelo IsoE (por tejido)
            combo_iso = qt.QComboBox()
            for m in iso_models:
                combo_iso.addItem(f"{m['id']}. {m['name']}")
            default_iso = default_isoe_model_for_roi(
                str(row.get("name", "")), is_tumor, [m["name"] for m in iso_models]
            )
            for j, m in enumerate(iso_models):
                if m["name"] == default_iso:
                    combo_iso.setCurrentIndex(j)
                    break
            self.utcpTable.setCellWidget(i, 2, combo_iso)
            # Modelo NTCP (solo OAR)
            if is_tumor:
                ntcp_item = qt.QTableWidgetItem("— (TCP)")
                ntcp_item.setFlags(ntcp_item.flags() & ~qt.Qt.ItemIsEditable)
                self.utcpTable.setItem(i, 3, ntcp_item)
                inc = qt.QTableWidgetItem()
                inc.setFlags(
                    inc.flags() & ~qt.Qt.ItemIsUserCheckable & ~qt.Qt.ItemIsEnabled
                )
                self.utcpTable.setItem(i, 4, inc)
            else:
                combo_ntcp = qt.QComboBox()
                for nm in ntcp_models:
                    if nm.get("kind") == "riback":
                        tag = " (no disponible)"
                    elif nm.get("provisional_para"):
                        tag = " (PROVISIONAL fuera de su endpoint)"
                    else:
                        tag = ""
                    combo_ntcp.addItem(nm.get("name", "LKB") + tag)
                combo_ntcp.insertItem(0, "(excluir de NTCP)")
                saved = row.get("utcp_ntcp_model")
                if saved:
                    for j, nm in enumerate(ntcp_models):
                        if nm.get("name") == saved:
                            combo_ntcp.setCurrentIndex(j + 1)
                            break
                self.utcpTable.setCellWidget(i, 3, combo_ntcp)
                inc = qt.QTableWidgetItem("Incluir")
                inc.setFlags(inc.flags() | qt.Qt.ItemIsUserCheckable)
                inc.setCheckState(
                    qt.Qt.Checked if row.get("utcp_include") else qt.Qt.Unchecked
                )
                self.utcpTable.setItem(i, 4, inc)
        self.utcpTable.blockSignals(False)

    def _utcpResolveTumor(self, config):
        """Resuelve la ROI tumoral por flag (nunca por indice)."""
        tumors = [r for r in config if r.get("tumor")]
        if not tumors:
            slicer.util.errorDisplay(
                "No hay ROI marcada como tumor en Dosimetry (columna Tumor). "
                "Marque exactamente una ROI tumoral para calcular TCP/UTCP."
            )
            return None
        if len(tumors) > 1:
            names = [r.get("name", "?") for r in tumors]
            sel, ok = qt.QInputDialog.getItem(
                slicer.util.mainWindow(),
                "Multiple ROI tumorales",
                "Hay %d ROIs marcadas como tumor (%s). Seleccione la tumoral "
                "para TCP:" % (len(tumors), ", ".join(names)),
                names,
                0,
                False,
            )
            if not ok:
                return None
            return next(r for r in tumors if r.get("name") == sel)
        return tumors[0]

    def _utcpPersistSelection(self):
        """Guarda utcp_isoe_model / utcp_ntcp_model / utcp_include en jsonc."""
        try:
            data = {}
            if os.path.exists(_ORGAN_CONFIG_PATH):
                data = load_jsonc(_ORGAN_CONFIG_PATH)
            rois = data.get("rois", [])
            by_name = {r.get("name"): r for r in rois}
            for i in range(self.utcpTable.rowCount):
                item = self.utcpTable.item(i, 0)
                if item is None:
                    continue
                name = item.text()
                r = by_name.get(name)
                if r is None:
                    continue
                combo_iso = self.utcpTable.cellWidget(i, 2)
                iso_models = self._utcp_models_all.get("isoeffective", [])
                if combo_iso is not None and 0 <= combo_iso.currentIndex < len(
                    iso_models
                ):
                    r["utcp_isoe_model"] = iso_models[combo_iso.currentIndex]["name"]
                combo_ntcp = self.utcpTable.cellWidget(i, 3)
                if combo_ntcp is not None:
                    ntcp_models = self._utcpNtcpCatalog()
                    idx = combo_ntcp.currentIndex - 1
                    r["utcp_ntcp_model"] = (
                        ntcp_models[idx].get("name")
                        if 0 <= idx < len(ntcp_models)
                        else None
                    )
                inc = self.utcpTable.item(i, 4)
                if inc is not None and (inc.flags() & qt.Qt.ItemIsUserCheckable):
                    r["utcp_include"] = inc.checkState() == qt.Qt.Checked
            d = os.path.dirname(_ORGAN_CONFIG_PATH)
            if not os.path.exists(d):
                os.makedirs(d)
            with open(_ORGAN_CONFIG_PATH, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            logging.warning(f"No se pudo persistir seleccion UTCP: {e}")

    def _utcpIsoParams14(self, model_name):
        for m in self._utcp_models_all.get("isoeffective", []):
            if m["name"] == model_name:
                return [
                    m["alpha_r"],
                    m["beta_r"],
                    m["alpha_boro"],
                    m["beta_boro"],
                    m["alpha_thn"],
                    m["beta_thn"],
                    m["alpha_fast"],
                    m["beta_fast"],
                    m["tof"],
                    m["tos"],
                    m["pf_g"],
                    m["ps_g"],
                    m["pf_bnct"],
                    m["ps_bnct"],
                ]
        raise ValueError("Modelo isoefectivo no encontrado: %r" % model_name)

    def _utcpCollectComps(self, meshes, mask_node):
        """Arma {nmesh: comps} con los tally crudos de cada mesh (suma _C2)."""
        ref_shape = np.array(slicer.util.arrayFromVolume(mask_node)).shape
        comps_map = {}
        for nmesh in sorted(meshes):
            volumes = self._getKermaVolumesForMesh(int(nmesh))
            if not volumes:
                raise ValueError(
                    "No hay componentes KERMA para mesh %d. Calcule la "
                    "dosimetria (pestaña Dosimetry) antes de correr UTCP." % nmesh
                )
            comps = {"boro": None, "thn": None, "fast": None, "gamma": None}
            for node, vol_name in volumes:
                arr = slicer.util.arrayFromVolume(node).astype(np.float64)
                if arr.shape != ref_shape:
                    from scipy.ndimage import zoom

                    factors = [r / a for r, a in zip(ref_shape, arr.shape)]
                    arr = zoom(arr, factors, order=1)
                if "Boron" in vol_name or "boro" in vol_name.lower():
                    comps["boro"] = (
                        arr if comps["boro"] is None else comps["boro"] + arr
                    )
                elif "Thermal" in vol_name or "thn" in vol_name.lower():
                    comps["thn"] = arr if comps["thn"] is None else comps["thn"] + arr
                elif "Fast" in vol_name:
                    comps["fast"] = (
                        arr if comps["fast"] is None else comps["fast"] + arr
                    )
                elif "Gamma" in vol_name or "photon" in vol_name.lower():
                    comps["gamma"] = (
                        arr if comps["gamma"] is None else comps["gamma"] + arr
                    )
            missing = [k for k, v in comps.items() if v is None]
            if missing:
                raise ValueError(
                    "Mesh %d sin componentes: %s. No se sustituyen "
                    "por otras distribuciones." % (nmesh, ", ".join(missing))
                )
            comps_map[int(nmesh)] = comps
        return comps_map

    def _utcpBuildInputs(self):
        """Prepara tumor/oars/comps_map/vol por voxel desde la UI."""
        self._utcpPersistSelection()
        config = self._dvh_config or []
        tumor_row = self._utcpResolveTumor(config)
        if tumor_row is None:
            return None
        mask_node = self.maskVolumeSelector.currentNode()
        if not mask_node:
            slicer.util.errorDisplay(
                "Seleccione un volumen de mascara (Mascaras_Combinadas_Gris)."
            )
            return None
        spacing = mask_node.GetSpacing()
        voxel_cc = float(spacing[0]) * float(spacing[1]) * float(spacing[2]) / 1000.0
        mask_arr = np.array(slicer.util.arrayFromVolume(mask_node))
        iso_models = self._utcp_models_all.get("isoeffective", [])
        ntcp_models = self._utcpNtcpCatalog()

        def vol_of(row):
            m = mask_arr == int(row.get("gray_level", 0))
            return float(m.sum()) * voxel_cc, m

        vol_t, mask_t = vol_of(tumor_row)
        if vol_t <= 0:
            slicer.util.errorDisplay(
                "La mascara tumoral '%s' esta vacia (gray_level %s)."
                % (tumor_row.get("name"), tumor_row.get("gray_level"))
            )
            return None
        iso_name = tumor_row.get("utcp_isoe_model") or default_isoe_model_for_roi(
            str(tumor_row.get("name", "")), True, [m["name"] for m in iso_models]
        )
        tumor = {
            "name": tumor_row.get("name", "TUMOR"),
            "mesh": int(tumor_row.get("nmesh", 1) or 1),
            "mask": mask_t,
            "volume_cc": vol_t,
            "bratio": float(tumor_row.get("ratio", 1.0) or 1.0),
            "isoe_params14": self._utcpIsoParams14(iso_name),
            "tcp_params": (
                builtin_tcp_models(self._utcp_models_all.get("tcp", [])) or [None]
            )[0],
        }
        oars = []
        for i in range(self.utcpTable.rowCount):
            item = self.utcpTable.item(i, 0)
            if item is None:
                continue
            row = next((r for r in config if r.get("name") == item.text()), None)
            if row is None or row.get("tumor"):
                continue
            inc = self.utcpTable.item(i, 4)
            include = bool(
                inc
                and (inc.flags() & qt.Qt.ItemIsUserCheckable)
                and inc.checkState() == qt.Qt.Checked
            )
            if not include:
                continue
            combo_ntcp = self.utcpTable.cellWidget(i, 3)
            idx = combo_ntcp.currentIndex - 1 if combo_ntcp is not None else -1
            if idx < 0 or idx >= len(ntcp_models):
                continue
            nm = dict(ntcp_models[idx])
            vol_o, mask_o = vol_of(row)
            iso_name = row.get("utcp_isoe_model") or default_isoe_model_for_roi(
                str(row.get("name", "")), False, [m["name"] for m in iso_models]
            )
            oars.append(
                {
                    "name": row.get("name", "?"),
                    "mesh": int(row.get("nmesh", 1) or 1),
                    "mask": mask_o,
                    "volume_cc": vol_o,
                    "bratio": float(row.get("ratio", 1.0) or 1.0),
                    "isoe_params14": self._utcpIsoParams14(iso_name),
                    "ntcp_params": nm,
                    "include": True,
                    "provisional": bool(
                        nm.get("provisional") or nm.get("provisional_para")
                    ),
                }
            )
        if not oars:
            slicer.util.messageBox(
                "No hay OARs incluidos en el producto UTCP: UTCP = TCP.",
                windowTitle="UTCP",
            )
        b10 = self.dvhB10Spin.value
        tir = self.dvhTirSpin.value
        meshes = {tumor["mesh"]} | {o["mesh"] for o in oars}
        try:
            comps_map = self._utcpCollectComps(meshes, mask_node)
        except ValueError as e:
            slicer.util.errorDisplay(str(e))
            return None
        return {
            "tumor": tumor,
            "oars": oars,
            "comps_map": comps_map,
            "b10": float(b10),
            "tir": float(tir),
        }

    def onUTCPCompute(self):
        if _UTCP_IMPORT_ERROR:
            slicer.util.errorDisplay(
                "Motor UTCP no disponible: %s" % _UTCP_IMPORT_ERROR
            )
            return
        inputs = self._utcpBuildInputs()
        if inputs is None:
            return
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            res = compute_utcp_point(
                inputs["comps_map"],
                inputs["tumor"],
                inputs["oars"],
                inputs["tir"],
                inputs["b10"],
            )
        except Exception as e:
            slicer.util.errorDisplay(f"Error en UTCP: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
            return
        finally:
            qt.QApplication.restoreOverrideCursor()
        lines = [
            f"t = {res['tir']:.2f} min,  B-10 = {inputs['b10']:.1f} ppm",
            f"TCP  = {res['tcp']:.4f}",
            f"UTCP = {res['utcp']:.4f}",
            f"({res['assumption']})",
            "",
        ]
        for name, info in res["ntcp"].items():
            if info.get("ntcp") is None:
                lines.append(
                    f"NTCP {name}: excluido ({info.get('error', 'no incluido')})"
                )
            else:
                lines.append(
                    f"NTCP {name} [{info.get('model', 'LKB')}] = "
                    f"{info['ntcp']:.4f}  (Dref = {info.get('dref_gy', float('nan')):.2f} Gy IsoE)"
                )
        for w in res["warnings"]:
            lines.append(f"[AVISO] {w}")
        slicer.util.messageBox("\n".join(lines), windowTitle="Resultado UTCP")

    def onUTCPOptimize(self):
        if _UTCP_IMPORT_ERROR:
            slicer.util.errorDisplay(
                "Motor UTCP no disponible: %s" % _UTCP_IMPORT_ERROR
            )
            return
        inputs = self._utcpBuildInputs()
        if inputs is None:
            return
        t_max = max(inputs["tir"], 1.0)
        tvals = np.linspace(1.0, t_max, 40)
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

            def _prog(i, total):
                # Mantiene la UI viva y muestra el avance del barrido.
                mainWindow = slicer.util.mainWindow()
                if mainWindow is not None:
                    mainWindow.statusBar().showMessage(
                        "UTCP: barrido temporal %d/%d" % (i + 1, total)
                    )
                slicer.app.processEvents()

            res = sweep_utcp(
                inputs["comps_map"],
                tvals[0],
                inputs["tumor"],
                inputs["oars"],
                tvals,
                inputs["b10"],
                progress=_prog,
            )
        except Exception as e:
            slicer.util.errorDisplay(f"Error en barrido UTCP: {str(e)}")
            import traceback

            logging.error(traceback.format_exc())
            return
        finally:
            qt.QApplication.restoreOverrideCursor()
        # Tabla + grafico
        try:
            self._utcpPlotSweep(res, inputs)
        except Exception as e:
            logging.warning(f"No se pudo graficar barrido UTCP: {e}")
        feasible = (
            f"\nMaximo factible con restricciones: t = {res['t_feasible']:.2f} min, "
            f"UTCP = {res['utcp_feasible_max']:.4f}"
            if res["t_feasible"] is not None
            else "\n(Sin restricciones dosimetricas configuradas)"
        )
        slicer.util.messageBox(
            f"Tiempo optimo (sin restricciones): t = {res['t_opt']:.2f} min\n"
            f"UTCP maximo = {res['utcp_max']:.4f}\n"
            f"TCP en t_opt = {res['tcp'][int(np.nanargmax(res['utcp']))]:.4f}"
            f"{feasible}\n\nCurvas guardadas en tabla/grafico 'UTCP'.",
            windowTitle="Optimizacion de tiempo de irradiacion",
        )

    def _utcpPlotSweep(self, res, inputs):
        table = slicer.mrmlScene.GetFirstNodeByName("Table_UTCP")
        if table is None:
            table = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", "Table_UTCP")
        table.RemoveAllColumns()

        def col(name, values):
            # AddColumn() sin argumentos devuelve vtkStringArray; para numericos
            # hay que pasar un vtkDoubleArray construido desde numpy.
            arr = numpy_to_vtk(np.asarray(values, dtype=np.float64), deep=1)
            arr.SetName(name)
            table.AddColumn(arr)

        col("tir_min", res["tir_values"])
        col("TCP", res["tcp"])
        col("UTCP", res["utcp"])
        for name, curve in res["ntcp"].items():
            col(f"NTCP_{name}", curve)
        plot = slicer.mrmlScene.GetFirstNodeByName("Chart_UTCP")
        if plot is None:
            plot = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLPlotChartNode", "Chart_UTCP"
            )
        # limpiar series previas del chart
        sn = slicer.mrmlScene.GetNodesByClass("vtkMRMLPlotSeriesNode")
        for i in range(sn.GetNumberOfItems()):
            s = sn.GetItemAsIndex(i)
            if s.GetChartNodeID() == plot.GetID():
                slicer.mrmlScene.RemoveNode(s)
        paleta = [
            (0.10, 0.10, 0.90),
            (0.85, 0.45, 0.05),
            (0.00, 0.60, 0.30),
            (0.55, 0.10, 0.75),
            (0.90, 0.10, 0.60),
            (0.00, 0.60, 0.70),
            (0.45, 0.45, 0.45),
            (0.75, 0.75, 0.00),
        ]
        series = [
            ("TCP", "TCP", (0.10, 0.10, 0.90)),
            ("UTCP", "UTCP", (0.90, 0.10, 0.10)),
        ]
        for k, (name, curve) in enumerate(res["ntcp"].items()):
            series.append((f"NTCP_{name}", f"NTCP {name}", paleta[k % len(paleta)]))
        for cname, title, rgb in series:
            s = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLPlotSeriesNode", f"Plot_UTCP_{title}"
            )
            s.SetAndObserveTableNodeID(table.GetID())
            s.SetXColumnName("tir_min")
            s.SetYColumnName(cname)
            s.SetLabel(title)
            s.SetPlotType(s.PlotTypeLine)
            s.SetColor(*rgb)
            s.SetAndObservePlotChartNodeID(plot.GetID())
        mw = slicer.util.mainWindow()
        if mw is not None:
            mw.statusBar().clearMessage()
        slicer.util.selectModule("Plots")

    # ------------------------------------------------------------------
    # Handlers: Mascaras
    # ------------------------------------------------------------------
    def _masksGetConfig(self):
        """Devuelve la config de ROIs: usa la cargada en la pestana Dosimetry o relee el .mat."""
        config = getattr(self, "_organ_config", None)
        if config:
            return config
        return BNCTLogic().loadConfigFromMat(_MAT_PATH)

    def _masksFillTable(self, result):
        """Rellena la tabla de resultados con el dict devuelto por la logica."""
        rows = result.get("rows", [])
        self.masksTable.setRowCount(len(rows))
        for i, row in enumerate(rows):
            self.masksTable.setItem(i, 0, qt.QTableWidgetItem(str(row["name"])))
            self.masksTable.setItem(i, 1, qt.QTableWidgetItem(str(row["gray_level"])))
            self.masksTable.setItem(i, 2, qt.QTableWidgetItem(str(row["nmesh"])))
            self.masksTable.setItem(i, 3, qt.QTableWidgetItem(str(row["nvox"])))
            estado = row["estado"]
            item = qt.QTableWidgetItem(estado)
            if estado != "OK":
                item.setForeground(qt.QBrush(qt.QColor(200, 30, 30)))
            self.masksTable.setItem(i, 4, item)

    def onGenerateMasks(self):
        """Genera los labelmaps por ROI y verifica superposicion."""
        import traceback

        try:
            node = self.structureSelector.currentNode()
            if node is None:
                slicer.util.errorDisplay(
                    "Selecciona un nodo de segmentacion o un volumen de estructura."
                )
                return
            config = self._masksGetConfig()
            logic = BNCTLogic()
            result = logic.generateMasks(
                node,
                config,
                createIndividual=self.chkMasksIndividual.checked,
                createCombined=self.chkMasksCombined.checked,
                skipEmpty=self.chkMasksSkipEmpty.checked,
                referenceNode=self.segReferenceSelector.currentNode(),
            )
            self._masksFillTable(result)
            self.masksReportText.plainText = logic.buildReport(result)
            try:
                if result.get("combined_gray_node") is not None:
                    self.voxelMaskSelector.setCurrentNode(result["combined_gray_node"])
                # Solo refrescar si hay una mascara gris seleccionada: si se genero
                # sin combinado, el selector puede haber quedado vacio y
                # onLoadVoxelLevels mostraria un error innecesario.
                if self.voxelMaskSelector.currentNode() is not None:
                    self.onLoadVoxelLevels()
            except Exception:
                _log(
                    "onGenerateMasks: no se pudo refrescar el tab Voxelizado:\n"
                    + traceback.format_exc()
                )
            slicer.util.infoDisplay(
                f"Mascaras generadas: {len(result.get('individual_nodes', []))} individuales\n"
                f"Segmentos procesados: {result.get('n_segments', 0)} / {len(result.get('rows', []))} filas\n"
                f"Voxeles con superposicion: {result.get('overlap_voxels', 0)}\n"
                f"Sin ROI en la config: {len(result.get('sin_config', []))}\n"
                f"Vacios: {len(result.get('empty_rois', []))}",
                windowTitle="Mascaras BNCT",
            )
        except Exception:
            _log("Error en onGenerateMasks:\n" + traceback.format_exc())
            slicer.util.errorDisplay("Error generando mascaras. Ver _bnct_debug.txt.")

    def onVerifyOverlap(self):
        """Verifica superposicion sin crear volumenes."""
        import traceback

        try:
            node = self.structureSelector.currentNode()
            if node is None:
                slicer.util.errorDisplay(
                    "Selecciona un nodo de segmentacion o un volumen de estructura."
                )
                return
            logic = BNCTLogic()
            result = logic.checkOverlap(
                node,
                self._masksGetConfig(),
                referenceNode=self.segReferenceSelector.currentNode(),
            )
            self._masksFillTable(result)
            report = logic.buildReport(result)
            self.masksReportText.plainText = report
            slicer.util.infoDisplay(report, windowTitle="Verificacion de Superposicion")
        except Exception:
            _log("Error en onVerifyOverlap:\n" + traceback.format_exc())
            slicer.util.errorDisplay(
                "Error verificando superposicion. Ver _bnct_debug.txt."
            )

    def onCleanupMasks(self):
        """Borra los labelmaps generados."""
        import traceback

        try:
            removed = BNCTLogic().cleanupMasks()
            self.masksTable.setRowCount(0)
            self.masksReportText.plainText = (
                f"Se eliminaron {removed} labelmaps de mascaras."
            )
            slicer.util.infoDisplay(
                f"Se eliminaron {removed} labelmaps.", windowTitle="Limpiar Mascaras"
            )
        except Exception:
            _log("Error en onCleanupMasks:\n" + traceback.format_exc())
            slicer.util.errorDisplay("Error limpiando mascaras. Ver _bnct_debug.txt.")


#
# BNCTLogic - Logica compartida para Dosimetry y DVH Analysis
#


class BNCTLogic(ScriptedLoadableModuleLogic):
    @staticmethod
    def load_jsonc(path):
        """Carga un archivo JSONC (JSON con comentarios // y /* */)."""
        import re

        with open(path) as f:
            text = f.read()
        text = re.sub(r"//.*?\n|/\*.*?\*/", "", text, flags=re.S)
        return json.loads(text)

    def loadModels(self):
        """Carga modelos radiobiologicos desde json/models.jsonc."""
        if not os.path.exists(_MODELS_PATH):
            raise FileNotFoundError(f"No se encuentra {_MODELS_PATH}")
        return self.load_jsonc(_MODELS_PATH)

    def loadConfigFromMat(self, path):
        """Carga configuracion KERMA desde archivo MAT (info_rois.mat)."""
        try:
            import h5py
        except ImportError:
            slicer.util.pip_install("h5py")
            import h5py
        config = []

        def read_ref(cell):
            if isinstance(cell, h5py.h5r.Reference):
                obj = f[cell]
            else:
                obj = cell
            data = obj[()]
            if isinstance(data, bytes):
                return data.decode("latin-1")
            if isinstance(data, np.ndarray) and data.dtype.kind == "u":
                return "".join(chr(int(c)) for c in data.ravel())
            if isinstance(data, np.ndarray):
                arr = data.ravel()
                return arr[0] if len(arr) == 1 else arr
            if isinstance(data, h5py.Group):
                return data
            return float(data) if np.isscalar(data) else data

        def read_rbe(rbe_group):
            if isinstance(rbe_group, h5py.Group):
                return {
                    "rbe_boron": float(read_ref(rbe_group["boron"][0, 0])),
                    "rbe_fast": float(read_ref(rbe_group["fast_neutron"][0, 0])),
                    "rbe_thermal": float(read_ref(rbe_group["nitrogen"][0, 0])),
                    "rbe_gamma": float(read_ref(rbe_group["photon"][0, 0])),
                }
            return {
                "rbe_boron": 1.3,
                "rbe_fast": 3.2,
                "rbe_thermal": 3.2,
                "rbe_gamma": 1.0,
            }

        with h5py.File(path, "r") as f:
            info = f["info_dvh_cnea"]
            num_rows = len(info["nroi"])
            for i in range(num_rows):
                rbe_grp = info["RBE"][i, 0]
                rbe_vals = (
                    read_rbe(rbe_grp)
                    if isinstance(rbe_grp, h5py.Group)
                    else read_rbe(rbe_grp)
                )
                row = {
                    "nroi": int(float(read_ref(info["nroi"][i, 0]))),
                    "nmesh": int(float(read_ref(info["nmesh"][i, 0]))),
                    "name": str(read_ref(info["roi_name"][i, 0])),
                    "rbe_boron": float(rbe_vals["rbe_boron"]),
                    "rbe_fast": float(rbe_vals["rbe_fast"]),
                    "rbe_thermal": float(rbe_vals["rbe_thermal"]),
                    "rbe_gamma": float(rbe_vals["rbe_gamma"]),
                    "ratio": float(read_ref(info["ratio"][i, 0])),
                    "tir": float(read_ref(info["tir"][i, 0])),
                    "tumor": bool(int(float(read_ref(info["tumor"][i, 0])))),
                    "gray_level": int(float(read_ref(info["gray_level"][i, 0]))),
                }
                config.append(row)
        return config

    def saveConfigToMat(self, path, config):
        """Guarda configuracion en formato MAT (HDF5)."""
        try:
            import h5py
        except ImportError:
            slicer.util.pip_install("h5py")
            import h5py
        import numpy as np
        from datetime import datetime

        n = len(config)
        with h5py.File(path, "w") as f:
            grp = f.create_group("info_dvh_cnea")
            for key in ["nroi", "nmesh", "gray_level", "tumor"]:
                grp.create_dataset(key, (n, 1), dtype="i4")
            for key in ["ratio", "tir", "Gy2MU", "vol", "Unit"]:
                grp.create_dataset(key, (n, 1), dtype="f8")
            for key in ["roi_name", "fecha", "X", "Y"]:
                grp.create_dataset(key, (n, 1), dtype=h5py.string_dtype())
            for key in ["Max", "Mean", "Min"]:
                grp.create_dataset(key, (n, 1), dtype="f8")
            for i, row in enumerate(config):
                grp["nroi"][i, 0] = row.get("nroi", i + 1)
                grp["nmesh"][i, 0] = row.get("nmesh", 1)
                grp["gray_level"][i, 0] = row.get("gray_level", 0)
                grp["tumor"][i, 0] = 1 if row.get("tumor", False) else 0
                grp["ratio"][i, 0] = row.get("ratio", 1.0)
                grp["tir"][i, 0] = row.get("tir", 20.0)
                grp["Gy2MU"][i, 0] = 0.0
                grp["vol"][i, 0] = 0.0
                grp["Unit"][i, 0] = 1.0
                grp["roi_name"][i, 0] = str(row.get("name", f"ROI_{i + 1}"))
                grp["fecha"][i, 0] = datetime.now().strftime("%d-%b-%Y")
                grp["X"][i, 0] = ""
                grp["Y"][i, 0] = ""
                grp["Max"][i, 0] = 0.0
                grp["Mean"][i, 0] = 0.0
                grp["Min"][i, 0] = 0.0
            rbe_dt = h5py.special_dtype(ref=h5py.Reference)
            rbe_dataset = grp.create_dataset("RBE", (n, 1), dtype=rbe_dt)
            for i in range(n):
                rbe_name = chr(65 + i) if i < 26 else f"rbe_{i}"
                rbe_grp = grp.create_group(f"#{rbe_name}")
                rbe_grp.create_dataset(
                    "boron", data=float(config[i].get("rbe_boron", 1.3))
                )
                rbe_grp.create_dataset(
                    "fast_neutron", data=float(config[i].get("rbe_fast", 3.2))
                )
                rbe_grp.create_dataset(
                    "nitrogen", data=float(config[i].get("rbe_thermal", 3.2))
                )
                rbe_grp.create_dataset(
                    "photon", data=float(config[i].get("rbe_gamma", 1.0))
                )
                rbe_dataset[i, 0] = rbe_grp.ref

    def saveConfigToJson(self, path, config):
        """Guarda configuracion en JSON plano."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(config, f, indent=2, ensure_ascii=False)

    def run(self, neutron_files, gamma_path, refNode, rbe_factors, suffix=""):
        """Procesa archivos PHITS y crea volúmenes KERMA en Slicer."""
        logging.info(f'Iniciando proceso de dosis multicampo (suffix="{suffix}")...')

        if not neutron_files:
            if not suffix:
                slicer.util.errorDisplay(
                    "No se encontraron archivos mesh-N.out (neutrones) en la carpeta."
                )
            return

        print(f"--- Archivos PHITS encontrados (suffix='{suffix}') ---")
        print(f"  KERMA ({len(neutron_files)}):")
        for f in neutron_files:
            print(f"    - {os.path.basename(f)}")
        print(
            f"  Gamma: {os.path.basename(gamma_path) if gamma_path else 'No encontrado'}"
        )

        final_origin = None
        final_spacing = None

        for i, npath in enumerate(neutron_files, start=1):
            print(f"--- Procesando {os.path.basename(npath)} (Field {i}) ---")

            n_data = self.parsePHITS(npath)

            if n_data["particle"] != "neutron":
                logging.warning(
                    f"  Saltando {os.path.basename(npath)}: partícula={n_data['particle']}"
                )
                continue

            if final_origin is None and not refNode:
                final_origin = n_data["origin"]
                final_spacing = n_data["spacing"]

            if refNode:
                n_data["user_offset"] = rbe_factors["offset"]
                for j in range(len(n_data["data"])):
                    n_data["data"][j] = self.interpolateDose(
                        n_data["data"][j], n_data, refNode
                    )
                final_origin = refNode.GetOrigin()
                final_spacing = refNode.GetSpacing()

            print(
                f"Field {i}: final_origin={final_origin}, final_spacing={final_spacing}"
            )

            comp_names = ["Boron", "FastNeutron", "ThermalNeutron"]
            for j, comp in enumerate(comp_names):
                vol_name = f"{comp}_Kerma{i}{suffix}"
                self.createDoseVolume(
                    n_data["data"][j],
                    final_origin,
                    final_spacing,
                    vol_name,
                    transformNodeID=refNode.GetTransformNodeID() if refNode else None,
                )

        has_gamma = False
        gamma_name = "Gamma_Photon"
        if suffix:
            gamma_name = f"Gamma_Photon{suffix}"
        if gamma_path and os.path.exists(gamma_path):
            print(f"--- Procesando gamma: {os.path.basename(gamma_path)} ---")
            g_data = self.parsePHITS(gamma_path)

            if refNode:
                g_data["user_offset"] = rbe_factors["offset"]
                g_data["data"][0] = self.interpolateDose(
                    g_data["data"][0], g_data, refNode
                )
                final_origin = refNode.GetOrigin()
                final_spacing = refNode.GetSpacing()
                print(
                    f"Gamma: final_origin={final_origin}, final_spacing={final_spacing}"
                )

            self.createDoseVolume(
                g_data["data"][0],
                final_origin,
                final_spacing,
                gamma_name,
                transformNodeID=refNode.GetTransformNodeID() if refNode else None,
            )
            has_gamma = True
        else:
            if not suffix:
                logging.warning(
                    "No se encontró archivo gamma. Solo se cargarán campos de neutrones."
                )

        self._organizeInSubjectHierarchy(len(neutron_files), has_gamma)

        print(f"--- Carga completa (suffix='{suffix}') ---")
        print(f"  KERMA components: {len(neutron_files)}")
        print(
            f"  Volúmenes creados: {len(neutron_files) * 3 + (1 if has_gamma else 0)}"
        )
        logging.info("Proceso multicampo finalizado.")

    def _getKermaVolumeNames(self):
        """Encuentra los nombres de volúmenes KERMA en la escena."""
        boro_nodes, fast_nodes, thermal_nodes = [], [], []
        gamma_node = None

        for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
            name = node.GetName()
            if name.startswith("Boron_Kerma"):
                boro_nodes.append(node)
            elif name.startswith("FastNeutron_Kerma"):
                fast_nodes.append(node)
            elif name.startswith("ThermalNeutron_Kerma"):
                thermal_nodes.append(node)
            elif name == "Gamma_Photon":
                gamma_node = node

        return boro_nodes, fast_nodes, thermal_nodes, gamma_node

    def _sumVolumes(self, nodes):
        """Suma los arrays de varios volúmenes (multi-campo)."""
        if not nodes:
            return None
        total = slicer.util.arrayFromVolume(nodes[0]).astype(np.float64)
        for n in nodes[1:]:
            total += slicer.util.arrayFromVolume(n).astype(np.float64)
        return total

    def _isBoronComponent(self, name):
        return "Boron" in name or "boro" in name.lower()

    def calculatePhysicalDose(
        self,
        nodes,
        b10_conc=15.0,
        ratio_tb=3.5,
        tir=1.0,
        scale=1.0,
        peso1=1.0,
        peso2=1.0,
    ):
        """Dosis Física = suma de componentes KERMA seleccionados con B10*ratio."""
        if not nodes:
            raise ValueError("No hay volúmenes KERMA seleccionados.")

        total = np.zeros_like(slicer.util.arrayFromVolume(nodes[0]).astype(np.float64))
        for node in nodes:
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if self._isBoronComponent(node.GetName()):
                arr = arr * b10_conc * ratio_tb
            if "_C2" in node.GetName():
                arr *= peso2
            else:
                arr *= peso1
            total += arr

        phys_dose = total * tir * scale

        origin = nodes[0].GetOrigin()
        spacing = nodes[0].GetSpacing()
        self.createDoseVolume(
            phys_dose,
            origin,
            spacing,
            "Dosis_Fisica_BNCT",
            transformNodeID=nodes[0].GetTransformNodeID(),
        )

        print(f"Dosis Física: max={np.max(phys_dose):.4f}, min={np.min(phys_dose):.4f}")

    def calculateRBEDose(
        self,
        nodes,
        rbe_values,
        b10_conc=15.0,
        ratio_tb=3.5,
        tir=1.0,
        scale=1.0,
        peso1=1.0,
        peso2=1.0,
    ):
        """Dosis RBE = suma de componentes KERMA ponderados por RBE."""
        if not nodes:
            raise ValueError("No hay volúmenes KERMA seleccionados.")
        if len(nodes) != len(rbe_values):
            raise ValueError("Cantidad de nodos y valores RBE no coincide.")

        total = np.zeros_like(slicer.util.arrayFromVolume(nodes[0]).astype(np.float64))
        for node, rbe in zip(nodes, rbe_values):
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if self._isBoronComponent(node.GetName()):
                arr = arr * b10_conc * ratio_tb
            if "_C2" in node.GetName():
                arr *= peso2
            else:
                arr *= peso1
            total += arr * rbe

        rbe_dose = total * tir * scale

        origin = nodes[0].GetOrigin()
        spacing = nodes[0].GetSpacing()
        self.createDoseVolume(
            rbe_dose,
            origin,
            spacing,
            "Dosis_RBE_BNCT",
            transformNodeID=nodes[0].GetTransformNodeID(),
        )

        print(f"Dosis RBE: max={np.max(rbe_dose):.4f}, min={np.min(rbe_dose):.4f}")

    def _leaCatchesideG(self, tof, tos, pf, ps, time_val):
        """Calcula el factor G de Lea-Catcheside para un tiempo dado."""
        xf = tof / time_val
        xs = tos / time_val
        Gf = 2 * xf * (1 - xf * (1 - np.exp(-1 / xf)))
        Gs = 2 * xs * (1 - xs * (1 - np.exp(-1 / xs)))
        return Gs - (pf * Gf + ps * Gs - ps * Gf)

    def calculateIsoEffectveDose(
        self,
        comp_arrays,
        params,
        bratio,
        tir=1.0,
        scale=1.0,
        origin=(0, 0, 0),
        spacing=(1, 1, 1),
        transformNodeID=None,
    ):
        """Dosis Isoefectiva = modelo LQ con Lea-Catcheside."""
        if comp_arrays["boro"] is None:
            raise ValueError("Componente Boron no encontrado.")

        alpha_r, beta_r = params[0], params[1]
        alpha_boro, beta_boro = params[2], params[3]
        alpha_thn, beta_thn = params[4], params[5]
        alpha_fast, beta_fast = params[6], params[7]
        tof, tos = params[8], params[9]
        pf_g, ps_g = params[10], params[11]
        pf_bnct, ps_bnct = params[12], params[13]

        alphabeta_r = alpha_r / beta_r if beta_r != 0 else 1e10

        alpha = np.array([alpha_boro, alpha_thn, alpha_fast, alpha_r, alpha_r])
        beta_mat = np.array([beta_boro, beta_thn, beta_fast, beta_r, beta_r])

        Gr = self._leaCatchesideG(tof, tos, pf_g, 1 - pf_g, 30.0)

        xf = tof / tir
        xs = tos / tir
        Gf = 2 * xf * (1 - xf * (1 - np.exp(-1 / np.maximum(xf, 1e-10))))
        Gs = 2 * xs * (1 - xs * (1 - np.exp(-1 / np.maximum(xs, 1e-10))))

        D_boro = comp_arrays["boro"] * bratio
        D_thn = (
            comp_arrays["thn"]
            if comp_arrays["thn"] is not None
            else np.zeros_like(D_boro)
        )
        D_fast = (
            comp_arrays["fast"]
            if comp_arrays["fast"] is not None
            else np.zeros_like(D_boro)
        )
        D_g = (
            comp_arrays["gamma"]
            if comp_arrays["gamma"] is not None
            else np.zeros_like(D_boro)
        )

        Dt = D_boro + D_thn + D_fast + D_g

        with np.errstate(divide="ignore", invalid="ignore"):
            fb = np.divide(D_boro, Dt, out=np.zeros_like(Dt), where=Dt > 0)
            fn = np.divide(D_thn + D_fast, Dt, out=np.zeros_like(Dt), where=Dt > 0)
            fg = np.divide(D_g, Dt, out=np.zeros_like(Dt), where=Dt > 0)

        denom_12 = fb + fn
        a_12 = np.divide(fb, denom_12, out=np.zeros_like(fb), where=denom_12 > 0)
        a_21 = np.divide(fn, denom_12, out=np.zeros_like(fn), where=denom_12 > 0)

        denom_13 = fb + fg
        a_13 = np.divide(fb, denom_13, out=np.zeros_like(fb), where=denom_13 > 0)
        a_31 = np.divide(fg, denom_13, out=np.zeros_like(fg), where=denom_13 > 0)

        denom_23 = fn + fg
        a_23 = np.divide(fn, denom_23, out=np.zeros_like(fn), where=denom_23 > 0)
        a_32 = np.divide(fg, denom_23, out=np.zeros_like(fg), where=denom_23 > 0)

        G_11 = pf_bnct * Gf + ps_bnct * Gs
        G_22 = G_11
        G_33 = pf_g * Gf + ps_g * Gs

        G_12 = Gs - (a_12 * pf_bnct + a_21 * pf_bnct) * (Gs - Gf)
        G_23 = Gs - (a_23 * pf_bnct + a_32 * pf_g) * (Gs - Gf)
        G_31 = Gs - (a_31 * pf_g + a_13 * pf_bnct) * (Gs - Gf)

        D_lineal = (
            alpha[0] * D_boro + alpha[1] * D_thn + alpha[2] * D_fast + alpha[3] * D_g
        ) * tir

        D_cuad = (
            G_11 * beta_mat[0] * D_boro**2
            + G_22 * beta_mat[1] * D_thn**2
            + G_22 * beta_mat[2] * D_fast**2
            + G_33 * beta_mat[3] * D_g**2
            + 2 * G_12 * np.sqrt(beta_mat[0] * beta_mat[1]) * D_boro * D_thn
            + 2 * G_12 * np.sqrt(beta_mat[0] * beta_mat[2]) * D_boro * D_fast
            + 2 * G_31 * np.sqrt(beta_mat[0] * beta_mat[3]) * D_boro * D_g
            + 2 * G_11 * np.sqrt(beta_mat[1] * beta_mat[2]) * D_thn * D_fast
            + 2 * G_23 * np.sqrt(beta_mat[1] * beta_mat[3]) * D_thn * D_g
            + 2 * G_23 * np.sqrt(beta_mat[2] * beta_mat[3]) * D_fast * D_g
        ) * tir**2

        DisoE = (
            0.5
            * (alphabeta_r / Gr)
            * (np.sqrt(1 + 4 * Gr * beta_r / alpha_r**2 * (D_lineal + D_cuad)) - 1)
            * scale
        )

        DisoE = np.nan_to_num(DisoE)

        self.createDoseVolume(
            DisoE,
            origin,
            spacing,
            "Dosis_IsoE_BNCT",
            transformNodeID=transformNodeID,
        )
        print(f"Dosis Isoefectiva: max={np.max(DisoE):.4f}, min={np.min(DisoE):.4f}")

    def interpolateDose(self, dose_array, source_data, target_node):
        """Remuestrea la dosis usando la matriz de transformación real de Slicer."""
        try:
            from scipy.interpolate import RegularGridInterpolator
        except ImportError:
            slicer.util.pip_install("scipy")
            from scipy.interpolate import RegularGridInterpolator

        nz, ny, nx = dose_array.shape
        ox, oy, oz = source_data["origin"]
        sx, sy, sz = source_data["spacing"]

        x_src = ox + np.arange(nx) * sx
        y_src = oy + np.arange(ny) * sy
        z_src = oz + np.arange(nz) * sz

        t_dims = target_node.GetImageData().GetDimensions()
        t_origin = target_node.GetOrigin()
        t_spacing = target_node.GetSpacing()
        ox_off, oy_off, oz_off = source_data.get("user_offset", (0, 0, 0))

        print(f"--- INTERPOLACIÓN BNCT ---")
        print(f"PHITS Grid Centers Bounds:")
        print(f"  X: [{x_src[0]:.2f}, {x_src[-1]:.2f}] (nx={nx}, sp={sx:.2f})")
        print(f"  Y: [{y_src[0]:.2f}, {y_src[-1]:.2f}] (ny={ny}, sp={sy:.2f})")
        print(f"  Z: [{z_src[0]:.2f}, {z_src[-1]:.2f}] (nz={nz}, sp={sz:.2f})")
        print(f"CT Reference (target):")
        print(f"  Dims: {t_dims}")
        print(f"  Origin: {t_origin}")
        print(f"  Spacing: {t_spacing}")
        print(f"  Offset (user): ({ox_off:.2f}, {oy_off:.2f}, {oz_off:.2f})")

        interp = RegularGridInterpolator(
            (z_src, y_src, x_src),
            dose_array,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

        i_1d = np.arange(t_dims[0])
        j_1d = np.arange(t_dims[1])
        jj, ii = np.meshgrid(j_1d, i_1d, indexing="ij")

        mat = vtk.vtkMatrix4x4()
        target_node.GetIJKToRASMatrix(mat)
        m_np = slicer.util.arrayFromVTKMatrix(mat)

        d_interp = np.zeros((t_dims[2], t_dims[1], t_dims[0]))

        for k_idx in range(t_dims[2]):
            pts_ijk = np.vstack(
                [ii.ravel(), jj.ravel(), np.full(ii.size, k_idx), np.ones(ii.size)]
            )
            pts_ras = m_np @ pts_ijk

            x_pts = pts_ras[0, :] - ox_off
            y_pts = pts_ras[1, :] - oy_off
            z_pts = pts_ras[2, :] - oz_off

            pts_interp = np.vstack([z_pts, y_pts, x_pts]).T
            d_interp[k_idx, :, :] = interp(pts_interp).reshape((t_dims[1], t_dims[0]))

            if (k_idx + 1) % 20 == 0 or k_idx == 0:
                print(f"  Slice {k_idx + 1}/{t_dims[2]} interpolated")

        d_interp = np.nan_to_num(d_interp)

        print(f"Dosis Interp Max: {np.max(d_interp):.6e}, Min: {np.min(d_interp):.6e}")
        return d_interp

    def parsePHITS(self, path):
        """Parser real de PHITS (Traducido de f_lectura_phits_neutron.m)."""
        logging.info(f"Analizando archivo: {path}")

        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        params = {}
        patterns = {
            "xmin": r"xmin\s*=\s*([\d\.\-\+eE]+)",
            "xmax": r"xmax\s*=\s*([\d\.\-\+eE]+)",
            "xdel": r"xdel\s*=\s*([\d\.\-\+eE]+)",
            "ymin": r"ymin\s*=\s*([\d\.\-\+eE]+)",
            "ymax": r"ymax\s*=\s*([\d\.\-\+eE]+)",
            "ydel": r"ydel\s*=\s*([\d\.\-\+eE]+)",
            "zmin": r"zmin\s*=\s*([\d\.\-\+eE]+)",
            "zmax": r"zmax\s*=\s*([\d\.\-\+eE]+)",
            "zdel": r"zdel\s*=\s*([\d\.\-\+eE]+)",
            "nx": r"nx\s*=\s*(\d+)",
            "ny": r"ny\s*=\s*(\d+)",
            "nz": r"nz\s*=\s*(\d+)",
        }

        for line in lines:
            for key, pattern in patterns.items():
                if key not in params:
                    match = re.search(pattern, line)
                    if match:
                        params[key] = float(match.group(1))

        nx, ny, nz = int(params["nx"]), int(params["ny"]), int(params["nz"])
        spacing = (params["xdel"] * 10.0, params["ydel"] * 10.0, params["zdel"] * 10.0)
        origin = (
            (params["xmin"] + params["xdel"] / 2.0) * 10.0,
            (params["ymin"] + params["ydel"] / 2.0) * 10.0,
            (params["zmin"] + params["zdel"] / 2.0) * 10.0,
        )

        particle = "unknown"
        for line in lines:
            m = re.search(r"^\s*part\s*=\s*(\w+)", line)
            if m:
                particle = m.group(1).lower()
                break

        newpage_indices = [
            i for i, line in enumerate(lines) if "newpage:" in line.lower()
        ]
        if not newpage_indices:
            newpage_indices = [
                i for i, line in enumerate(lines) if "hc" in line.lower()
            ]

        if not newpage_indices:
            raise ValueError(
                "No se encontraron bloques de datos (newpage/hc) en el archivo."
            )

        total_pages = len(newpage_indices)
        num_msets = total_pages // nz
        if num_msets > 3:
            num_msets = 3
        elif num_msets == 0 and total_pages > 0:
            num_msets = 1

        import math

        m_matlab = int(math.ceil(nx * ny / 10.0))
        m1 = m_matlab - 1

        msets_data = []

        def robust_float(s):
            try:
                s_clean = s.replace("f", ".").replace("d", "e").replace("D", "e")
                s_clean = re.sub(r"(\d)([\+\-])(\d)", r"\1E\2\3", s_clean)
                return float(s_clean)
            except:
                return 0.0

        for m in range(num_msets):
            mset_array = np.zeros((nz, ny, nx))
            logging.info(f"Cargando mset {m + 1} al estilo MATLAB...")

            for k in range(nz):
                page_idx = m * nz + k
                if page_idx >= len(newpage_indices):
                    break

                idx = newpage_indices[page_idx]

                salto = 1
                if "newpage:" in lines[idx].lower():
                    for offset in range(1, 50):
                        if (idx + offset) < len(lines) and "hc:" in lines[
                            idx + offset
                        ].lower():
                            salto = offset + 1
                            break
                    else:
                        salto = 25

                A_lines = lines[idx + salto : idx + salto + m_matlab]

                all_values = []
                for line in A_lines:
                    vals = [robust_float(v) for v in line.split()]
                    all_values.extend(vals)
                    if len(vals) < 10:
                        all_values.extend([0.0] * (10 - len(vals)))

                B_matlab = np.zeros((nx, ny))
                for j in range(ny):
                    start_ptr = j * nx
                    B_matlab[:, j] = all_values[start_ptr : start_ptr + nx]

                mset_array[k, :, :] = np.flipud(B_matlab.T)

            msets_data.append(mset_array)

        return {
            "data": msets_data,
            "origin": origin,
            "spacing": spacing,
            "shape": (nz, ny, nx),
            "particle": particle,
        }

    def createDoseVolume(self, array, origin, spacing, name, transformNodeID=None):
        """Crea un volumen de dosis con la geometria local dada.

        `origin`/`spacing` son los valores LOCALES del volumen de referencia
        (la interpolacion tambien usa la matriz IJK->RAS local), por eso el
        volumen de dosis debe colgar del mismo transform padre que la CT para
        que no quede desalineado.
        """
        volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", name)

        volumeNode.SetOrigin(origin)
        volumeNode.SetSpacing(spacing)
        slicer.util.updateVolumeFromArray(volumeNode, array)

        volumeNode.SetOrigin(origin)
        volumeNode.SetSpacing(spacing)

        if transformNodeID:
            volumeNode.SetAndObserveTransformNodeID(transformNodeID)

        print(f"--- Volumen Creado: {name} ---")
        print(f"  Array shape: {array.shape}")
        print(f"  Origin: {volumeNode.GetOrigin()}")
        print(f"  Spacing: {volumeNode.GetSpacing()}")

        volumeNode.CreateDefaultDisplayNodes()
        displayNode = volumeNode.GetDisplayNode()

        if displayNode:
            rainbowColorNode = slicer.mrmlScene.GetNodeByID(
                "vtkMRMLColorTableNodeFileRainbow.txt"
            )
            if not rainbowColorNode:
                rainbowColorNode = slicer.util.getNode("Rainbow")

            if rainbowColorNode:
                colorNodeID = rainbowColorNode.GetID()
                if hasattr(displayNode, "SetAndObserveColorNodeID"):
                    displayNode.SetAndObserveColorNodeID(colorNodeID)
                elif hasattr(displayNode, "SetAndUseColorNodeID"):
                    displayNode.SetAndUseColorNodeID(colorNodeID)
                elif hasattr(displayNode, "SetColorNodeID"):
                    displayNode.SetColorNodeID(colorNodeID)
            else:
                logging.warning(
                    f"No se pudo encontrar la escala de colores Rainbow para {name}"
                )

            displayNode.AutoWindowLevelOn()
            displayNode.SetThreshold(
                np.max(array) * 0.01 if np.max(array) > 0 else 0, np.max(array)
            )
            displayNode.SetApplyThreshold(True)

        slicer.util.setSliceViewerLayers(background=volumeNode, fit=True)

        logging.info(f"Volumen '{name}' creado.")
        return volumeNode

    def _cleanupOldVolumes(self):
        prefixes = [
            "Boron_Kerma",
            "FastNeutron_Kerma",
            "ThermalNeutron_Kerma",
            "Gamma_Photon",
            "Raw_Boro_PHITS",
            "Raw_Fast_PHITS",
            "Raw_Thermal_PHITS",
            "Raw_Gamma_PHITS",
            "Dosis_Final_BNCT",
            "Dosis_Fisica_BNCT",
            "Dosis_RBE_BNCT",
            "Dosis_IsoE_BNCT",
        ]
        node_col = slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode")
        for i in range(node_col.GetNumberOfItems()):
            node = node_col.GetItemAsObject(i)
            if node and any(node.GetName().startswith(p) for p in prefixes):
                slicer.mrmlScene.RemoveNode(node)

    def _organizeInSubjectHierarchy(self, field_count, has_gamma):
        shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(
            slicer.mrmlScene
        )
        if not shNode:
            return
        sceneId = shNode.GetSceneItemID()

        def find_or_create_folder(parent_id, folder_name):
            child_ids = vtk.vtkIdList()
            shNode.GetItemChildren(parent_id, child_ids)
            for i in range(child_ids.GetNumberOfIds()):
                child = child_ids.GetId(i)
                if shNode.GetItemName(child) == folder_name:
                    return child
            return shNode.CreateFolderItem(parent_id, folder_name)

        bnctId = find_or_create_folder(sceneId, "BNCT_Dosimetry")
        kermaId = find_or_create_folder(bnctId, "Kerma_Components")

        for i in range(1, field_count + 1):
            kerma_folder_id = find_or_create_folder(kermaId, f"Kerma_{i}")
            for comp in ["Boron", "FastNeutron", "ThermalNeutron"]:
                vol_name = f"{comp}_Kerma{i}"
                vol_node = slicer.mrmlScene.GetFirstNodeByName(vol_name)
                if vol_node:
                    item_id = shNode.GetItemByDataNode(vol_node)
                    if item_id:
                        shNode.SetItemParent(item_id, kerma_folder_id)

        if has_gamma:
            gammaId = find_or_create_folder(bnctId, "Gamma")
            vol_node = slicer.mrmlScene.GetFirstNodeByName("Gamma_Photon")
            if vol_node:
                item_id = shNode.GetItemByDataNode(vol_node)
                if item_id:
                    shNode.SetItemParent(item_id, gammaId)

    def _calcDoseForMesh(
        self,
        nmesh,
        dose_type,
        b10,
        tir,
        ratio,
        rbe_b,
        rbe_f,
        rbe_t,
        rbe_g,
        model_params,
    ):
        volumes = self._getKermaVolumesForMesh(nmesh)
        if not volumes:
            return None
        comp = {"boro": None, "fast": None, "thermal": None, "gamma": None}
        for node, vol_name in volumes:
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if "Boron" in vol_name or "boro" in vol_name.lower():
                comp["boro"] = arr if comp["boro"] is None else comp["boro"] + arr
            elif "Fast" in vol_name:
                comp["fast"] = arr if comp["fast"] is None else comp["fast"] + arr
            elif "Thermal" in vol_name or "thn" in vol_name.lower():
                comp["thermal"] = (
                    arr if comp["thermal"] is None else comp["thermal"] + arr
                )
            elif "Gamma" in vol_name or "photon" in vol_name.lower():
                comp["gamma"] = arr if comp["gamma"] is None else comp["gamma"] + arr
        if comp["boro"] is None:
            return None
        if dose_type == "Dosis Fisica":
            return (
                comp["boro"] * b10 * ratio
                + (comp["fast"] if comp["fast"] is not None else 0)
                + (comp["thermal"] if comp["thermal"] is not None else 0)
                + (comp["gamma"] if comp["gamma"] is not None else 0)
            ) * tir
        elif dose_type == "Dosis RBE":
            return (
                comp["boro"] * b10 * ratio * rbe_b
                + (comp["fast"] if comp["fast"] is not None else 0) * rbe_f
                + (comp["thermal"] if comp["thermal"] is not None else 0) * rbe_t
                + (comp["gamma"] if comp["gamma"] is not None else 0) * rbe_g
            ) * tir
        elif dose_type == "Dosis Isoefectiva" and model_params is not None:
            return self._calcIsoDose(comp, model_params, ratio, tir)
        return None

    def _calcIsoDose(self, comp, model_params, ratio, tir):
        alpha_r, beta_r = model_params[0], model_params[1]
        alpha_boro, beta_boro = model_params[2], model_params[3]
        alpha_thn, beta_thn = model_params[4], model_params[5]
        alpha_fast, beta_fast = model_params[6], model_params[7]
        tof, tos = model_params[8], model_params[9]
        pf_g, ps_g = model_params[10], model_params[11]
        pf_bnct, ps_bnct = model_params[12], model_params[13]
        alphabeta_r = alpha_r / beta_r if beta_r != 0 else 1e10
        D_boro = comp["boro"] * ratio
        D_thn = (
            comp["thermal"] if comp["thermal"] is not None else np.zeros_like(D_boro)
        )
        D_fast = comp["fast"] if comp["fast"] is not None else np.zeros_like(D_boro)
        D_g = comp["gamma"] if comp["gamma"] is not None else np.zeros_like(D_boro)
        Dt = D_boro + D_thn + D_fast + D_g
        with np.errstate(divide="ignore", invalid="ignore"):
            fb = np.divide(D_boro, Dt, out=np.zeros_like(Dt), where=Dt > 0)
            fn = np.divide(D_thn + D_fast, Dt, out=np.zeros_like(Dt), where=Dt > 0)
            fg = np.divide(D_g, Dt, out=np.zeros_like(Dt), where=Dt > 0)

        def lc_g(tof, tos, pf, ps, t):
            xf, xs = tof / t, tos / t
            Gf = 2 * xf * (1 - xf * (1 - np.exp(-1 / max(xf, 1e-10))))
            Gs = 2 * xs * (1 - xs * (1 - np.exp(-1 / max(xs, 1e-10))))
            return Gs - (pf * Gf + ps * Gs - ps * Gf)

        Gr = lc_g(tof, tos, pf_g, 1 - pf_g, 30.0)
        G_11 = lc_g(tof, tos, pf_bnct, ps_bnct, tir)
        G_22 = G_11
        G_33 = lc_g(tof, tos, pf_g, ps_g, tir)
        with np.errstate(divide="ignore", invalid="ignore"):
            a_12 = np.divide(fb, fb + fn, out=np.zeros_like(fb), where=(fb + fn) > 0)
            a_21 = np.divide(fn, fb + fn, out=np.zeros_like(fn), where=(fb + fn) > 0)
            a_13 = np.divide(fb, fb + fg, out=np.zeros_like(fb), where=(fb + fg) > 0)
            a_31 = np.divide(fg, fb + fg, out=np.zeros_like(fg), where=(fb + fg) > 0)
            a_23 = np.divide(fn, fn + fg, out=np.zeros_like(fn), where=(fn + fg) > 0)
            a_32 = np.divide(fg, fn + fg, out=np.zeros_like(fg), where=(fn + fg) > 0)
        G_12 = lc_g(tof, tos, a_12 * pf_bnct + a_21 * pf_bnct, 0, tir)
        G_23 = lc_g(tof, tos, a_23 * pf_bnct + a_32 * pf_g, 0, tir)
        G_31 = lc_g(tof, tos, a_31 * pf_g + a_13 * pf_bnct, 0, tir)
        D_lineal = (
            alpha_boro * D_boro
            + alpha_thn * D_thn
            + alpha_fast * D_fast
            + alpha_r * D_g
        ) * tir
        D_cuad = (
            G_11 * beta_boro * D_boro**2
            + G_22 * beta_thn * D_thn**2
            + G_22 * beta_fast * D_fast**2
            + G_33 * beta_r * D_g**2
            + 2 * G_12 * np.sqrt(beta_boro * beta_thn) * D_boro * D_thn
            + 2 * G_12 * np.sqrt(beta_boro * beta_fast) * D_boro * D_fast
            + 2 * G_31 * np.sqrt(beta_boro * beta_r) * D_boro * D_g
            + 2 * G_11 * np.sqrt(beta_thn * beta_fast) * D_thn * D_fast
            + 2 * G_23 * np.sqrt(beta_thn * beta_r) * D_thn * D_g
            + 2 * G_23 * np.sqrt(beta_fast * beta_r) * D_fast * D_g
        ) * tir**2
        dose = (
            0.5
            * (alphabeta_r / Gr)
            * (np.sqrt(1 + 4 * Gr * beta_r / alpha_r**2 * (D_lineal + D_cuad)) - 1)
        )
        return np.nan_to_num(dose)

    def _cleanupOldDVH(self):
        for cls, prefix in [
            ("vtkMRMLPlotSeriesNode", "DVH_"),
            ("vtkMRMLTableNode", "Table_DVH_"),
        ]:
            coll = slicer.mrmlScene.GetNodesByClass(cls)
            for i in range(coll.GetNumberOfItems()):
                node = coll.GetItemAsObject(i)
                if node and node.GetName().startswith(prefix):
                    slicer.mrmlScene.RemoveNode(node)
        old = slicer.mrmlScene.GetFirstNodeByName("Chart_DVH")
        if old:
            slicer.mrmlScene.RemoveNode(old)

    def calculateDVH(
        self, organs, dose_type, b10_conc, tir, model_params=None, mask_volume=None
    ):
        import re

        _log(f"\n=== DVH START: {len(organs)} organs, {dose_type} ===")

        if mask_volume is None:
            raise ValueError("Debe seleccionar un volumen de mascara.")

        ref_node = None
        for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
            name = node.GetName()
            if any(x in name for x in ["Kerma", "Gamma", "Dosis", "DVH", "Temp"]):
                continue
            try:
                arr = slicer.util.arrayFromVolume(node)
                if arr.ndim == 3:
                    ref_node = node
                    break
            except:
                continue
        if ref_node is None:
            raise ValueError("No se encontro imagen de referencia.")
        ct_shape = slicer.util.arrayFromVolume(ref_node).shape

        mask_array = slicer.util.arrayFromVolume(mask_volume).astype(np.float64)

        self._cleanupOldDVH()
        plotChartNode = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLPlotChartNode", "Chart_DVH"
        )

        all_names = []
        any_dose_computed = False
        global_max_dose = 0.0
        spacing = ref_node.GetSpacing()
        voxel_vol_cc = spacing[0] * spacing[1] * spacing[2] / 1000.0

        for org in organs:
            nmesh = int(org["nmesh"])
            ratio = org["ratio"]
            rbe_b = org["rbe_boron"]
            rbe_f = org["rbe_fast"]
            rbe_t = org["rbe_thermal"]
            rbe_g = org["rbe_gamma"]

            print(
                f"[DVH] Procesando {org['name']} (mesh {nmesh}, gray_level {org.get('gray_level', 0)})..."
            )
            dose = self._calcDoseForMesh(
                nmesh,
                dose_type,
                b10_conc,
                tir,
                ratio,
                rbe_b,
                rbe_f,
                rbe_t,
                rbe_g,
                model_params,
            )
            if dose is None:
                print(f"[DVH]  -> SKIP: no hay kerma volumes para mesh {nmesh}")
                continue
            if dose.shape != ct_shape:
                try:
                    from scipy.ndimage import zoom

                    factors = [
                        max(1, ct_shape[0]) / max(1, dose.shape[0]),
                        max(1, ct_shape[1]) / max(1, dose.shape[1]),
                        max(1, ct_shape[2]) / max(1, dose.shape[2]),
                    ]
                    dose = zoom(dose, factors, order=1)
                except Exception as e:
                    print(f"[DVH]  -> SKIP: zoom fallo - {e}")
                    continue
            any_dose_computed = True

            gl = org.get("gray_level", 0)
            mask = mask_array == gl
            roi_doses = dose[mask]
            if len(roi_doses) == 0:
                print(
                    f"[DVH]  -> SKIP: gray_level={gl} sin voxeles en mascara para '{org['name']}'"
                )
                continue

            roi_doses = roi_doses[~np.isnan(roi_doses)]
            if len(roi_doses) == 0:
                print(f"[DVH]  -> SKIP: todos NaN para '{org['name']}'")
                continue

            sorted_doses = np.sort(roi_doses, axis=None)
            n = len(sorted_doses)
            dvh_vol = np.linspace(100.0, 0.0, n)

            sorted_doses = np.insert(sorted_doses, 0, 0.0)
            dvh_vol = np.insert(dvh_vol, 0, 100.0)
            n = len(sorted_doses)

            vol_cc = voxel_vol_cc * (n - 1)
            Dmax = float(sorted_doses[-1])
            Dmin = float(sorted_doses[0])
            Dmean = float(np.mean(roi_doses))
            if Dmax > global_max_dose:
                global_max_dose = Dmax

            _log(
                f"  OK {org['name']}: {n} vox, vol={vol_cc:.4f}cc, "
                f"Dmin={Dmin:.6e} Dmax={Dmax:.6e} Dmedia={Dmean:.6e}"
            )
            _log(f"    primeras 5 dosis: {sorted_doses[:5].tolist()}")
            _log(f"    ultimas 5 dosis: {sorted_doses[-5:].tolist()}")

            safe_name = re.sub(r"[^a-zA-Z0-9]", "_", org["name"])
            all_names.append(org["name"])

            tn = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLTableNode", f"Table_DVH_{safe_name}"
            )
            table_array = np.column_stack((sorted_doses, dvh_vol))
            slicer.util.updateTableFromArray(tn, table_array, ["Dose", "VolumePct"])

            tn.GetTable().GetColumn(0).SetName("Dose")
            tn.GetTable().GetColumn(1).SetName("VolumePct")
            tn.Modified()

            print(
                f"[DVH] Tabla '{tn.GetName()}': {tn.GetTable().GetNumberOfRows()} filas"
            )
            print(
                f"[DVH] Rango dosis: {sorted_doses[0]:.6e} - {sorted_doses[-1]:.6e} Gy"
            )

            ps = slicer.mrmlScene.AddNewNodeByClass(
                "vtkMRMLPlotSeriesNode", f"DVH_{safe_name}"
            )
            ps.SetAndObserveTableNodeID(tn.GetID())
            ps.SetXColumnName("Dose")
            ps.SetYColumnName("VolumePct")
            ps.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
            ps.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleSolid)
            ps.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleNone)
            ps.SetUniqueColor()
            plotChartNode.AddAndObservePlotSeriesNodeID(ps.GetID())

        if not any_dose_computed:
            raise ValueError("No se pudo calcular dosis para ningun organo.")
        if not all_names:
            raise ValueError("No se encontraron voxeles para ningun organo.")

        _log(f"  global_max_dose = {global_max_dose:.6e}")
        _log(f"  Nombres: {all_names}")

        dose_tag = dose_type.replace("Dosis ", "")
        plotChartNode.SetTitle(
            f"DVH - {', '.join(all_names[:3])}"
            + ("..." if len(all_names) > 3 else f" ({dose_tag})")
        )
        plotChartNode.SetXAxisTitle("Dosis [Gy]")
        plotChartNode.SetYAxisTitle("Volumen [%]")
        plotChartNode.SetLegendVisibility(True)

        xmin = 0.0
        xmax = global_max_dose * 1.05 if global_max_dose > 0 else 1.0
        plotChartNode.SetXAxisRange(xmin, xmax)
        plotChartNode.SetYAxisRange(0, 105)
        _log(f"  Axis ranges: X=[{xmin}, {xmax:.6e}] Y=[0, 105]")

        layoutManager = slicer.app.layoutManager()
        layoutWithPlot = slicer.modules.plots.logic().GetLayoutWithPlot(
            layoutManager.layout
        )
        layoutManager.setLayout(layoutWithPlot)
        plotWidget = layoutManager.plotWidget(0)
        if plotWidget:
            pvn = plotWidget.mrmlPlotViewNode()
            if pvn:
                pvn.SetPlotChartNodeID(plotChartNode.GetID())

        plotChartNode.Modified()
        if plotWidget:
            plotWidget.update()

        _log(f"=== DVH END: {len(all_names)} organs ===")

        slicer.util.infoDisplay(
            f"DVH calculado: {len(all_names)} organos\n"
            f"Dosis max: {global_max_dose:.4e} Gy\n"
            f"Ver tabla _dvh_debug.txt para datos crudos",
            windowTitle="DVH Completo",
        )

    # ------------------------------------------------------------------
    # Mascaras: generacion de labelmaps por ROI
    # ------------------------------------------------------------------
    def generateMasks(
        self,
        structureNode,
        config,
        createIndividual=True,
        createCombined=True,
        skipEmpty=True,
        referenceNode=None,
    ):
        """Genera un labelmap binario por ROI y un labelmap combinado.

        `structureNode` puede ser un vtkMRMLSegmentationNode (ruta principal:
        cada segmento es una ROI, emparejada por NOMBRE) o un
        vtkMRMLScalarVolumeNode (fallback: cada ROI son los voxeles cuyo valor
        coincide con su gray_level).
        `referenceNode` define la geometria (origin/spacing) de los labelmaps
        de salida; si es None se busca un volumen escalar 3D en la escena.
        Retorna un dict con el resultado por ROI y la verificacion de superposicion.
        """
        import traceback

        if structureNode is None:
            raise ValueError(
                "No hay volumen ni segmentacion de estructura seleccionada."
            )

        # Detectar el tipo de nodo: los segmentos del Segment Editor son la via principal
        try:
            is_seg = bool(structureNode.IsA("vtkMRMLSegmentationNode"))
        except Exception:
            is_seg = False

        if is_seg:
            return self._generateMasksFromSegmentation(
                structureNode,
                config,
                createIndividual=createIndividual,
                createCombined=createCombined,
                skipEmpty=skipEmpty,
                referenceNode=referenceNode,
            )

        # Limpiar labelmaps de una corrida anterior (ruta volumen/labelmap).
        self._removeOldMaskNodes(keepNodes=(structureNode, referenceNode))

        arr = slicer.util.arrayFromVolume(structureNode)
        if arr is None or np.asarray(arr).ndim != 3:
            raise ValueError("El volumen de estructura debe ser 3D.")
        arr = np.asarray(arr)

        # Comparacion robusta: el volumen puede ser float, se redondea una sola vez
        arr_f = np.asarray(arr, dtype=np.float64)
        arr_int = np.rint(arr_f).astype(np.int64)

        # Diagnostico: avisar si el volumen no parece un labelmap de gray_level.
        # Si el volumen es una imagen de intensidades (no enteros de etiqueta),
        # arr_int == gray_level solo coincidira en puntos aislados y las mascaras
        # resultantes seran basura.
        uniq_all = np.unique(arr_int)
        nonzero_uniq = uniq_all[uniq_all != 0]
        is_integer_valued = bool(np.all(arr_f == arr_int))
        looks_like_labelmap = is_integer_valued and nonzero_uniq.size <= 64
        if not looks_like_labelmap:
            _log(
                "generateMasks: ADVERTENCIA: el volumen de estructura NO parece un labelmap de gray_level "
                f"(enteros={is_integer_valued}, valores distintos={uniq_all.size}, "
                f"rango=[{int(arr_int.min())},{int(arr_int.max())}]). "
                "arr_int == gray_level solo coincidira en puntos aislados. "
                "Usa un nodo de Segmentacion o un labelmap con los gray_level exactos. "
                f"np.unique(arr_int)[:20]={uniq_all[:20].tolist()}"
            )

        gray_levels = [int(r.get("gray_level", 0)) for r in config]
        max_gl = max(gray_levels) if gray_levels else 0
        combined_dtype = np.uint8 if max_gl <= 255 else np.uint16
        if combined_dtype is not np.uint8:
            _log(
                f"generateMasks: gray_level maximo {max_gl} > 255, combinado en uint16"
            )

        combined = np.zeros(arr.shape, dtype=combined_dtype)
        counts = np.zeros(arr.shape, dtype=np.uint16)

        rows = []
        individual_nodes = []
        skipped_gl0 = []
        # Recolectar entradas para escribir el combinado ordenado: primero los grandes (mayor nvox),
        # despues los pequenos (menor nvox) para que sobrelayen y se vean (p.ej. pulmon dentro de BODY).
        combined_entries = []

        for roi in config:
            gl = int(roi.get("gray_level", 0))
            name = str(roi.get("name", f"ROI_{gl}"))
            nmesh = int(roi.get("nmesh", 0))

            # gray_level 0 es aire: no es una ROI
            if gl == 0:
                skipped_gl0.append(name)
                _log(f"generateMasks: ROI '{name}' con gray_level 0 (aire), se omite")
                continue

            mask = arr_int == gl
            nvox = int(np.count_nonzero(mask))
            estado = "OK" if nvox > 0 else "VACIA"

            counts += mask.astype(np.uint16)
            combined_entries.append((mask, gl, name, nvox))

            if createIndividual and (nvox > 0 or not skipEmpty):
                safe_name = re.sub(r"\s+", "_", name.strip())
                safe_name = re.sub(r"[^0-9A-Za-z_.-]", "_", safe_name) or f"ROI_{gl}"
                color = self._get_distinct_color(gl % 100, 100)
                # Usar el transform de structureNode si existe
                transformID = (
                    structureNode.GetTransformNodeID() if structureNode else None
                )
                node = self._createLabelmapVolume(
                    mask.astype(np.uint8),
                    structureNode,
                    f"{_PREFIX_INDIVIDUAL}{safe_name}",
                    pixel_colors={1: color},
                    transformNodeID=transformID,
                )
                individual_nodes.append(node)

            rows.append(
                {
                    "name": name,
                    "gray_level": gl,
                    "nmesh": nmesh,
                    "nvox": nvox,
                    "estado": estado,
                }
            )

        combined_node = None
        combined_gray_node = None
        if createCombined:
            # Escribir el combinado ordenado: primero los grandes (mayor nvox),
            # despues los pequenos (menor nvox) para que sobrelayen y se vean.
            combined_entries.sort(key=lambda e: e[3], reverse=True)
            for mask, gl, name, nvox in combined_entries:
                if nvox > 0:
                    combined[mask] = gl
            combined_colors = {}
            for row in rows:
                if row["nvox"] > 0:
                    combined_colors[row["gray_level"]] = self._get_distinct_color(
                        row["gray_level"] % 100, 100
                    )
            transformID = structureNode.GetTransformNodeID() if structureNode else None
            combined_node = self._createLabelmapVolume(
                combined,
                structureNode,
                _NAME_COMBINED,
                pixel_colors=combined_colors,
                transformNodeID=transformID,
            )
            # Crear nodo gris combinado para voxelization
            combined_gray_node = self._createLabelmapVolume(
                combined,
                structureNode,
                _NAME_COMBINED_GRAY,
                grayscale=True,
                transformNodeID=transformID,
            )

        # --- Verificacion de superposicion ---
        overlap_voxels = int(np.count_nonzero(counts > 1))
        overlap_rois = []
        if overlap_voxels > 0:
            overlap_rois = self._roisInOverlap(arr_int, rows)

        # --- gray_level duplicado en la configuracion ---
        by_gl = {}
        for row in rows:
            by_gl.setdefault(row["gray_level"], []).append(row["name"])
        duplicates = [
            (gl, names) for gl, names in sorted(by_gl.items()) if len(names) > 1
        ]
        for gl, names in duplicates:
            _log(f"generateMasks: CONFLICTO gray_level {gl} compartido por {names}")

        # --- Valores del volumen que no corresponden a ninguna ROI ---
        unmapped_values = []
        try:
            uniq = np.unique(arr_int)
            gl_set = set(by_gl.keys())
            unmapped_values = [
                int(v) for v in uniq if int(v) != 0 and int(v) not in gl_set
            ]
        except Exception:
            _log(
                "generateMasks: fallo el calculo de valores no mapeados:\n"
                + traceback.format_exc()
            )

        empty_rois = [row["name"] for row in rows if row["nvox"] == 0]

        result = {
            "rows": rows,
            "overlap_voxels": overlap_voxels,
            "overlap_rois": overlap_rois,
            "duplicates": duplicates,
            "unmapped_values": unmapped_values,
            "empty_rois": empty_rois,
            "combined_node": combined_node,
            "combined_gray_node": combined_gray_node,
            "individual_nodes": individual_nodes,
            "total_voxels": int(arr.size),
            "volume_name": structureNode.GetName(),
            "dims": f"{arr.shape[2]} x {arr.shape[1]} x {arr.shape[0]}",
            "n_rois_config": len(config) + len(skipped_gl0),
            "skipped_gl0": skipped_gl0,
            "mode": "volume",
            "sin_config": [],
            "n_segments": 0,
        }

        _log(
            f"generateMasks: total={result['total_voxels']} voxeles, "
            f"superposicion={overlap_voxels}, duplicados={len(duplicates)}, "
            f"no_mapeados={len(unmapped_values)}, vacias={len(empty_rois)}"
        )
        if individual_nodes:
            self._groupMaskNodesInFolder(individual_nodes)
        return result

    def _generateMasksFromSegmentation(
        self,
        structureNode,
        config,
        createIndividual=True,
        createCombined=True,
        skipEmpty=True,
        referenceNode=None,
    ):
        """Genera labelmaps a partir de los segmentos de un vtkMRMLSegmentationNode.

        El emparejamiento segmento -> ROI de la config es POR NOMBRE del segmento,
        comparado con `_normName` (mayusculas, sin acentos, separadores a espacio).

        Los gray_level provienen EXCLUSIVAMENTE de `json/organ_config.jsonc`. Los
        segmentos que no tienen ROI en la config se exportan igualmente como
        labelmaps individuales marcados 'SIN CONFIG', se excluyen del labelmap
        combinado y NO reciben gray_level (quedan con gray_level 0); el modulo no
        inventa valores para ellos.
        """
        import traceback

        # Indice de la config por nombre normalizado, y por nombre "compacto"
        # (sin espacios) para tolerar variantes como 'RT OPTIC NERVE' vs 'RTOPTIC NERVE'.
        def _compact(s):
            return _normName(s).replace(" ", "")

        cfg_by_name = {}
        for row in config:
            cfg_by_name[_normName(row.get("name", ""))] = row
            cfg_by_name.setdefault(_compact(row.get("name", "")), row)

        # El aire vive en la clave 'aire' del jsonc (no en 'rois'): fila pseudo
        # con gray_level 0 para que el aire quede como fondo (universo U=1).
        _AIR_ROW = {"name": "AIRE", "gray_level": 0, "universo": 1, "nmesh": -1}
        for _alias in ("AIRE", "AIR", "EXTERIOR", "AIREEXTERIOR"):
            cfg_by_name.setdefault(_alias, _AIR_ROW)

        def _matchConfigRow(seg_name):
            """Devuelve la fila de config para un segmento, o None."""
            key = _normName(seg_name)
            if key in cfg_by_name:
                return cfg_by_name[key]
            return cfg_by_name.get(_compact(seg_name))

        # --- Volumen de referencia: define la geometria de salida ---
        ref_node = referenceNode
        if ref_node is None:
            for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
                name = node.GetName() or ""
                if any(x in name for x in ["Kerma", "Gamma", "Dosis", "DVH", "Temp"]):
                    continue
                # Un vtkMRMLLabelMapVolumeNode ES un vtkMRMLScalarVolumeNode: si no se
                # descarta, ref_node puede terminar siendo un labelmap del propio
                # modulo (p.ej. Mascaras_Combinadas_Gris de una corrida anterior),
                # que luego queda protegido por _removeOldMaskNodes y hace que el
                # combinado nuevo se cree como _1, _2... acumulando duplicados.
                if node.IsA("vtkMRMLLabelMapVolumeNode"):
                    continue
                if (
                    name.startswith(_PREFIX_INDIVIDUAL)
                    or name == _NAME_COMBINED
                    or name.startswith(_NAME_COMBINED + "_")
                    or name == _NAME_COMBINED_GRAY
                    or name.startswith(_NAME_COMBINED_GRAY + "_")
                ):
                    continue
                try:
                    arr_ref = slicer.util.arrayFromVolume(node)
                    if arr_ref is not None and np.asarray(arr_ref).ndim == 3:
                        ref_node = node
                        break
                except Exception:
                    continue
        if ref_node is None:
            raise ValueError(
                "Selecciona un volumen de referencia para definir la geometria de salida."
            )

        # Limpiar labelmaps de una corrida anterior para no dejar duplicados
        # (_1, _2...) que confundian al tab Voxelizado.
        self._removeOldMaskNodes(keepNodes=(structureNode, ref_node))

        # --- Representacion binaria: asegurarla antes de extraer arrays ---
        if hasattr(structureNode, "CreateBinaryLabelmapRepresentation"):
            try:
                structureNode.CreateBinaryLabelmapRepresentation()
            except Exception:
                _log(
                    "generateMasks: CreateBinaryLabelmapRepresentation fallo:\n"
                    + traceback.format_exc()
                )

        segmentation = structureNode.GetSegmentation()
        n_seg = segmentation.GetNumberOfSegments()

        # Lista de (seg_id, seg_name) sin tocar aun los arrays
        seg_items = []
        for i in range(n_seg):
            seg_id = segmentation.GetNthSegmentID(i)
            seg = segmentation.GetSegment(seg_id)
            seg_name = seg.GetName() if seg else seg_id
            seg_items.append((seg_id, seg_name))

        # Los gray_level salen UNICAMENTE de la config (organ_config.jsonc).
        # No se inventan valores para segmentos sin ROI en la config.
        gray_levels = [int(r.get("gray_level", 0)) for r in config]
        max_gl = max(gray_levels) if gray_levels else 0
        combined_dtype = np.uint16 if max_gl > 255 else np.uint8
        _log(
            f"generateMasks: gray_levels desde organ_config.jsonc = "
            f"{sorted(set(g for g in gray_levels if g > 0))}"
        )

        # --- Geometria de salida (la del volumen de referencia) ---
        ref_arr = np.asarray(slicer.util.arrayFromVolume(ref_node))
        if ref_arr.ndim != 3:
            raise ValueError("El volumen de referencia debe ser 3D.")
        ref_shape = ref_arr.shape
        combined = np.zeros(ref_shape, dtype=combined_dtype)
        counts = np.zeros(ref_shape, dtype=np.uint16)

        # API disponible segun version de Slicer
        has_bin_api = hasattr(slicer.util, "arrayFromSegmentBinaryLabelmap")
        has_lbl_api = hasattr(slicer.util, "arrayFromSegmentLabelmap")
        if not has_bin_api and not has_lbl_api:
            raise RuntimeError(
                "Esta version de Slicer no soporta "
                "arrayFromSegmentBinaryLabelmap. Actualiza a Slicer 5.x."
            )
        if not has_bin_api:
            _log(
                "generateMasks: arrayFromSegmentBinaryLabelmap no disponible, "
                "se usa arrayFromSegmentLabelmap (> 0)"
            )

        rows = []
        individual_nodes = []
        sin_config = []
        masks = []
        # Entradas para escribir el combinado al final: (mask, gl, seg_name, nvox, is_body).
        # BODY se escribe al final y como resta de los demas segmentos.
        seg_combined_entries = []

        for idx, (seg_id, seg_name) in enumerate(seg_items):
            cfg_row = _matchConfigRow(seg_name)
            if cfg_row is None:
                sin_config.append(seg_name)
                _log(f"generateMasks: segmento '{seg_name}' sin ROI en la config")
                _log(
                    f"generateMasks: segmento '{seg_name}' (SIN CONFIG) no se incluye en el combinado"
                )

            # Extraer el array binario del segmento en la geometria de referencia
            arr = None
            try:
                if has_bin_api:
                    arr = slicer.util.arrayFromSegmentBinaryLabelmap(
                        structureNode, seg_id, ref_node
                    )
                arr = np.asarray(arr)
            except Exception:
                _log(
                    f"generateMasks: arrayFromSegmentBinaryLabelmap fallo para '{seg_name}', "
                    "se intenta arrayFromSegmentLabelmap:\n" + traceback.format_exc()
                )
                try:
                    arr = slicer.util.arrayFromSegmentLabelmap(
                        structureNode, seg_id, ref_node
                    )
                    arr = np.asarray(arr)
                except Exception:
                    _log(
                        f"generateMasks: fallo extrayendo el segmento '{seg_name}':\n"
                        + traceback.format_exc()
                    )
                    rows.append(
                        {
                            "name": seg_name,
                            "gray_level": 0,
                            "nmesh": int(cfg_row.get("nmesh", -1)) if cfg_row else -1,
                            "nvox": 0,
                            "estado": "ERROR",
                        }
                    )
                    continue

            # Todos los arrays deben compartir la geometria del volumen de referencia
            if arr.ndim != 3 or arr.shape != ref_shape:
                got = tuple(arr.shape)
                _log(
                    f"generateMasks: segmento '{seg_name}' con geometria {got} != {ref_shape}, se omite"
                )
                rows.append(
                    {
                        "name": seg_name,
                        "gray_level": 0,
                        "nmesh": int(cfg_row.get("nmesh", -1)) if cfg_row else -1,
                        "nvox": 0,
                        "estado": "GEOMETRIA",
                    }
                )
                continue

            mask = arr > 0
            nvox = int(np.count_nonzero(mask))

            # Valor en el labelmap combinado: SIEMPRE desde organ_config.jsonc.
            # Un segmento sin ROI en la config no recibe gray_level inventado.
            if cfg_row is not None:
                gl = int(cfg_row.get("gray_level", 0))
                if gl > 0:
                    estado = "OK" if nvox > 0 else "VACIA"
                else:
                    estado = "SIN GRAY_LEVEL"
            else:
                gl = 0
                estado = "SIN CONFIG"

            is_body = ("BODY" in _normName(seg_name)) or (
                "CUERPO" in _normName(seg_name)
            )
            if createCombined and nvox > 0 and cfg_row is not None and gl > 0:
                seg_combined_entries.append((mask, gl, seg_name, nvox, is_body))
            counts += mask.astype(np.uint16)
            seg_color = structureNode.GetSegmentation().GetSegment(seg_id).GetColor()
            masks.append((mask, gl, seg_name, seg_color))

            if createIndividual and (nvox > 0 or not skipEmpty):
                safe_name = re.sub(r"\s+", "_", str(seg_name).strip())
                safe_name = re.sub(r"[^0-9A-Za-z_.-]", "_", safe_name) or f"SEG_{idx}"
                seg_color = (
                    structureNode.GetSegmentation().GetSegment(seg_id).GetColor()
                )
                # Usar el transform de ref_node si existe, sino el de structureNode
                transformID = None
                if ref_node and ref_node.GetTransformNodeID():
                    transformID = ref_node.GetTransformNodeID()
                elif structureNode and structureNode.GetTransformNodeID():
                    transformID = structureNode.GetTransformNodeID()
                node = self._createLabelmapVolume(
                    mask.astype(np.uint8),
                    ref_node,
                    f"{_PREFIX_INDIVIDUAL}{safe_name}",
                    pixel_colors={1: seg_color},
                    transformNodeID=transformID,
                )
                individual_nodes.append(node)

            rows.append(
                {
                    "name": seg_name,
                    "gray_level": gl,
                    "nmesh": int(cfg_row.get("nmesh", -1))
                    if cfg_row is not None
                    else -1,
                    "nvox": nvox,
                    "estado": estado,
                }
            )

        combined_node = None
        combined_gray_node = None
        if createCombined:
            body_entries = [e for e in seg_combined_entries if e[4]]
            non_body_entries = [e for e in seg_combined_entries if not e[4]]
            if body_entries:
                # BODY se escribe AL FINAL y como resta de los demas segmentos:
                # solo los voxeles de BODY no cubiertos por ningun otro segmento.
                union_others = np.zeros(ref_shape, dtype=bool)
                for mask, gl, seg_name, nvox, _is_body in non_body_entries:
                    union_others |= mask
                # Primero los no-body; grandes primero para que los pequenos sobrelayen.
                non_body_entries.sort(key=lambda e: e[3], reverse=True)
                for mask, gl, seg_name, nvox, _is_body in non_body_entries:
                    combined[mask] = gl
                # BODY al final, recortado a lo no cubierto por los demas.
                for mask, gl, seg_name, nvox, _is_body in body_entries:
                    body_effective = mask & ~union_others
                    combined[body_effective] = gl
                    _log(
                        f"generateMasks: segmento BODY '{seg_name}' escrito como resta de los demas "
                        f"({int(np.count_nonzero(body_effective))} de {nvox} voxeles)"
                    )
            else:
                # Sin segmento BODY: comportamiento original, escribir en el orden del bucle.
                for mask, gl, seg_name, nvox, _is_body in seg_combined_entries:
                    combined[mask] = gl
            combined_colors = {}
            for mask_arr, gl, seg_name, seg_color in masks:
                combined_colors[gl] = seg_color
            # Calcular transformNodeID: usar ref_node si tiene transform, sino structureNode
            transformID = None
            if ref_node and ref_node.GetTransformNodeID():
                transformID = ref_node.GetTransformNodeID()
            elif structureNode and structureNode.GetTransformNodeID():
                transformID = structureNode.GetTransformNodeID()
            combined_node = self._createLabelmapVolume(
                combined,
                ref_node,
                _NAME_COMBINED,
                pixel_colors=combined_colors,
                transformNodeID=transformID,
            )
            # Crear nodo gris combinado para voxelization
            combined_gray_node = self._createLabelmapVolume(
                combined,
                ref_node,
                _NAME_COMBINED_GRAY,
                grayscale=True,
                transformNodeID=transformID,
            )

        # --- Verificacion de superposicion: en segmentacion SI puede haber solape ---
        overlap_voxels = int(np.count_nonzero(counts > 1))
        overlap_rois = []
        if overlap_voxels > 0:
            overlap_region = counts > 1
            for mask, gl, seg_name, seg_color in masks:
                if int(np.count_nonzero(mask & overlap_region)) > 0:
                    overlap_rois.append(f"{seg_name} (valor {gl})")

        # --- Valores repetidos en el labelmap combinado ---
        by_gl = {}
        for row in rows:
            if row["gray_level"] > 0:
                by_gl.setdefault(row["gray_level"], []).append(row["name"])
        duplicates = [
            (gl, names) for gl, names in sorted(by_gl.items()) if len(names) > 1
        ]
        for gl, names in duplicates:
            _log(f"generateMasks: CONFLICTO valor {gl} compartido por {names}")

        empty_rois = [
            row["name"]
            for row in rows
            if row["nvox"] == 0
            and row.get("estado") in ("VACIA", "SIN CONFIG", "SIN GRAY_LEVEL")
        ]

        result = {
            "rows": rows,
            "overlap_voxels": overlap_voxels,
            "overlap_rois": overlap_rois,
            "duplicates": duplicates,
            "unmapped_values": [],
            "empty_rois": empty_rois,
            "combined_node": combined_node,
            "combined_gray_node": combined_gray_node,
            "individual_nodes": individual_nodes,
            "total_voxels": int(ref_arr.size),
            "volume_name": structureNode.GetName(),
            "dims": f"{ref_shape[2]} x {ref_shape[1]} x {ref_shape[0]}",
            "n_rois_config": len(config),
            "skipped_gl0": [],
            "mode": "segmentation",
            "sin_config": sin_config,
            "n_segments": n_seg,
        }

        _log(
            f"generateMasks[segmentation]: segmentos={n_seg}, filas={len(rows)}, "
            f"superposicion={overlap_voxels}, sin_config={len(sin_config)}, "
            f"vacios={len(empty_rois)}, ref='{ref_node.GetName()}'"
        )
        if individual_nodes:
            self._groupMaskNodesInFolder(individual_nodes)
        return result

    def _roisInOverlap(self, arr_int, rows):
        """Devuelve las ROIs que aportan voxeles a la region superpuesta (solo si hay conflicto)."""
        gls = [row["gray_level"] for row in rows]
        if not gls:
            return []
        counts = np.zeros(arr_int.shape, dtype=np.uint16)
        masks = []
        for row in rows:
            mask = arr_int == row["gray_level"]
            counts += mask.astype(np.uint16)
            masks.append(mask)
        overlap = counts > 1
        names = []
        for row, mask in zip(rows, masks):
            if int(np.count_nonzero(mask & overlap)) > 0:
                names.append(f"{row['name']} (gray_level {row['gray_level']})")
        return names

    def checkOverlap(self, structureNode, config, referenceNode=None):
        """Verifica la superposicion SIN crear volumenes (solo lectura)."""
        if structureNode is None:
            raise ValueError("No hay estructura seleccionada.")

        # La validacion 3D solo aplica al modo volumen; un nodo de segmentacion no es un array
        is_seg = False
        try:
            is_seg = bool(structureNode.IsA("vtkMRMLSegmentationNode"))
        except Exception:
            is_seg = False
        if not is_seg:
            arr = slicer.util.arrayFromVolume(structureNode)
            if arr is None or np.asarray(arr).ndim != 3:
                raise ValueError("El volumen de estructura debe ser 3D.")

        # No se crean nodos: createIndividual=False y createCombined=False
        result = self.generateMasks(
            structureNode,
            config,
            createIndividual=False,
            createCombined=False,
            skipEmpty=False,
            referenceNode=referenceNode,
        )
        result["created_nodes"] = 0
        return result

    # ------------------------------------------------------------------
    # Mascaras: nodos labelmap
    # ------------------------------------------------------------------
    def _groupMaskNodesInFolder(self, nodes, folderName=None):
        """Agrupa nodos de mascaras en una carpeta del Subject Hierarchy (Data module)."""
        if folderName is None:
            folderName = _NAME_MASKS_FOLDER
        nodes = [n for n in (nodes or []) if n is not None]
        if not nodes:
            return None
        try:
            shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(
                slicer.mrmlScene
            )
            if shNode is None:
                return None
            folderId = shNode.GetItemByName(folderName)
            if folderId and shNode.GetItemLevel(folderId) != "Folder":
                folderId = 0
            if not folderId:
                folderId = shNode.CreateFolderItem(shNode.GetSceneItemID(), folderName)
            for node in nodes:
                itemId = shNode.GetItemByDataNode(node)
                if itemId:
                    shNode.SetItemParent(itemId, folderId)
                else:
                    shNode.CreateItem(folderId, node)
            shNode.SetItemExpanded(folderId, True)
            _log(
                f"Mascaras agrupadas en carpeta Subject Hierarchy: '{folderName}' ({len(nodes)} nodos)"
            )
            return folderId
        except Exception as e:
            _log(f"No se pudo agrupar las mascaras en carpeta '{folderName}': {e}")
            return None

    def _get_distinct_color(self, index, total):
        """Genera un color RGB distintivo (0-1) cíclico para distinguir ROIs."""
        import colorsys

        hue = (index * 0.61803398875) % 1.0  # golden ratio para buena distribucion
        r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 0.9)
        return (r, g, b)

    def _removeOldMaskNodes(self, keepNodes=()):
        """Borra los labelmaps generados por el modulo en una corrida anterior.

        Evita nodos duplicados tipo 'Mascaras_Combinadas_Gris_1', '_2' que hacian
        que el tab Voxelizado siguiera leyendo una mascara vieja (con gray_level
        inventados). Nunca borra los nodos pasados en `keepNodes` (entrada/CT).
        """
        keep = set()
        for n in keepNodes:
            if n is not None:
                keep.add(n.GetID())
        removed = []
        try:
            for node in list(
                slicer.mrmlScene.GetNodesByClass("vtkMRMLLabelMapVolumeNode")
            ):
                if node.GetID() in keep:
                    continue
                name = node.GetName() or ""
                if (
                    name.startswith(_PREFIX_INDIVIDUAL)
                    or name == _NAME_COMBINED
                    or name == _NAME_COMBINED_GRAY
                    or name.startswith(_NAME_COMBINED + "_")
                    or name.startswith(_NAME_COMBINED_GRAY + "_")
                ):
                    removed.append(name)
                    slicer.mrmlScene.RemoveNode(node)
        except Exception as e:
            _log(f"_removeOldMaskNodes: error: {e}")
        if removed:
            _log(
                f"_removeOldMaskNodes: eliminados {len(removed)} nodos previos: {removed}"
            )
        return removed

    def _createLabelmapVolume(
        self,
        mask_array,
        referenceNode,
        name,
        pixel_colors=None,
        transformNodeID=None,
        grayscale=False,
    ):
        """Crea un labelmap con la MISMA geometria que el volumen de referencia.

        Preserva el dtype entero de entrada (uint8 o uint16) para no truncar
        gray_level > 255. Solo fuerza uint8 si el array no es entero.

        `pixel_colors`: dict {valor: (r, g, b, [a])} para asignar colores
        distintivos a cada valor del labelmap.
        `transformNodeID`: ID del nodo transform al que pertenece el volumen de referencia.
            Si se proporciona, el labelmap se coloca bajo ese mismo nodo transform
            en lugar de multiplicar la matriz IJKToRAS.
        `grayscale`: si True, crea una tabla de colores en escala de grises en lugar
            de colores distintivos por ROI.
        """
        mask_array = np.asarray(mask_array)
        if mask_array.dtype.kind not in ("u", "i"):
            _log(
                f"_createLabelmapVolume: dtype {mask_array.dtype} no entero para "
                f"'{name}', se convierte a uint8"
            )
            mask_array = mask_array.astype(np.uint8)
        if mask_array.ndim == 2:
            _log(
                f"_createLabelmapVolume: array 2D para '{name}', se agrega una dimension"
            )
            mask_array = mask_array.reshape((1,) + mask_array.shape)
        mask_array = np.ascontiguousarray(mask_array)

        if referenceNode is None:
            raise ValueError(
                f"_createLabelmapVolume: referenceNode es None para '{name}'"
            )
        node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", name)
        slicer.util.updateVolumeFromArray(node, mask_array)
        node.SetOrigin(referenceNode.GetOrigin())
        node.SetSpacing(referenceNode.GetSpacing())
        # GetIJKToRASMatrix requiere un vtkMatrix4x4 para llenar (no devuelve la matriz)
        mat = vtk.vtkMatrix4x4()
        referenceNode.GetIJKToRASMatrix(mat)
        # Si se proporcionó un nodo transform explícito, usar SetAndObserveTransformNodeID
        # en lugar de multiplicar la matriz IJK->RAS, para que el labelmap quede alineado
        # con la jerarquía de transforms del volumen de referencia/Ct.
        if transformNodeID:
            # Mantener la geometria IJK->RAS del volumen de referencia (sin transform)
            # y luego colocar el labelmap bajo el mismo nodo transform que la CT/estructura
            node.SetIJKToRASMatrix(mat)
            node.SetAndObserveTransformNodeID(transformNodeID)
        else:
            # Comportamiento anterior: intentar aplicar transform multipliéndolo a IJKToRAS
            transformNodeID = referenceNode.GetTransformNodeID()
            transformNode = None
            if transformNodeID:
                transformNode = slicer.mrmlScene.GetNodeByID(transformNodeID)
            if transformNode:
                try:
                    tmat = vtk.vtkMatrix4x4()
                    transformNode.GetMatrixTransformToWorld(tmat)
                    finalMat = vtk.vtkMatrix4x4()
                    vtk.vtkMatrix4x4.Multiply4x4(tmat, mat, finalMat)
                    node.SetIJKToRASMatrix(finalMat)
                    _log(
                        f"_createLabelmapVolume: aplicado transform de '{referenceNode.GetName()}' a '{name}'"
                    )
                except Exception as e:
                    _log(
                        f"_createLabelmapVolume: no se pudo aplicar transform de '{referenceNode.GetName()}': {e}"
                    )
                    node.SetIJKToRASMatrix(mat)
            else:
                # No hay transform, usar la matriz IJKToRAS original
                node.SetIJKToRASMatrix(mat)
        node.CreateDefaultDisplayNodes()
        display = node.GetDisplayNode()
        if display:
            if pixel_colors and not grayscale:
                # Crear un color table node personalizado con los colores asignados por valor.
                # SetPixelValueColor no es confiable en Slicer 5.8.1 para labelmaps.
                import colorsys as _colorsys

                colorNode = slicer.mrmlScene.AddNewNodeByClass(
                    "vtkMRMLColorTableNode", name + "_ColorTable"
                )
                colorNode.SetTypeToUser()
                # Asignar suficientes slots antes de SetColor() (exigido por vtkMRMLColorTableNode)
                max_val = max(pixel_colors.keys()) if pixel_colors else 0
                colorNode.SetNumberOfColors(max_val + 1)
                # Valor 0 = aire/fondo (transparente o negro)
                colorNode.SetColor(0, 0.0, 0.0, 0.0, 0.0)
                for valor, color in pixel_colors.items():
                    r, g, b = color[0], color[1], color[2]
                    a = color[3] if len(color) > 3 else 1.0
                    colorNode.SetColor(int(valor), r, g, b, a)
                display.SetAndObserveColorNodeID(colorNode.GetID())
            elif grayscale:
                # Tabla de colores en escala de grises: valor 0 transparente, 1..max_val en grises
                max_val = int(np.max(mask_array))
                if max_val < 1:
                    max_val = 1
                colorNode = slicer.mrmlScene.AddNewNodeByClass(
                    "vtkMRMLColorTableNode", name + "_ColorTable_Gray"
                )
                colorNode.SetTypeToUser()
                colorNode.SetNumberOfColors(max_val + 1)
                # Valor 0 = aire/fondo (transparente)
                colorNode.SetColor(0, 0.0, 0.0, 0.0, 0.0)
                for v in range(1, max_val + 1):
                    gray_val = v / max_val
                    colorNode.SetColor(v, gray_val, gray_val, gray_val, 1.0)
                display.SetAndObserveColorNodeID(colorNode.GetID())
                # Establecer rango escalar para que la escala de grises sea visible
                display.SetScalarRange(0, max_val)
            else:
                display.SetAndObserveColorNodeID("vtkMRMLColorTableNodeLabels")
        return node

    def cleanupMasks(self):
        """Elimina de la escena los labelmaps generados por este modulo."""
        import traceback

        removed = 0
        try:
            collection = slicer.mrmlScene.GetNodesByClass("vtkMRMLLabelMapVolumeNode")
            nodes = [
                collection.GetItemAsObject(i)
                for i in range(collection.GetNumberOfItems())
            ]
            for node in nodes:
                if node is None:
                    continue
                name = node.GetName() or ""
                if (
                    name.startswith(_PREFIX_INDIVIDUAL)
                    or name == _NAME_COMBINED
                    or name == _NAME_COMBINED_GRAY
                ):
                    slicer.mrmlScene.RemoveNode(node)
                    removed += 1
        except Exception:
            _log("Error en cleanupMasks:\n" + traceback.format_exc())
        _log(f"cleanupMasks: {removed} nodos eliminados")
        return removed

    # ------------------------------------------------------------------
    # Centrado del BODY en el origen RAS
    # ------------------------------------------------------------------
    def _findBodySegment(self, segmentationNode):
        """Devuelve (segmentId, segmentName) del segmento BODY/CUERPO.

        Reusa el criterio de la generacion de mascaras: el nombre normalizado
        (`_normName`) contiene BODY o CUERPO. Se prefiere la coincidencia
        exacta sobre la coincidencia por subcadena.
        """
        if segmentationNode is None:
            return (None, None)
        try:
            if not segmentationNode.IsA("vtkMRMLSegmentationNode"):
                return (None, None)
            segmentation = segmentationNode.GetSegmentation()
            if segmentation is None:
                return (None, None)
            n_seg = segmentation.GetNumberOfSegments()
        except Exception:
            _log("_findBodySegment: no se pudo leer la segmentacion")
            return (None, None)

        fallback = (None, None)
        for i in range(n_seg):
            try:
                seg_id = segmentation.GetNthSegmentID(i)
                seg = segmentation.GetSegment(seg_id)
                seg_name = seg.GetName() if seg else seg_id
            except Exception:
                continue
            norm = _normName(seg_name)
            if norm in ("BODY", "CUERPO"):
                return (seg_id, seg_name)
            if ("BODY" in norm) or ("CUERPO" in norm):
                if fallback[0] is None:
                    fallback = (seg_id, seg_name)
        return fallback

    def _bodyTightBoundsRAS(self, array, referenceNode):
        """BBox ajustado (tight) del BODY mapeado a RAS.

        `array` es un array binario en la geometria local de `referenceNode`
        (orden k, j, i = Z, Y, X). Devuelve (bounds, center) con
        bounds = [xmin, xmax, ymin, ymax, zmin, zmax] y center = [cx, cy, cz].
        Devuelve (None, None) si el array es None, no es 3D o esta vacio.
        """
        if array is None or referenceNode is None:
            return (None, None)
        arr = np.asarray(array)
        if arr.ndim != 3:
            _log("_bodyTightBoundsRAS: el array del BODY no es 3D")
            return (None, None)
        mask = arr > 0
        if not bool(mask.any()):
            return (None, None)

        ks = np.where(mask.any(axis=(1, 2)))[0]
        js = np.where(mask.any(axis=(0, 2)))[0]
        iss = np.where(mask.any(axis=(0, 1)))[0]
        if ks.size == 0 or js.size == 0 or iss.size == 0:
            return (None, None)

        # El voxel i ocupa [i-0.5, i+0.5] en el espacio IJK: se usan las
        # caras externas del bloque de voxeles para el extent fisico real.
        i_vals = (float(iss[0]) - 0.5, float(iss[-1]) + 0.5)
        j_vals = (float(js[0]) - 0.5, float(js[-1]) + 0.5)
        k_vals = (float(ks[0]) - 0.5, float(ks[-1]) + 0.5)

        corners = np.array(
            [[i, j, k] for i in i_vals for j in j_vals for k in k_vals],
            dtype=np.float64,
        )

        mat = vtk.vtkMatrix4x4()
        referenceNode.GetIJKToRASMatrix(mat)
        m_np = np.asarray(slicer.util.arrayFromVTKMatrix(mat), dtype=np.float64)
        m_np = m_np.reshape((4, 4))
        homo = np.vstack([corners.T, np.ones((1, corners.shape[0]))])
        ras = m_np @ homo

        bounds = [
            float(ras[0].min()),
            float(ras[0].max()),
            float(ras[1].min()),
            float(ras[1].max()),
            float(ras[2].min()),
            float(ras[2].max()),
        ]
        center = [
            (bounds[0] + bounds[1]) / 2.0,
            (bounds[2] + bounds[3]) / 2.0,
            (bounds[4] + bounds[5]) / 2.0,
        ]
        return (bounds, center)

    def _findBodyLabelmap(self):
        """Labelmap BODY del modulo (Mascara_BODY o nombre que contenga BODY)."""
        import traceback

        try:
            collection = slicer.mrmlScene.GetNodesByClass(
                "vtkMRMLLabelMapVolumeNode"
            )
            candidates = []
            for i in range(collection.GetNumberOfItems()):
                node = collection.GetItemAsObject(i)
                if node is None:
                    continue
                norm = _normName(node.GetName())
                if ("BODY" in norm) or ("CUERPO" in norm):
                    candidates.append(node)
            if not candidates:
                return None
            exact = f"{_PREFIX_INDIVIDUAL}BODY"
            for node in candidates:
                if (node.GetName() or "") == exact:
                    return node
            return candidates[0]
        except Exception:
            _log("_findBodyLabelmap: error:\n" + traceback.format_exc())
            return None

    def _findReferenceVolume(self):
        """Primer volumen escalar 3D de la escena que puede servir de geometria."""
        try:
            collection = slicer.mrmlScene.GetNodesByClass(
                "vtkMRMLScalarVolumeNode"
            )
            for i in range(collection.GetNumberOfItems()):
                node = collection.GetItemAsObject(i)
                if node is None:
                    continue
                name = node.GetName() or ""
                if any(x in name for x in ["Kerma", "Gamma", "Dosis", "DVH", "Temp"]):
                    continue
                if node.IsA("vtkMRMLLabelMapVolumeNode"):
                    continue
                if name.startswith(_PREFIX_INDIVIDUAL) or name.startswith(
                    _NAME_COMBINED
                ):
                    continue
                try:
                    arr = slicer.util.arrayFromVolume(node)
                    if arr is not None and np.asarray(arr).ndim == 3:
                        return node
                except Exception:
                    continue
        except Exception:
            pass
        return None

    def _sameIJKToRAS(self, nodeA, nodeB, tol=1e-4):
        """True si dos nodos comparten la misma matriz IJK->RAS local."""
        if nodeA is None or nodeB is None:
            return False
        try:
            ma = vtk.vtkMatrix4x4()
            mb = vtk.vtkMatrix4x4()
            nodeA.GetIJKToRASMatrix(ma)
            nodeB.GetIJKToRASMatrix(mb)
            a = np.asarray(
                slicer.util.arrayFromVTKMatrix(ma), dtype=np.float64
            ).reshape((4, 4))
            b = np.asarray(
                slicer.util.arrayFromVTKMatrix(mb), dtype=np.float64
            ).reshape((4, 4))
            return bool(np.allclose(a, b, atol=tol))
        except Exception:
            return False

    def _transformChain(self, node):
        """Cadena de transforms de un nodo: [padre directo, ..., raiz]."""
        chain = []
        if node is None:
            return chain
        try:
            tid = node.GetTransformNodeID()
        except Exception:
            return chain
        seen = set()
        while tid:
            if tid in seen:
                break
            seen.add(tid)
            tnode = slicer.mrmlScene.GetNodeByID(tid)
            if tnode is None:
                break
            chain.append(tnode)
            try:
                tid = tnode.GetTransformNodeID()
            except Exception:
                break
        return chain

    def _chainMatrix(self, chain, overrideNode=None, overrideMat=None):
        """Matriz 4x4 (numpy) de una cadena de transforms [directo, ..., raiz].

        Devuelve M_raiz * ... * M_directo (la geometria local se aplica
        despues, por la izquierda, con este factor).
        `overrideNode`/`overrideMat` permiten sustituir la matriz de un
        transform concreto de la cadena.
        """
        m = np.eye(4, dtype=np.float64)
        for tnode in reversed(chain):
            if overrideNode is not None and tnode.GetID() == overrideNode.GetID():
                tm_np = np.asarray(overrideMat, dtype=np.float64).reshape((4, 4))
            else:
                tm = vtk.vtkMatrix4x4()
                try:
                    tnode.GetMatrixTransformToParent(tm)
                except Exception:
                    continue
                tm_np = np.asarray(
                    slicer.util.arrayFromVTKMatrix(tm), dtype=np.float64
                ).reshape((4, 4))
            m = m @ tm_np
        return m

    def _boundsFromLocal(self, bounds_local, mat_np):
        """AABB en mundo de un AABB local mapeado por `mat_np` (4x4)."""
        corners = np.array(
            [
                [x, y, z, 1.0]
                for x in (bounds_local[0], bounds_local[1])
                for y in (bounds_local[2], bounds_local[3])
                for z in (bounds_local[4], bounds_local[5])
            ],
            dtype=np.float64,
        ).T
        ras = np.asarray(mat_np, dtype=np.float64).reshape((4, 4)) @ corners
        return [
            float(ras[0].min()),
            float(ras[0].max()),
            float(ras[1].min()),
            float(ras[1].max()),
            float(ras[2].min()),
            float(ras[2].max()),
        ]

    @staticmethod
    def _volumeDimensions(node):
        """(nx, ny, nz) del volumen; None si no se puede determinar.

        Usa primero `GetImageData().GetDimensions()` y, como fallback, la
        forma del array (`slicer.util.arrayFromVolume`), que viene en orden
        (nz, ny, nx).
        """
        if node is None:
            return None
        try:
            img = node.GetImageData()
            if img is not None:
                dims = img.GetDimensions()
                if dims and all(int(v) > 0 for v in dims):
                    return (int(dims[0]), int(dims[1]), int(dims[2]))
        except Exception:
            pass
        try:
            arr = slicer.util.arrayFromVolume(node)
            if arr is not None:
                shape = np.asarray(arr).shape
                if len(shape) == 3:
                    return (int(shape[2]), int(shape[1]), int(shape[0]))
        except Exception:
            pass
        return None

    def _isAxisAlignedIJKToRAS(self, node, tol=1e-6):
        """True si la matriz IJK->RAS del nodo es diagonal (ejes alineados).

        La formula de centrado asume que el eje indice coincide con el eje
        RAS. Si la matriz no es diagonal se loguea un WARNING explicito y se
        devuelve False (el llamador sigue adelante con `target_local`).
        """
        try:
            mat = vtk.vtkMatrix4x4()
            node.GetIJKToRASMatrix(mat)
            m = np.asarray(
                slicer.util.arrayFromVTKMatrix(mat), dtype=np.float64
            ).reshape((4, 4))
        except Exception:
            _log(
                "centerBody: WARNING no se pudo leer la matriz IJK->RAS; se "
                "asume que los ejes indice/RAS estan alineados"
            )
            return False
        lin = m[:3, :3]
        offdiag = lin - np.diag(np.diag(lin))
        if not np.allclose(offdiag, 0.0, atol=tol):
            _log(
                "centerBody: WARNING la matriz IJK->RAS de la imagen de "
                "referencia no es diagonal (el eje indice no coincide con "
                "RAS); la formula de centrado asume ejes alineados y se "
                "aplica de todos modos"
            )
            return False
        return True

    def centerBody(self, referenceNode=None, structureNode=None):
        """Centra la imagen de referencia en el origen RAS (0,0,0).

        Crea (o reutiliza) un `vtkMRMLLinearTransformNode` con una traslacion
        pura y lo encadena como hijo del transform padre existente (nunca se
        pisa el padre a ciegas).

        El offset NO sale del tight bbox del BODY sino de las DIMENSIONES y el
        SPACING de la imagen de referencia:

            shift_x = N_x * _CENTER_X_COEF
            shift_y = N_y * dy / 2
            shift_z = N_z * dz / 2
            target_local = origin + shift

        donde `origin` es `referenceNode.GetOrigin()` (frame RAS local de la
        imagen). `target_local` es el punto que debe quedar en el origen RAS.
        Si la imagen tiene un transform padre, la traslacion se compensa con
        la inversa de la parte lineal de ese padre.

        El BODY (segmento BODY/CUERPO o labelmap `Mascara_BODY`) se sigue
        detectando para incluirlo entre los targets y reportar sus bounds
        finales, pero NO interviene en el calculo del offset. Si no se
        encuentra, el centrado de la imagen sigue funcionando igual.

        Devuelve un dict con 'offset', 'bounds', 'image_bounds', 'dims',
        'spacing', 'n_nodes', 'transform_name' y 'segment_name', o None si no
        se pudo determinar la geometria de la imagen de referencia.
        """
        import traceback

        try:
            body_array = None
            body_geom_node = None
            body_node = None
            segment_name = None
            source = None

            # --- 0) Imagen de referencia: define dimensiones, spacing, origin
            #         y la cadena de transforms canonica ---
            ref = referenceNode
            if ref is None:
                ref = self._findReferenceVolume()
                if ref is not None:
                    _log(
                        "centerBody: no hay imagen de referencia seleccionada; "
                        f"se usa '{ref.GetName()}'"
                    )
            if ref is None:
                _log(
                    "centerBody: no hay imagen de referencia en la escena; "
                    "no se puede calcular el centrado"
                )
                return None

            dims = self._volumeDimensions(ref)
            if dims is None:
                _log(
                    "centerBody: no se pudieron obtener las dimensiones de la "
                    f"imagen de referencia '{ref.GetName()}'"
                )
                return None
            spacing = ref.GetSpacing()
            origin = ref.GetOrigin()
            nx, ny, nz = (int(dims[0]), int(dims[1]), int(dims[2]))
            dx, dy, dz = (
                float(spacing[0]),
                float(spacing[1]),
                float(spacing[2]),
            )

            # La formula asume que el eje indice coincide con el eje RAS.
            self._isAxisAlignedIJKToRAS(ref)

            # --- 1) BODY desde la segmentacion (solo para targets y bounds) ---
            if structureNode is not None:
                seg_id, seg_name = self._findBodySegment(structureNode)
                if seg_id is not None:
                    arr = None
                    try:
                        arr = slicer.util.arrayFromSegmentBinaryLabelmap(
                            structureNode, seg_id, ref
                        )
                    except Exception:
                        arr = None
                    if arr is None:
                        try:
                            arr = slicer.util.arrayFromSegmentLabelmap(
                                structureNode, seg_id, ref
                            )
                        except Exception:
                            arr = None
                    if arr is not None:
                        body_array = np.asarray(arr)
                        body_geom_node = ref
                        body_node = ref
                        segment_name = seg_name
                        source = "segmentacion"

            # --- 2) fallback: labelmap BODY del modulo ---
            if body_array is None:
                body_lm = self._findBodyLabelmap()
                if body_lm is not None:
                    arr = slicer.util.arrayFromVolume(body_lm)
                    if arr is not None:
                        body_array = np.asarray(arr)
                        body_geom_node = body_lm
                        body_node = body_lm
                        segment_name = body_lm.GetName()
                        source = "labelmap"
                        # Si la CT de referencia comparte la misma geometria
                        # local, se usa como ancla: la CT define el transform
                        # padre que se preserva al encadenar el centrado.
                        if self._sameIJKToRAS(body_lm, ref):
                            body_geom_node = ref

            if body_array is None:
                _log(
                    "centerBody: no se encontro el BODY; se centra igual la "
                    "imagen de referencia"
                )

            # --- 3) Offset a partir de dimensiones / spacing / origin ---
            shift_x = nx * _CENTER_X_COEF
            shift_y = ny * dy / 2.0
            shift_z = nz * dz / 2.0
            target_local = (
                float(origin[0]) + shift_x,
                float(origin[1]) + shift_y,
                float(origin[2]) + shift_z,
            )

            # --- Nodo transform de centrado (reutilizable) ---
            transformNode = slicer.mrmlScene.GetFirstNodeByName(_NAME_CENTER_TRANSFORM)
            if transformNode is not None and not transformNode.IsA(
                "vtkMRMLLinearTransformNode"
            ):
                slicer.mrmlScene.RemoveNode(transformNode)
                transformNode = None
            if transformNode is None:
                transformNode = slicer.mrmlScene.AddNewNodeByClass(
                    "vtkMRMLLinearTransformNode", _NAME_CENTER_TRANSFORM
                )

            chain = self._transformChain(ref)
            ct_index = None
            for i, tnode in enumerate(chain):
                if tnode.GetID() == transformNode.GetID():
                    ct_index = i
                    break

            if ct_index is None:
                # El transform se inserta entre el nodo y su padre directo:
                # el padre canonico se conserva como padre del transform.
                canonical_parent = chain[0] if chain else None
                above = chain
                below = []
            else:
                canonical_parent = (
                    chain[ct_index + 1] if ct_index + 1 < len(chain) else None
                )
                above = chain[ct_index + 1 :]
                below = chain[:ct_index]

            # Punto objetivo en mundo, sin el efecto del transform de centrado
            A = self._chainMatrix(above)
            B = self._chainMatrix(below)
            target_vec = np.array(
                [target_local[0], target_local[1], target_local[2], 1.0],
                dtype=np.float64,
            )
            target_world = A @ (B @ target_vec)
            try:
                inv_lin = np.linalg.inv(A[:3, :3])
            except np.linalg.LinAlgError:
                _log(
                    "centerBody: la parte lineal del transform padre no es "
                    "invertible; no se aplica el centrado"
                )
                return None
            # Traslacion en el marco local del padre canonico
            offset = -(inv_lin @ target_world[:3])

            mat = vtk.vtkMatrix4x4()
            mat.Identity()
            mat.SetElement(0, 3, float(offset[0]))
            mat.SetElement(1, 3, float(offset[1]))
            mat.SetElement(2, 3, float(offset[2]))

            if canonical_parent is not None:
                transformNode.SetAndObserveTransformNodeID(canonical_parent.GetID())
            else:
                transformNode.SetAndObserveTransformNodeID(None)
            transformNode.SetMatrixTransformToParent(mat)

            # --- Targets a mover ---
            targets = []
            if structureNode is not None:
                targets.append(structureNode)
            if ref is not None:
                targets.append(ref)
            body_lm = self._findBodyLabelmap()
            if body_lm is not None and source is not None:
                targets.append(body_lm)
            patient = slicer.mrmlScene.GetFirstNodeByName("Paciente")
            if patient is not None:
                targets.append(patient)
            try:
                collection = slicer.mrmlScene.GetNodesByClass(
                    "vtkMRMLLabelMapVolumeNode"
                )
                for i in range(collection.GetNumberOfItems()):
                    node = collection.GetItemAsObject(i)
                    if node is None:
                        continue
                    name = node.GetName() or ""
                    if name.startswith(_PREFIX_INDIVIDUAL) or name.startswith(
                        _NAME_COMBINED
                    ):
                        targets.append(node)
            except Exception:
                pass

            parent_id = (
                canonical_parent.GetID() if canonical_parent is not None else None
            )
            ct_id = transformNode.GetID()
            applied = 0
            skipped = []
            seen_ids = set()
            for node in targets:
                if node is None:
                    continue
                try:
                    nid = node.GetID()
                except Exception:
                    continue
                if nid in seen_ids:
                    continue
                seen_ids.add(nid)
                name = node.GetName() or nid
                # Nunca mover la geometria fija del simulador PHITS
                if name in ("Conos_Pared_Reactor", "Cono"):
                    continue
                node_chain = self._transformChain(node)
                if any(t.GetID() == ct_id for t in node_chain):
                    # Ya estaba bajo el transform de centrado: se mueve solo
                    applied += 1
                    continue
                prev_parent_id = node.GetTransformNodeID()
                if prev_parent_id != parent_id:
                    skipped.append(name)
                    continue
                if prev_parent_id:
                    try:
                        transformNode.SetNodeReferenceID(
                            f"{_NAME_CENTER_TRANSFORM}_prev_{nid}", prev_parent_id
                        )
                    except Exception:
                        pass
                node.SetAndObserveTransformNodeID(ct_id)
                applied += 1

            if skipped:
                _log(
                    "centerBody: nodos omitidos por transform padre incompatible: "
                    + ", ".join(skipped)
                )

            # --- Bounds finales (RAS mundo) post-transform ---
            mat_np = np.asarray(
                slicer.util.arrayFromVTKMatrix(mat), dtype=np.float64
            ).reshape((4, 4))

            # Imagen de referencia (post-transform)
            image_bounds_local = [
                float(origin[0]) - dx / 2.0,
                float(origin[0]) + (nx - 0.5) * dx,
                float(origin[1]) - dy / 2.0,
                float(origin[1]) + (ny - 0.5) * dy,
                float(origin[2]) - dz / 2.0,
                float(origin[2]) + (nz - 0.5) * dz,
            ]
            image_chain = self._chainMatrix(
                self._transformChain(ref),
                overrideNode=transformNode,
                overrideMat=mat_np,
            )
            image_bounds = self._boundsFromLocal(image_bounds_local, image_chain)

            # BODY (post-transform). El tight bbox se usa SOLO para reportar
            # los bounds finales: ya no participa en el calculo del offset.
            bounds_local = None
            final_bounds = None
            if body_array is not None and body_geom_node is not None:
                bounds_local, _ = self._bodyTightBoundsRAS(
                    body_array, body_geom_node
                )
                if bounds_local is not None:
                    body_chain = self._chainMatrix(
                        self._transformChain(body_node),
                        overrideNode=transformNode,
                        overrideMat=mat_np,
                    )
                    final_bounds = self._boundsFromLocal(bounds_local, body_chain)

            _log(
                f"centerBody: imagen '{ref.GetName()}' dims={dims} "
                f"spacing={(dx, dy, dz)} offset="
                f"{tuple(float(v) for v in offset)} nodos movidos={applied}"
            )

            return {
                "offset": (
                    float(offset[0]),
                    float(offset[1]),
                    float(offset[2]),
                ),
                "offset_world": (
                    float(-target_world[0]),
                    float(-target_world[1]),
                    float(-target_world[2]),
                ),
                "bounds": final_bounds,
                "bounds_local": bounds_local,
                "image_bounds": image_bounds,
                "image_bounds_local": image_bounds_local,
                "dims": (nx, ny, nz),
                "spacing": (dx, dy, dz),
                "origin": (
                    float(origin[0]),
                    float(origin[1]),
                    float(origin[2]),
                ),
                "shift": (float(shift_x), float(shift_y), float(shift_z)),
                "target_local": target_local,
                "target_world": (
                    float(target_world[0]),
                    float(target_world[1]),
                    float(target_world[2]),
                ),
                "body_found": bool(body_array is not None),
                "n_nodes": applied,
                "transform_name": _NAME_CENTER_TRANSFORM,
                "reference_name": ref.GetName() or "",
                "segment_name": segment_name or "",
                "source": source or "",
                "skipped": skipped,
            }
        except Exception:
            _log("centerBody: error:\n" + traceback.format_exc())
            return None

    def uncenterBody(self):
        """Quita el transform de centrado a los nodos que lo usan y lo elimina.

        Devuelve el numero de nodos liberados.
        """
        import traceback

        transformNode = slicer.mrmlScene.GetFirstNodeByName(_NAME_CENTER_TRANSFORM)
        if transformNode is None:
            _log("uncenterBody: no existe el transform de centrado")
            return 0

        ct_id = transformNode.GetID()
        released = 0
        try:
            nodes = slicer.mrmlScene.GetNodes()
            for i in range(nodes.GetNumberOfItems()):
                node = nodes.GetItemAsObject(i)
                if node is None:
                    continue
                try:
                    if node.GetTransformNodeID() != ct_id:
                        continue
                except Exception:
                    continue
                # Restaurar el transform padre original (guardado al centrar)
                prev_id = None
                try:
                    prev_id = transformNode.GetNodeReferenceID(
                        f"{_NAME_CENTER_TRANSFORM}_prev_{node.GetID()}"
                    )
                except Exception:
                    prev_id = None
                node.SetAndObserveTransformNodeID(prev_id if prev_id else None)
                released += 1
            slicer.mrmlScene.RemoveNode(transformNode)
        except Exception:
            _log("uncenterBody: error:\n" + traceback.format_exc())
        _log(f"uncenterBody: {released} nodos liberados")
        return released

    def _modelToWorldPolyData(self, model):
        """Devuelve el polydata del modelo transformado a world (o None)."""
        try:
            poly = model.GetPolyData()
            if poly is None or poly.GetNumberOfPoints() == 0:
                return None
            parent = model.GetParentTransformNode()
            if parent is None:
                return poly
            transform = vtk.vtkGeneralTransform()
            slicer.vtkMRMLTransformNode.GetTransformBetweenNodes(
                parent, None, transform
            )
            transformer = vtk.vtkTransformPolyDataFilter()
            transformer.SetTransform(transform)
            transformer.SetInputData(poly)
            transformer.Update()
            return transformer.GetOutput()
        except Exception as e:
            _log(f"_modelToWorldPolyData: error con '{model.GetName()}': {e}")
            return None

    def detectCollisions(self, model1, model2):
        """Detecta colisiones entre dos vtkMRMLModelNode. Retorna dict o None."""
        try:
            if model1 is None or model2 is None:
                return None
            poly1 = self._modelToWorldPolyData(model1)
            poly2 = self._modelToWorldPolyData(model2)
            if poly1 is None or poly2 is None:
                _log(
                    "detectCollisions: no se pudo obtener polydata de uno o ambos modelos"
                )
                return None

            tri1 = vtk.vtkTriangleFilter()
            tri1.SetInputData(poly1)
            tri1.Update()
            tri2 = vtk.vtkTriangleFilter()
            tri2.SetInputData(poly2)
            tri2.Update()

            detector = vtk.vtkCollisionDetectionFilter()
            detector.SetInputData(0, tri1.GetOutput())
            detector.SetInputData(1, tri2.GetOutput())
            m1 = vtk.vtkMatrix4x4()
            m2 = vtk.vtkMatrix4x4()
            detector.SetMatrix(0, m1)
            detector.SetMatrix(1, m2)
            detector.SetCollisionModeToAllContacts()
            detector.Update()

            contacts = int(detector.GetNumberOfContacts())
            _log(
                f"detectCollisions: '{model1.GetName()}' vs '{model2.GetName()}' -> {contacts} contactos"
            )
            return {"contacts": contacts, "collision": contacts > 0}
        except Exception as e:
            _log(f"detectCollisions: error: {e}")
            return None

    def _mapVoxelOrientation(self, arr, orientation):
        """Reordena el array de Slicer (k,j,i) al orden y orientacion que espera PHITS."""
        M = np.transpose(arr, (1, 2, 0))  # Slicer (k,j,i) -> (Y, X, Z)
        if orientation == "Axial":
            M = np.transpose(M, (1, 0, 2))
            M = M[:, ::-1, :]
            M = np.flip(M, axis=2)
        elif orientation == "Sagital":
            M = np.transpose(M, (2, 1, 0))
            M = M[:, ::-1, :]
            M = np.flip(M, axis=2)
        elif orientation == "Coronal":
            M = np.transpose(M, (2, 0, 1))
            M = M[:, ::-1, :]
            M = np.flip(M, axis=2)
        else:
            raise ValueError(f"Orientacion no valida: {orientation}")
        return M

    def _rle_encode(self, values):
        """Run-length encoding estilo PHITS: valor y -(L-1) cuando la corrida es L > 1."""
        arr = np.asarray(values).ravel()
        if arr.size == 0:
            return []
        change = np.flatnonzero(arr[1:] != arr[:-1]) + 1
        starts = np.concatenate(([0], change))
        ends = np.concatenate((change, [arr.size]))
        lengths = ends - starts
        vals = arr[starts]
        out = []
        for v, L in zip(vals, lengths):
            out.append(int(v))
            if L > 1:
                out.append(-(int(L) - 1))
        return out

    def resampleMaskVolume(self, maskNode, factorX, factorY, factorZ):
        """Submuestrea un labelmap (nearest-neighbour) y deja el nodo nuevo en la escena.

        Los factores deben cumplir 0 < f <= 1. Con f = 1 en los tres ejes
        devuelve el mismo nodo de entrada, sin crear nada nuevo.
        """
        try:
            if maskNode is None:
                return None
            factors = [float(factorX), float(factorY), float(factorZ)]
            for f in factors:
                if f <= 0.0 or f > 1.0:
                    _log(
                        f"resampleMaskVolume: factor invalido {f} (debe ser 0 < f <= 1)"
                    )
                    return None
            if all(abs(f - 1.0) < 1e-9 for f in factors):
                return maskNode

            arr = slicer.util.arrayFromVolume(maskNode)
            if arr is None:
                _log("resampleMaskVolume: arrayFromVolume devolvio None")
                return None
            if arr.ndim != 3:
                _log(
                    f"resampleMaskVolume: array con {arr.ndim} dimensiones, se esperaba 3"
                )
                return None

            # El array de Slicer es (k, j, i) = (Z, Y, X)
            axes_factors = [factors[2], factors[1], factors[0]]
            spacing = maskNode.GetSpacing()  # (sx, sy, sz)
            axes_spacing = [spacing[2], spacing[1], spacing[0]]  # (sz, sy, sx)

            new_shape = [
                max(1, int(round(arr.shape[a] * axes_factors[a]))) for a in range(3)
            ]
            idx = []
            for a in range(3):
                f = axes_factors[a]
                n = np.arange(new_shape[a], dtype=np.float64)
                i = np.rint((n + 0.5) / f - 0.5).astype(np.int64)
                np.clip(i, 0, arr.shape[a] - 1, out=i)
                idx.append(i)
            resampled = arr[np.ix_(idx[0], idx[1], idx[2])]

            name = (
                f"{maskNode.GetName()}_resample_"
                f"{factors[0]:g}x{factors[1]:g}x{factors[2]:g}"
            )
            old = slicer.mrmlScene.GetFirstNodeByName(name)
            if old is not None:
                slicer.mrmlScene.RemoveNode(old)

            node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLabelMapVolumeNode", name)
            slicer.util.updateVolumeFromArray(node, np.ascontiguousarray(resampled))
            node.SetOrigin(maskNode.GetOrigin())

            new_spacing_zyx = [axes_spacing[a] / axes_factors[a] for a in range(3)]
            node.SetSpacing(new_spacing_zyx[2], new_spacing_zyx[1], new_spacing_zyx[0])

            # Copiar la orientacion IJK->RAS escalando cada columna por el nuevo spacing
            mat = vtk.vtkMatrix4x4()
            maskNode.GetIJKToRASMatrix(mat)
            out = vtk.vtkMatrix4x4()
            for c in range(3):
                col = [mat.GetElement(r, c) for r in range(3)]
                norm = float(np.linalg.norm(col))
                scale = (new_spacing_zyx[2 - c] / norm) if norm > 0 else 1.0
                for r in range(3):
                    out.SetElement(r, c, col[r] * scale)
            for r in range(4):
                out.SetElement(r, 3, mat.GetElement(r, 3))
            node.SetIJKToRASMatrix(out)

            tid = maskNode.GetTransformNodeID()
            if tid:
                node.SetAndObserveTransformNodeID(tid)

            node.CreateDefaultDisplayNodes()
            try:
                src_disp = maskNode.GetDisplayNode()
                dst_disp = node.GetDisplayNode()
                if src_disp is not None and dst_disp is not None:
                    dst_disp.SetAndObserveColorNodeID(src_disp.GetColorNodeID())
            except Exception:
                pass

            _log(
                f"resampleMaskVolume: '{maskNode.GetName()}' {tuple(arr.shape)} -> "
                f"{tuple(resampled.shape)} factores=({factors[0]:g},{factors[1]:g},"
                f"{factors[2]:g}) spacing={tuple(node.GetSpacing())}"
            )
            return node
        except Exception as e:
            _log(f"resampleMaskVolume: error: {e}")
            return None

    def generateVoxelDat(self, maskNode, universo_by_gray, orientation, outputPath):
        """Genera el voxel.dat de PHITS desde la mascara gris. Retorna dict o None."""
        try:
            arr = slicer.util.arrayFromVolume(maskNode)
            if arr is None:
                _log("generateVoxelDat: arrayFromVolume devolvio None")
                return None

            arr_int = np.rint(np.asarray(arr, dtype=np.float64)).astype(np.int64)
            M = self._mapVoxelOrientation(arr_int, orientation)

            universo_by_gray = {int(k): int(v) for k, v in universo_by_gray.items()}
            present = [int(v) for v in np.unique(M)]
            missing = [v for v in present if v not in universo_by_gray]
            if missing:
                _log(f"generateVoxelDat: gray_levels sin universo asignado: {missing}")
                return None

            mapped = np.zeros(M.shape, dtype=np.int32)
            for gl, uni in universo_by_gray.items():
                mapped[M == gl] = uni

            flat = mapped.flatten(order="F")
            encoded = self._rle_encode(flat)

            out_dir = os.path.dirname(outputPath)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)

            per_line = 21
            with open(outputPath, "w") as f:
                for i in range(0, len(encoded), per_line):
                    chunk = encoded[i : i + per_line]
                    f.write("     " + " ".join(str(x) for x in chunk) + " \n")

            universos = sorted(set(universo_by_gray[gl] for gl in present))
            _log(
                f"generateVoxelDat: {outputPath} shape={M.shape} orient={orientation} "
                f"voxeles={int(flat.size)} valores_rle={len(encoded)} universos={universos}"
            )
            return {
                "path": outputPath,
                "values": len(encoded),
                "voxels": int(flat.size),
                "shape": tuple(int(s) for s in M.shape),
                "orientation": orientation,
                "universos": universos,
            }
        except Exception as e:
            _log(f"generateVoxelDat: error: {e}")
            return None

    # ------------------------------------------------------------------
    # Mascaras: informe
    # ------------------------------------------------------------------
    def buildReport(self, result):
        """Construye el texto del informe de verificacion a partir del dict de resultados."""
        rows = result.get("rows", [])
        lines = []
        lines.append("=== VERIFICACION DE MASCARAS BNCT ===")
        lines.append(f"Volumen: {result.get('volume_name', '-')}")
        lines.append(f"Dimensiones (X x Y x Z): {result.get('dims', '-')}")
        lines.append(f"Voxeles totales: {result.get('total_voxels', 0)}")
        lines.append(
            f"ROIs en configuracion: {result.get('n_rois_config', len(rows))}   "
            f"evaluadas: {len(rows)}"
        )
        if result.get("skipped_gl0"):
            lines.append(
                f"ROIs omitidas (gray_level 0 = aire): {len(result['skipped_gl0'])}"
            )
        lines.append("")

        # 1) Superposicion real de voxeles
        ov = result.get("overlap_voxels", 0)
        if ov == 0:
            lines.append("1) Superposicion de voxeles: 0 (disjuntas verificadas)")
        else:
            lines.append(f"1) Superposicion de voxeles: {ov} -> CONFLICTO")
            for item in result.get("overlap_rois", [])[:15]:
                lines.append(f"   - {item}")

        # 2) gray_level duplicado en la configuracion
        dups = result.get("duplicates", [])
        if not dups:
            lines.append("2) gray_level duplicado: 0 (todos los valores son unicos)")
        else:
            lines.append(f"2) gray_level duplicado: {len(dups)} -> CONFLICTO REAL")
            for gl, names in dups:
                lines.append(
                    f"   - gray_level {gl}: {', '.join(str(n) for n in names)}"
                )

        # 3) Voxeles no mapeados
        unmapped = result.get("unmapped_values", [])
        if not unmapped:
            lines.append(
                "3) Voxeles no mapeados: 0 (todos los valores del volumen pertenecen a una ROI)"
            )
        else:
            lines.append(
                f"3) Valores del volumen sin ROI: {len(unmapped)} valores distintos"
            )
            lines.append(
                f"   {', '.join(str(v) for v in unmapped[:15])}"
                + (" ..." if len(unmapped) > 15 else "")
            )

        # 4) ROIs vacias
        empty = result.get("empty_rois", [])
        if not empty:
            lines.append("4) ROIs vacias: 0")
        else:
            lines.append(f"4) ROIs vacias: {len(empty)} (sin voxeles en el volumen)")
            for name in empty:
                gl = next((r["gray_level"] for r in rows if r["name"] == name), "?")
                lines.append(f"   - {name} (gray_level {gl})")

        # Resumen por ROI
        lines.append("")
        n_ok = sum(1 for r in rows if r.get("estado") == "OK")
        lines.append(
            f"Resumen: {n_ok} ROIs con voxeles, "
            f"{len(empty)} vacias, {len(dups)} conflictos de gray_level, "
            f"{ov} voxeles superpuestos"
        )
        return "\n".join(lines)


# Fin de BNCT.py
