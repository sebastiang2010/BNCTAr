import os, json
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np
import re

#
# BNCTDosimetry
#

class BNCTDosimetry(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "BNCT Dosimetry"
        self.parent.categories = ["BNCT"]
        self.parent.dependencies = []
        self.parent.contributors = ["Antigravity AI / BNCTAr Team"]
        self.parent.helpText = """
Este módulo permite la carga de archivos de dosis PHITS y el cálculo de dosis pesada para BNCT.
"""
        self.parent.acknowledgementText = "Desarrollado para la integración de BNCTAr en 3D Slicer."

#
# BNCTDosimetryWidget
#

class BNCTDosimetryWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        
        # Limpiar la consola al cargar el módulo
        try:
            slicer.app.pythonConsole().clear()
        except:
            pass

        # --- Interfaz de Usuario ---
        parametersCollapsibleButton = ctk.ctkCollapsibleButton()
        parametersCollapsibleButton.text = "Configuración de Dosis BNCT"
        self.layout.addWidget(parametersCollapsibleButton)
        parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

        # Selección de carpeta PHITS Output
        self.phitsDirSelector = ctk.ctkPathLineEdit()
        self.phitsDirSelector.filters = ctk.ctkPathLineEdit.Dirs
        self.phitsDirSelector.nameFilters = []
        parametersFormLayout.addRow("Carpeta PHITS Output:", self.phitsDirSelector)

        # Selección de imagen de referencia
        self.referenceSelector = slicer.qMRMLNodeComboBox()
        self.referenceSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.referenceSelector.selectNodeUponCreation = True
        self.referenceSelector.addEnabled = False
        self.referenceSelector.removeEnabled = False
        self.referenceSelector.noneEnabled = True
        self.referenceSelector.showHidden = False
        self.referenceSelector.setMRMLScene(slicer.mrmlScene)
        self.referenceSelector.connect('currentNodeChanged(vtkMRMLNode*)', self._onReferenceNodeChanged)
        parametersFormLayout.addRow("Imagen de Referencia (CT/MRI):", self.referenceSelector)

        # Concentración de Boro
        self.boronConcentration = qt.QDoubleSpinBox()
        self.boronConcentration.setRange(0, 1000)
        self.boronConcentration.value = 15.0
        parametersFormLayout.addRow("Conc. Boro en Sangre (ppm):", self.boronConcentration)

        self.factorTIR = qt.QDoubleSpinBox()
        self.factorTIR.value = 1.0
        parametersFormLayout.addRow("Factor TIR:", self.factorTIR)

        self.scaleFactor = qt.QDoubleSpinBox()
        self.scaleFactor.setDecimals(15)
        self.scaleFactor.setRange(0, 1e20)
        self.scaleFactor.value = 1.0
        parametersFormLayout.addRow("Factor de Escala (S):", self.scaleFactor)

        # Offsets de coordenadas
        self.offsetX = qt.QDoubleSpinBox()
        self.offsetX.setRange(-5000, 5000)
        self.offsetX.value = 0.0
        parametersFormLayout.addRow("Offset X (mm):", self.offsetX)

        self.offsetY = qt.QDoubleSpinBox()
        self.offsetY.setRange(-5000, 5000)
        self.offsetY.value = 0.0
        parametersFormLayout.addRow("Offset Y (mm):", self.offsetY)

        self.offsetZ = qt.QDoubleSpinBox()
        self.offsetZ.setRange(-5000, 5000)
        self.offsetZ.value = 0.0
        parametersFormLayout.addRow("Offset Z (mm):", self.offsetZ)

        # --- Sección de campos múltiples ---
        self.fieldCountCombo = qt.QComboBox()
        self.fieldCountCombo.addItems(["1 Campo", "2 Campos"])
        parametersFormLayout.addRow("Campos:", self.fieldCountCombo)

        self.peso1Spin = qt.QDoubleSpinBox()
        self.peso1Spin.setRange(0, 100)
        self.peso1Spin.setDecimals(4)
        self.peso1Spin.value = 1.0
        parametersFormLayout.addRow("Peso Campo 1:", self.peso1Spin)

        self.peso2Spin = qt.QDoubleSpinBox()
        self.peso2Spin.setRange(0, 100)
        self.peso2Spin.setDecimals(4)
        self.peso2Spin.value = 1.0
        self.peso2Spin.visible = False
        parametersFormLayout.addRow("Peso Campo 2:", self.peso2Spin)

        self.fieldCountCombo.connect('currentIndexChanged(int)', self._onFieldCountChanged)

        # Botón para interpolar dosis
        self.applyButton = qt.QPushButton("Interpolar Dosis")
        self.applyButton.toolTip = "Alinea centros, carga los archivos y genera el volumen de dosis en Slicer."
        self.applyButton.enabled = True
        parametersFormLayout.addRow(self.applyButton)

        self.physDoseButton = qt.QPushButton("Calcular Dosis Física")
        self.physDoseButton.toolTip = "Selecciona componentes KERMA y calcula dosis física con Tissue/Blood Ratio."
        self.physDoseButton.enabled = True
        parametersFormLayout.addRow(self.physDoseButton)

        self.rbeDoseButton = qt.QPushButton("Calcular Dosis Pesada por RBE")
        self.rbeDoseButton.toolTip = "Selecciona componentes, asigna RBE y calcula dosis pesada."
        self.rbeDoseButton.enabled = True
        parametersFormLayout.addRow(self.rbeDoseButton)

        self.isoDoseButton = qt.QPushButton("Calcular Dosis Isoefectiva")
        self.isoDoseButton.toolTip = "Calcula dosis isoefectiva usando modelo radiobiológico del archivo models.jsonc."
        self.isoDoseButton.enabled = True
        parametersFormLayout.addRow(self.isoDoseButton)

        # Conexiones
        self.applyButton.connect('clicked(bool)', self.onApplyButton)
        self.physDoseButton.connect('clicked(bool)', self.onPhysDose)
        self.rbeDoseButton.connect('clicked(bool)', self.onRBEDose)
        self.isoDoseButton.connect('clicked(bool)', self.onIsoDose)

        # --- Sección Configuración KERMA por ROI ---
        configCollapsibleButton = ctk.ctkCollapsibleButton()
        configCollapsibleButton.text = "Configuración KERMA por ROI"
        configCollapsibleButton.collapsed = False
        self.layout.addWidget(configCollapsibleButton)
        configFormLayout = qt.QFormLayout(configCollapsibleButton)

        btnLayout = qt.QHBoxLayout()
        self.btnSaveMat = qt.QPushButton("Guardar MAT")
        self.btnSaveJson = qt.QPushButton("Guardar JSON")
        self.btnReload = qt.QPushButton("Recargar Default")
        self.btnAddRow = qt.QPushButton("+ Fila")
        self.btnDelRow = qt.QPushButton("- Fila")
        btnLayout.addWidget(self.btnSaveMat)
        btnLayout.addWidget(self.btnSaveJson)
        btnLayout.addWidget(self.btnReload)
        btnLayout.addStretch()
        btnLayout.addWidget(self.btnAddRow)
        btnLayout.addWidget(self.btnDelRow)
        configFormLayout.addRow(btnLayout)

        self.configTable = qt.QTableWidget()
        self.configTable.setColumnCount(10)
        self.configTable.setHorizontalHeaderLabels(
            ["nroi", "nmesh", "Nombre", "RBE_B", "RBE_F", "RBE_T", "RBE_G",
             "Ratio", "Tumor", "GrayLvl"])
        self.configTable.horizontalHeader().setStretchLastSection(True)
        self.configTable.setAlternatingRowColors(True)
        self.configTable.setMinimumHeight(250)
        configFormLayout.addRow(self.configTable)

        self.btnSaveMat.connect('clicked(bool)', self.onSaveMat)
        self.btnSaveJson.connect('clicked(bool)', self.onSaveJson)
        self.btnReload.connect('clicked(bool)', self.onReloadDefault)
        self.btnAddRow.connect('clicked(bool)', self.onAddRow)
        self.btnDelRow.connect('clicked(bool)', self.onRemoveRow)
        self.configTable.cellChanged.connect(self._onConfigCellChanged)

        self._loadDefaultConfig()

        # Cargar modelos radiobiológicos
        self._models_data = None
        try:
            logic = BNCTDosimetryLogic()
            self._models_data = logic.loadModels()
            models_count = len(self._models_data.get('isoeffective', []))
            print(f"Modelos radiobiológicos cargados: {models_count} isoefectivos")
        except Exception as e:
            logging.warning(f"No se pudieron cargar los modelos: {e}")

        # Cargar configuración de órganos
        self._organ_config = self._loadOrganConfig()
        if self._organ_config:
            print(f"Config órganos cargada: {len(self._organ_config)} ROIs")
        else:
            logging.warning("No hay configuración de órganos disponible.")

        # Añadir estirador al final
        self.layout.addStretch(1)

    def _computeCenterOffset(self):
        """Calcula el offset entre el centro de la malla PHITS y la imagen de referencia.
        Retorna (offset_x, offset_y, offset_z) o None si faltan datos."""
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

        bounds = [0]*6
        refNode.GetBounds(bounds)
        center_ras = np.array([
            (bounds[0] + bounds[1]) / 2.0,
            (bounds[2] + bounds[3]) / 2.0,
            (bounds[4] + bounds[5]) / 2.0
        ])

        logic = BNCTDosimetryLogic()
        phits = logic.parsePHITS(n_files[0])
        origin = np.array(phits['origin'])
        shape = np.array(phits['shape'])
        spacing = np.array(phits['spacing'])
        center_phits = origin + (shape[::-1] - 1) * spacing / 2.0

        offset = center_ras - center_phits
        return (offset[0], offset[1], offset[2])

    def _onReferenceNodeChanged(self, node):
        """Auto-centrado al seleccionar imagen de referencia."""
        if not node:
            return
        offset = self._computeCenterOffset()
        if offset is not None:
            self.offsetX.value = offset[0]
            self.offsetY.value = offset[1]
            self.offsetZ.value = offset[2]
            logging.info(f"Auto-centrado: offset={offset}")

    def _onFieldCountChanged(self, idx):
        """Muestra/oculta peso del campo 2 segun seleccion."""
        self.peso2Spin.visible = (idx == 1)  # 1 = "2 Campos"

    def _loadDefaultConfig(self):
        jsonc_path = os.path.join(os.path.dirname(__file__), '..', 'json', 'organ_config.jsonc')
        if os.path.exists(jsonc_path):
            try:
                logic = BNCTDosimetryLogic()
                data = logic.load_jsonc(jsonc_path)
                config = data.get('rois', [])
                if config:
                    self._populateConfigTable(config)
                    print(f"Config KERMA: {len(config)} ROIs desde organ_config.jsonc")
                    return
            except Exception as e:
                logging.warning(f"No se pudo cargar organ_config.jsonc: {e}")
        mat_path = self._defaultMatPath()
        if os.path.exists(mat_path):
            try:
                logic = BNCTDosimetryLogic()
                config = logic.loadConfigFromMat(mat_path)
                self._populateConfigTable(config)
                self._syncOrganConfigJson()
                print(f"Config KERMA: {len(config)} ROIs desde {mat_path}")
            except Exception as e:
                logging.warning(f"Fallback MAT falló: {e}")

    def _defaultMatPath(self):
        return os.path.join(os.path.dirname(__file__), 'mat', 'info_rois.mat')

    def _populateConfigTable(self, config):
        self.configTable.blockSignals(True)
        self.configTable.setRowCount(len(config))
        for i, row in enumerate(config):
            self.configTable.setItem(i, 0, qt.QTableWidgetItem(str(int(row['nroi']))))
            self.configTable.setItem(i, 1, qt.QTableWidgetItem(str(int(row['nmesh']))))
            self.configTable.setItem(i, 2, qt.QTableWidgetItem(str(row['name'])))
            self.configTable.setItem(i, 3, qt.QTableWidgetItem(f"{row['rbe_boron']:.2f}"))
            self.configTable.setItem(i, 4, qt.QTableWidgetItem(f"{row['rbe_fast']:.2f}"))
            self.configTable.setItem(i, 5, qt.QTableWidgetItem(f"{row['rbe_thermal']:.2f}"))
            self.configTable.setItem(i, 6, qt.QTableWidgetItem(f"{row['rbe_gamma']:.2f}"))
            self.configTable.setItem(i, 7, qt.QTableWidgetItem(f"{row['ratio']:.1f}"))
            chk = qt.QTableWidgetItem()
            chk.setCheckState(qt.Qt.Checked if row.get('tumor', False) else qt.Qt.Unchecked)
            self.configTable.setItem(i, 8, chk)
            self.configTable.setItem(i, 9, qt.QTableWidgetItem(str(int(row.get('gray_level', 0)))))
        self.configTable.resizeColumnsToContents()
        self.configTable.blockSignals(False)

    def _readConfigFromTable(self):
        config = []
        for i in range(self.configTable.rowCount):
            def val(col):
                item = self.configTable.item(i, col)
                return item.text() if item else ''
            row = {
                'nroi': int(float(val(0))) if val(0) else 0,
                'nmesh': int(float(val(1))) if val(1) else 0,
                'name': val(2),
                'rbe_boron': float(val(3)) if val(3) else 0.0,
                'rbe_fast': float(val(4)) if val(4) else 0.0,
                'rbe_thermal': float(val(5)) if val(5) else 0.0,
                'rbe_gamma': float(val(6)) if val(6) else 0.0,
                'ratio': float(val(7)) if val(7) else 0.0,
                'tumor': self.configTable.item(i, 8).checkState() == qt.Qt.Checked if self.configTable.item(i, 8) else False,
                'gray_level': int(float(val(9))) if val(9) else 0,
            }
            config.append(row)
        return config

    def _syncOrganConfigJson(self):
        config = self._readConfigFromTable()
        path = os.path.join(os.path.dirname(__file__), '..', 'json', 'organ_config.jsonc')
        try:
            d = os.path.dirname(path)
            if not os.path.exists(d):
                os.makedirs(d)
            with open(path, 'w', encoding='utf-8') as f:
                json.dump({'rois': config}, f, indent=2, ensure_ascii=False)
        except Exception as e:
            logging.warning(f"No se pudo sincronizar organ_config.jsonc: {e}")

    def onSaveMat(self):
        config = self._readConfigFromTable()
        path = self._saveDialogPath("Guardar MAT", "MAT files (*.mat)")
        if not path:
            return
        logic = BNCTDosimetryLogic()
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
        logic = BNCTDosimetryLogic()
        try:
            logic.saveConfigToJson(path, config)
            print(f"Config guardada en {path}")
        except Exception as e:
            slicer.util.errorDisplay(f"Error al guardar JSON: {e}")

    def onReloadDefault(self):
        if os.path.exists(self._defaultMatPath()):
            logic = BNCTDosimetryLogic()
            try:
                config = logic.loadConfigFromMat(self._defaultMatPath())
                self._populateConfigTable(config)
                self._syncOrganConfigJson()
                print(f"Config recargada: {len(config)} ROIs.")
                return
            except ImportError:
                slicer.util.infoDisplay("h5py no encontrado. Instalando...")
                slicer.util.pip_install("h5py")
                config = logic.loadConfigFromMat(self._defaultMatPath())
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

    def _discoverPhitsFiles(self, directory):
        neutron_files = []
        gamma_path = None
        out_files = sorted([f for f in os.listdir(directory)
                            if f.endswith('.out') and '_err' not in f.lower() and '_error' not in f.lower()])
        for f in out_files:
            low = f.lower()
            if 'photon' in low or 'gamma' in low:
                gamma_path = os.path.join(directory, f)
            elif re.match(r'mesh-\d+\.out', low):
                neutron_files.append(os.path.join(directory, f))
        def mesh_number(p):
            m = re.search(r'mesh-(\d+)', os.path.basename(p).lower())
            return int(m.group(1)) if m else 9999
        neutron_files.sort(key=mesh_number)
        return neutron_files, gamma_path

    def onApplyButton(self):
        logic = BNCTDosimetryLogic()
        
        phitsDir = self.phitsDirSelector.currentPath
        refNode = self.referenceSelector.currentNode()
        
        if not phitsDir or not os.path.isdir(phitsDir):
            slicer.util.errorDisplay("Debe seleccionar una carpeta PHITS válida.")
            return
        
        n_files, g_path = self._discoverPhitsFiles(phitsDir)
        if not n_files:
            slicer.util.errorDisplay("No se encontraron archivos mesh-N.out (neutrones) en la carpeta.")
            return
        if not g_path:
            logging.warning("No se encontró archivo gamma. Solo se cargarán neutrones.")

        # Auto-centrado antes de interpolar
        if refNode:
            offset = self._computeCenterOffset()
            if offset is not None:
                self.offsetX.value = offset[0]
                self.offsetY.value = offset[1]
                self.offsetZ.value = offset[2]
                logging.info(f"Auto-centrado en Interpolar Dosis: offset={offset}")

        print(f"--- Archivos PHITS encontrados ---")
        print(f"  Neutrones ({len(n_files)}):")
        for f in n_files:
            print(f"    - {os.path.basename(f)}")
        print(f"  Gamma:     {os.path.basename(g_path) if g_path else 'No encontrado'}")
        print(f"  Carpeta:   {phitsDir}")

        rbe_factors = {
            'offset': (self.offsetX.value, self.offsetY.value, self.offsetZ.value)
        }
        
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic._cleanupOldVolumes()
            # Campo 1
            logic.run(n_files, g_path, refNode, rbe_factors, suffix="")
            # Campo 2 si esta seleccionado
            if self.fieldCountCombo.currentIndex == 1:
                dir_c2 = os.path.join(phitsDir, "campo_2")
                if os.path.isdir(dir_c2):
                    n_files_2, g_path_2 = self._discoverPhitsFiles(dir_c2)
                    if n_files_2:
                        logic.run(n_files_2, g_path_2, refNode, rbe_factors, suffix="_C2")
                    else:
                        print("  campo_2 vacio o sin archivos validos")
                else:
                    print("  No existe campo_2/ en la carpeta PHITS")
        except Exception as e:
            slicer.util.errorDisplay(f"Error durante el cálculo: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def _loadOrganConfig(self):
        """Carga configuración de órganos desde json/organ_config.jsonc
        o fallback a info_rois.mat."""
        # 1. Intentar JSONC (con comentarios)
        jsonc_path = os.path.join(os.path.dirname(__file__), '..', 'json', 'organ_config.jsonc')
        if os.path.exists(jsonc_path):
            try:
                logic = BNCTDosimetryLogic()
                data = logic.load_jsonc(jsonc_path)
                config = data.get('rois', [])
                if config:
                    print(f"Config órganos cargada desde organ_config.jsonc ({len(config)} ROIs)")
                    return config
            except Exception as e:
                logging.warning(f"Fallback a info_rois.mat: {e}")

        # 2. Fallback: leer info_rois.mat directamente
        mat_path = os.path.join(os.path.dirname(__file__), 'mat', 'info_rois.mat')
        if os.path.exists(mat_path):
            try:
                logic = BNCTDosimetryLogic()
                config = logic.loadConfigFromMat(mat_path)
                print(f"Config órganos cargada desde {mat_path} ({len(config)} ROIs)")
                return config
            except Exception as e:
                logging.warning(f"Fallback a info_rois.mat falló: {e}")

        return None

    def _getKermaVolumesForMesh(self, mesh_number):
        """Retorna lista de (node, nombre_corto) con los KERMA volumes que corresponden
        al mesh_number dado. Incluye campo_2 (_C2) si existe."""
        import re
        exclude = ["Dosis_Fisica_BNCT", "Dosis_RBE_BNCT", "Dosis_Final_BNCT", "Dosis_IsoE_BNCT"]
        volumes = []
        for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
            name = node.GetName()
            if name in exclude:
                continue
            if name == 'Gamma_Photon' or name == 'Gamma_Photon_C2':
                volumes.append((node, name))
                continue
            if 'Kerma' in name:
                match = re.search(r'Kerma(\d+)(_C2)?$', name)
                if match and int(match.group(1)) == mesh_number:
                    volumes.append((node, name))
        return volumes

    def _buildOrganDialog(self, title, ask_rbe=False):
        """Dialog con dropdown de órganos desde el cache de configuración.
        Filtra volúmenes por nmesh del órgano seleccionado.
        Retorna (organ_row, rbe_mapped, volumes) o None si cancel."""
        config = self._organ_config
        if not config:
            slicer.util.errorDisplay(
                "No hay configuración de órganos disponible. "
                "Abra el módulo 'BNCT DVH Analysis' y cargue/configure los ROIs primero.")
            return None

        dialog = qt.QDialog(self.parent)
        dialog.setWindowTitle(title)
        dialog.setMinimumWidth(440)
        layout = qt.QVBoxLayout(dialog)
        layout.setSpacing(6)

        # Dropdown de órganos
        organ_layout = qt.QHBoxLayout()
        organ_label = qt.QLabel("Órgano:")
        organ_combo = qt.QComboBox()
        organ_combo.setMinimumWidth(200)
        for row in config:
            organ_combo.addItem(f"{row['name']}  [mesh {int(row['nmesh'])}]")
        organ_layout.addWidget(organ_label)
        organ_layout.addWidget(organ_combo, 1)
        layout.addLayout(organ_layout)

        # Label de mesh info
        mesh_info = qt.QLabel()
        mesh_info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(mesh_info)

        # Valores (read-only)
        values_group = qt.QGroupBox("Parámetros")
        values_layout = qt.QFormLayout(values_group)

        ratio_display = qt.QLabel()
        values_layout.addRow("Ratio T/S:", ratio_display)

        rbe_displays = {}
        if ask_rbe:
            for key, label in [('rbe_boron', 'RBE Boro'), ('rbe_fast', 'RBE Fast'),
                               ('rbe_thermal', 'RBE Thermal'), ('rbe_gamma', 'RBE Gamma')]:
                lbl = qt.QLabel()
                values_layout.addRow(f"{label}:", lbl)
                rbe_displays[key] = lbl

        layout.addWidget(values_group)

        def update_values():
            idx = organ_combo.currentIndex
            if idx < 0 or idx >= len(config):
                return
            row = config[idx]
            ratio_display.text = f"{row['ratio']:.2f}"
            nmesh = int(row['nmesh'])
            mesh_info.text = f"Mesh: {nmesh}"
            if ask_rbe:
                rbe_displays['rbe_boron'].text = f"{row['rbe_boron']:.2f}"
                rbe_displays['rbe_fast'].text = f"{row['rbe_fast']:.2f}"
                rbe_displays['rbe_thermal'].text = f"{row['rbe_thermal']:.2f}"
                rbe_displays['rbe_gamma'].text = f"{row['rbe_gamma']:.2f}"

        organ_combo.connect('currentIndexChanged(int)', lambda idx: update_values())
        update_values()

        # Botones
        btn_layout = qt.QHBoxLayout()
        ok_btn = qt.QPushButton("OK")
        cancel_btn = qt.QPushButton("Cancel")
        btn_layout.addStretch()
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

        ok_btn.connect('clicked()', dialog.accept)
        cancel_btn.connect('clicked()', dialog.reject)

        if dialog.exec() != qt.QDialog.Accepted:
            return None

        idx = organ_combo.currentIndex
        organ_row = config[idx]
        nmesh = int(organ_row['nmesh'])

        # Filtrar volúmenes por mesh
        volumes = self._getKermaVolumesForMesh(nmesh)
        if not volumes:
            slicer.util.errorDisplay(
                f"No se encontraron componentes KERMA para mesh {nmesh}. "
                f"Ejecute 'Interpolar Dosis' primero.")
            return None

        # Mapear RBE a cada volumen según su tipo
        rbe_mapped = []
        for node, vol_name in volumes:
            if 'Boron' in vol_name or 'boro' in vol_name.lower():
                rbe_mapped.append(organ_row['rbe_boron'])
            elif 'Fast' in vol_name:
                rbe_mapped.append(organ_row['rbe_fast'])
            elif 'Thermal' in vol_name or 'thn' in vol_name.lower():
                rbe_mapped.append(organ_row['rbe_thermal'])
            elif 'Gamma' in vol_name or 'photon' in vol_name.lower():
                rbe_mapped.append(organ_row['rbe_gamma'])
            else:
                rbe_mapped.append(1.0)

        return (organ_row, rbe_mapped, volumes)

    def onPhysDose(self):
        """Calcula Dosis Física usando órgano de la tabla KERMA."""
        result = self._buildOrganDialog("Dosis Física", ask_rbe=False)
        if result is None:
            return
        organ_row, _, volumes = result
        nodes = [node for node, _ in volumes]

        logic = BNCTDosimetryLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculatePhysicalDose(
                nodes=nodes,
                b10_conc=self.boronConcentration.value,
                ratio_tb=organ_row['ratio'],
                tir=self.factorTIR.value,
                scale=self.scaleFactor.value,
                peso1=self.peso1Spin.value,
                peso2=self.peso2Spin.value
            )
            slicer.util.infoDisplay("Dosis Física calculada correctamente.")
        except Exception as e:
            slicer.util.errorDisplay(f"Error: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def onRBEDose(self):
        """Calcula Dosis RBE usando órgano de la tabla KERMA."""
        result = self._buildOrganDialog("Dosis Pesada por RBE", ask_rbe=True)
        if result is None:
            return
        organ_row, rbe_mapped, volumes = result
        nodes = [node for node, _ in volumes]

        logic = BNCTDosimetryLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculateRBEDose(
                nodes=nodes,
                rbe_values=rbe_mapped,
                b10_conc=self.boronConcentration.value,
                ratio_tb=organ_row['ratio'],
                tir=self.factorTIR.value,
                scale=self.scaleFactor.value,
                peso1=self.peso1Spin.value,
                peso2=self.peso2Spin.value
            )
            slicer.util.infoDisplay("Dosis RBE calculada correctamente.")
        except Exception as e:
            slicer.util.errorDisplay(f"Error: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

    def _buildIsoDialog(self):
        """Dialog con dropdown de órgano (del cache de configuración, filtra por nmesh)
        + dropdown de modelo isoefectivo (del JSONC) + Tumor/blood B10 ratio.
        Retorna (params_arr, bratio, volumes, organ_row) o None si cancel."""
        config = self._organ_config
        if not config:
            slicer.util.errorDisplay(
                "No hay configuración de órganos disponible. "
                "Abra el módulo 'BNCT DVH Analysis' y cargue/configure los ROIs primero.")
            return None

        if not self._models_data:
            slicer.util.errorDisplay(
                "No se pudieron cargar los modelos radiobiológicos. "
                "Verifique que existe json/models.jsonc.")
            return None

        models = self._models_data.get('isoeffective', [])
        if not models:
            slicer.util.errorDisplay("No hay modelos isoefectivos en models.jsonc.")
            return None

        dialog = qt.QDialog(self.parent)
        dialog.setWindowTitle("Dosis Isoefectiva")
        dialog.setMinimumWidth(520)
        layout = qt.QVBoxLayout(dialog)
        layout.setSpacing(6)

        # Dropdown de órganos (desde tabla KERMA)
        organ_layout = qt.QHBoxLayout()
        organ_label = qt.QLabel("Órgano:")
        organ_combo = qt.QComboBox()
        organ_combo.setMinimumWidth(200)
        for row in config:
            organ_combo.addItem(f"{row['name']}  [mesh {int(row['nmesh'])}]")
        organ_layout.addWidget(organ_label)
        organ_layout.addWidget(organ_combo, 1)
        layout.addLayout(organ_layout)

        # Mesh info
        mesh_info = qt.QLabel()
        mesh_info.setStyleSheet("color: #888; font-size: 10px;")
        layout.addWidget(mesh_info)

        # Dropdown de modelos (del JSONC)
        model_layout = qt.QHBoxLayout()
        model_label = qt.QLabel("Modelo:")
        model_combo = qt.QComboBox()
        model_combo.setMinimumWidth(300)
        for m in models:
            model_combo.addItem(f"{m['id']}. {m['name']}")
        model_layout.addWidget(model_label)
        model_layout.addWidget(model_combo, 1)
        layout.addLayout(model_layout)

        # Info del modelo seleccionado
        info_label = qt.QLabel()
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(info_label)

        # Parámetros del modelo (read-only)
        params_group = qt.QGroupBox("Parámetros del modelo")
        params_layout = qt.QFormLayout(params_group)
        param_displays = {}
        param_keys = [
            ('alpha_r', 'α_r'), ('beta_r', 'β_r'),
            ('alpha_boro', 'α_boro'), ('beta_boro', 'β_boro'),
            ('alpha_thn', 'α_thn'), ('beta_thn', 'β_thn'),
            ('alpha_fast', 'α_fast'), ('beta_fast', 'β_fast'),
            ('tof', 'tof'), ('tos', 'tos'),
            ('pf_g', 'pf_g'), ('ps_g', 'ps_g'),
            ('pf_bnct', 'pf_bnct'), ('ps_bnct', 'ps_bnct')
        ]
        for key, display in param_keys:
            lbl = qt.QLabel()
            params_layout.addRow(f"{display}:", lbl)
            param_displays[key] = lbl
        layout.addWidget(params_group)

        # Tumor/blood B10 ratio
        bratio_layout = qt.QHBoxLayout()
        bratio_label = qt.QLabel("Tumor/blood B10 ratio:")
        bratio_spin = qt.QDoubleSpinBox()
        bratio_spin.setRange(0, 100)
        bratio_spin.setDecimals(2)
        bratio_spin.value = 3.5
        bratio_spin.setFixedWidth(120)
        bratio_layout.addWidget(bratio_label)
        bratio_layout.addWidget(bratio_spin)
        bratio_layout.addStretch()
        layout.addLayout(bratio_layout)

        def update_model():
            idx = model_combo.currentIndex
            if idx < 0 or idx >= len(models):
                return
            m = models[idx]
            info_label.text = m.get('ref', '')
            for key in param_keys:
                k = key[0]
                val = m.get(k, 0)
                if abs(val) < 1e-10:
                    param_displays[k].text = "0"
                else:
                    param_displays[k].text = f"{val:.6g}"

        def update_organ():
            idx = organ_combo.currentIndex
            if idx < 0 or idx >= len(config):
                return
            row = config[idx]
            mesh_info.text = f"Mesh: {int(row['nmesh'])}"

        organ_combo.connect('currentIndexChanged(int)', lambda i: update_organ())
        model_combo.connect('currentIndexChanged(int)', lambda i: update_model())
        update_organ()
        update_model()

        # Botones
        btn_layout = qt.QHBoxLayout()
        ok_btn = qt.QPushButton("OK")
        cancel_btn = qt.QPushButton("Cancel")
        btn_layout.addStretch()
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

        ok_btn.connect('clicked()', dialog.accept)
        cancel_btn.connect('clicked()', dialog.reject)

        if dialog.exec() != qt.QDialog.Accepted:
            return None

        org_idx = organ_combo.currentIndex
        mdl_idx = model_combo.currentIndex
        organ_row = config[org_idx]
        nmesh = int(organ_row['nmesh'])
        model_params = models[mdl_idx]
        bratio = bratio_spin.value

        # Filtrar volúmenes por mesh
        volumes = self._getKermaVolumesForMesh(nmesh)
        if not volumes:
            slicer.util.errorDisplay(
                f"No se encontraron componentes KERMA para mesh {nmesh}. "
                f"Ejecute 'Interpolar Dosis' primero.")
            return None

        # Construir array de 14 parámetros como en MATLAB
        params_arr = [
            model_params['alpha_r'], model_params['beta_r'],
            model_params['alpha_boro'], model_params['beta_boro'],
            model_params['alpha_thn'], model_params['beta_thn'],
            model_params['alpha_fast'], model_params['beta_fast'],
            model_params['tof'], model_params['tos'],
            model_params['pf_g'], model_params['ps_g'],
            model_params['pf_bnct'], model_params['ps_bnct']
        ]
        return (params_arr, bratio, volumes, organ_row)

    def onIsoDose(self):
        """Calcula Dosis Isoefectiva usando órgano + modelo radiobiológico."""
        result = self._buildIsoDialog()
        if result is None:
            return
        params_arr, bratio, volumes, organ_row = result
        nodes = [node for node, _ in volumes]
        # Mapear componentes: boro, thn, fast, gamma
        comp_arrays = {'boro': None, 'thn': None, 'fast': None, 'gamma': None}
        for node, vol_name in volumes:
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if 'Boron' in vol_name or 'boro' in vol_name.lower():
                if comp_arrays['boro'] is None:
                    comp_arrays['boro'] = arr
                else:
                    comp_arrays['boro'] += arr
            elif 'Thermal' in vol_name or 'thn' in vol_name.lower():
                if comp_arrays['thn'] is None:
                    comp_arrays['thn'] = arr
                else:
                    comp_arrays['thn'] += arr
            elif 'Fast' in vol_name:
                if comp_arrays['fast'] is None:
                    comp_arrays['fast'] = arr
                else:
                    comp_arrays['fast'] += arr
            elif 'Gamma' in vol_name or 'photon' in vol_name.lower():
                if comp_arrays['gamma'] is None:
                    comp_arrays['gamma'] = arr
                else:
                    comp_arrays['gamma'] += arr

        logic = BNCTDosimetryLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculateIsoEffectveDose(
                comp_arrays=comp_arrays,
                params=params_arr,
                bratio=bratio,
                tir=self.factorTIR.value,
                scale=self.scaleFactor.value,
                origin=nodes[0].GetOrigin(),
                spacing=nodes[0].GetSpacing()
            )
            slicer.util.infoDisplay("Dosis Isoefectiva calculada correctamente.")
        except Exception as e:
            slicer.util.errorDisplay(f"Error: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
        finally:
            qt.QApplication.restoreOverrideCursor()

#
# BNCTDosimetryLogic
#

class BNCTDosimetryLogic(ScriptedLoadableModuleLogic):

    @staticmethod
    def load_jsonc(path):
        """Carga un archivo JSONC (JSON con comentarios // y /* */)."""
        import re
        with open(path) as f:
            text = f.read()
        text = re.sub(r'//.*?\n|/\*.*?\*/', '', text, flags=re.S)
        return json.loads(text)

    @staticmethod
    def loadModels():
        """Carga los modelos radiobiológicos desde json/models.jsonc.
        Retorna dict con claves 'isoeffective', 'tcp', 'ntcp'."""
        models_path = os.path.join(os.path.dirname(__file__), '..', 'json', 'models.jsonc')
        if not os.path.exists(models_path):
            raise FileNotFoundError(f"No se encuentra {models_path}")
        return BNCTDosimetryLogic.load_jsonc(models_path)

    def loadConfigFromMat(self, path):
        """Carga configuración KERMA desde archivo MAT (info_rois.mat).
        Para fallback cuando json/organ_config.json no existe."""
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
                return data.decode('latin-1')
            if isinstance(data, np.ndarray) and data.dtype.kind == 'u':
                return ''.join(chr(int(c)) for c in data.ravel())
            if isinstance(data, np.ndarray):
                arr = data.ravel()
                return arr[0] if len(arr) == 1 else arr
            if isinstance(data, h5py.Group):
                return data
            return float(data) if np.isscalar(data) else data

        def read_rbe(rbe_group):
            if isinstance(rbe_group, h5py.Group):
                return {
                    'rbe_boron': float(read_ref(rbe_group['boron'][0, 0])),
                    'rbe_fast': float(read_ref(rbe_group['fast_neutron'][0, 0])),
                    'rbe_thermal': float(read_ref(rbe_group['nitrogen'][0, 0])),
                    'rbe_gamma': float(read_ref(rbe_group['photon'][0, 0])),
                }
            return {'rbe_boron': 1.3, 'rbe_fast': 3.2, 'rbe_thermal': 3.2, 'rbe_gamma': 1.0}

        with h5py.File(path, 'r') as f:
            info = f['info_dvh_cnea']
            num_rows = len(info['nroi'])
            for i in range(num_rows):
                rbe_grp = info['RBE'][i, 0]
                rbe_vals = read_rbe(rbe_grp) if isinstance(rbe_grp, h5py.Group) else read_rbe(rbe_grp)
                row = {
                    'nroi': int(float(read_ref(info['nroi'][i, 0]))),
                    'nmesh': int(float(read_ref(info['nmesh'][i, 0]))),
                    'name': str(read_ref(info['roi_name'][i, 0])),
                    'rbe_boron': float(rbe_vals['rbe_boron']),
                    'rbe_fast': float(rbe_vals['rbe_fast']),
                    'rbe_thermal': float(rbe_vals['rbe_thermal']),
                    'rbe_gamma': float(rbe_vals['rbe_gamma']),
                    'ratio': float(read_ref(info['ratio'][i, 0])),
                    'tir': float(read_ref(info['tir'][i, 0])),
                    'tumor': bool(int(float(read_ref(info['tumor'][i, 0])))),
                    'gray_level': int(float(read_ref(info['gray_level'][i, 0]))),
                }
                config.append(row)
        return config

    def saveConfigToMat(self, path, config):
        try:
            import h5py
        except ImportError:
            slicer.util.pip_install("h5py")
            import h5py
        import numpy as np
        from datetime import datetime
        n = len(config)
        with h5py.File(path, 'w') as f:
            grp = f.create_group('info_dvh_cnea')
            for key in ['nroi', 'nmesh', 'gray_level', 'tumor']:
                grp.create_dataset(key, (n, 1), dtype='i4')
            for key in ['ratio', 'tir', 'Gy2MU', 'vol', 'Unit']:
                grp.create_dataset(key, (n, 1), dtype='f8')
            for key in ['roi_name', 'fecha', 'X', 'Y']:
                grp.create_dataset(key, (n, 1), dtype=h5py.string_dtype())
            for key in ['Max', 'Mean', 'Min']:
                grp.create_dataset(key, (n, 1), dtype='f8')
            for i, row in enumerate(config):
                grp['nroi'][i, 0] = row.get('nroi', i+1)
                grp['nmesh'][i, 0] = row.get('nmesh', 1)
                grp['gray_level'][i, 0] = row.get('gray_level', 0)
                grp['tumor'][i, 0] = 1 if row.get('tumor', False) else 0
                grp['ratio'][i, 0] = row.get('ratio', 1.0)
                grp['tir'][i, 0] = row.get('tir', 20.0)
                grp['Gy2MU'][i, 0] = 0.0
                grp['vol'][i, 0] = 0.0
                grp['Unit'][i, 0] = 1.0
                grp['roi_name'][i, 0] = str(row.get('name', f'ROI_{i+1}'))
                grp['fecha'][i, 0] = datetime.now().strftime('%d-%b-%Y')
                grp['X'][i, 0] = ''
                grp['Y'][i, 0] = ''
                grp['Max'][i, 0] = 0.0
                grp['Mean'][i, 0] = 0.0
                grp['Min'][i, 0] = 0.0
            rbe_dt = h5py.special_dtype(ref=h5py.Reference)
            rbe_dataset = grp.create_dataset('RBE', (n, 1), dtype=rbe_dt)
            for i in range(n):
                rbe_name = chr(65 + i) if i < 26 else f'rbe_{i}'
                rbe_grp = grp.create_group(f'#{rbe_name}')
                rbe_grp.create_dataset('boron', data=float(config[i].get('rbe_boron', 1.3)))
                rbe_grp.create_dataset('fast_neutron', data=float(config[i].get('rbe_fast', 3.2)))
                rbe_grp.create_dataset('nitrogen', data=float(config[i].get('rbe_thermal', 3.2)))
                rbe_grp.create_dataset('photon', data=float(config[i].get('rbe_gamma', 1.0)))
                rbe_dataset[i, 0] = rbe_grp.ref

    def saveConfigToJson(self, path, config):
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)

    def run(self, neutron_files, gamma_path, refNode, rbe_factors, suffix=""):
        logging.info(f'Iniciando proceso de dosis multicampo (suffix="{suffix}")...')

        if not neutron_files:
            if not suffix:
                slicer.util.errorDisplay("No se encontraron archivos mesh-N.out (neutrones) en la carpeta.")
            return

        print(f"--- Archivos PHITS encontrados (suffix='{suffix}') ---")
        print(f"  KERMA ({len(neutron_files)}):")
        for f in neutron_files:
            print(f"    - {os.path.basename(f)}")
        print(f"  Gamma: {os.path.basename(gamma_path) if gamma_path else 'No encontrado'}")

        final_origin = None
        final_spacing = None

        for i, npath in enumerate(neutron_files, start=1):
            print(f"--- Procesando {os.path.basename(npath)} (Field {i}) ---")

            n_data = self.parsePHITS(npath)

            if n_data['particle'] != 'neutron':
                logging.warning(f"  Saltando {os.path.basename(npath)}: partícula={n_data['particle']}")
                continue

            if final_origin is None and not refNode:
                final_origin = n_data['origin']
                final_spacing = n_data['spacing']

            if refNode:
                n_data['user_offset'] = rbe_factors['offset']
                for j in range(len(n_data['data'])):
                    n_data['data'][j] = self.interpolateDose(n_data['data'][j], n_data, refNode)
                final_origin = refNode.GetOrigin()
                final_spacing = refNode.GetSpacing()

            print(f"Field {i}: final_origin={final_origin}, final_spacing={final_spacing}")
            
            comp_names = ["Boron", "FastNeutron", "ThermalNeutron"]
            for j, comp in enumerate(comp_names):
                vol_name = f"{comp}_Kerma{i}{suffix}"
                self.createDoseVolume(n_data['data'][j], final_origin, final_spacing, vol_name)

        has_gamma = False
        gamma_name = "Gamma_Photon"
        if suffix:
            gamma_name = f"Gamma_Photon{suffix}"
        if gamma_path and os.path.exists(gamma_path):
            print(f"--- Procesando gamma: {os.path.basename(gamma_path)} ---")
            g_data = self.parsePHITS(gamma_path)

            if refNode:
                g_data['user_offset'] = rbe_factors['offset']
                g_data['data'][0] = self.interpolateDose(g_data['data'][0], g_data, refNode)
                final_origin = refNode.GetOrigin()
                final_spacing = refNode.GetSpacing()
                print(f"Gamma: final_origin={final_origin}, final_spacing={final_spacing}")

            self.createDoseVolume(g_data['data'][0], final_origin, final_spacing, gamma_name)
            has_gamma = True
        else:
            if not suffix:
                logging.warning("No se encontró archivo gamma. Solo se cargarán campos de neutrones.")

        self._organizeInSubjectHierarchy(len(neutron_files), has_gamma)

        print(f"--- Carga completa (suffix='{suffix}') ---")
        print(f"  KERMA components: {len(neutron_files)}")
        print(f"  Volúmenes creados: {len(neutron_files) * 3 + (1 if has_gamma else 0)}")
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
        return 'Boron' in name or 'boro' in name.lower()

    def calculatePhysicalDose(self, nodes, b10_conc=15.0, ratio_tb=3.5, tir=1.0, scale=1.0,
                               peso1=1.0, peso2=1.0):
        """Dosis Física = suma de componentes KERMA seleccionados con B10*ratio.
        peso1, peso2: pesos para campo_1 y campo_2 (detectado por sufijo _C2)."""
        if not nodes:
            raise ValueError("No hay volúmenes KERMA seleccionados.")

        total = np.zeros_like(slicer.util.arrayFromVolume(nodes[0]).astype(np.float64))
        for node in nodes:
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if self._isBoronComponent(node.GetName()):
                arr = arr * b10_conc * ratio_tb
            # Aplicar peso del campo
            if '_C2' in node.GetName():
                arr *= peso2
            else:
                arr *= peso1
            total += arr

        phys_dose = total * tir * scale

        origin = nodes[0].GetOrigin()
        spacing = nodes[0].GetSpacing()
        self.createDoseVolume(phys_dose, origin, spacing, "Dosis_Fisica_BNCT")

        print(f"Dosis Física: max={np.max(phys_dose):.4f}, min={np.min(phys_dose):.4f}")

    def calculateRBEDose(self, nodes, rbe_values, b10_conc=15.0, ratio_tb=3.5, tir=1.0, scale=1.0,
                          peso1=1.0, peso2=1.0):
        """Dosis RBE = suma de componentes KERMA ponderados por RBE.
        peso1, peso2: pesos para campo_1 y campo_2."""
        if not nodes:
            raise ValueError("No hay volúmenes KERMA seleccionados.")
        if len(nodes) != len(rbe_values):
            raise ValueError("Cantidad de nodos y valores RBE no coincide.")

        total = np.zeros_like(slicer.util.arrayFromVolume(nodes[0]).astype(np.float64))
        for node, rbe in zip(nodes, rbe_values):
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if self._isBoronComponent(node.GetName()):
                arr = arr * b10_conc * ratio_tb
            # Aplicar peso del campo
            if '_C2' in node.GetName():
                arr *= peso2
            else:
                arr *= peso1
            total += arr * rbe

        rbe_dose = total * tir * scale

        origin = nodes[0].GetOrigin()
        spacing = nodes[0].GetSpacing()
        self.createDoseVolume(rbe_dose, origin, spacing, "Dosis_RBE_BNCT")

        print(f"Dosis RBE: max={np.max(rbe_dose):.4f}, min={np.min(rbe_dose):.4f}")

    def _leaCatchesideG(self, tof, tos, pf, ps, time_val):
        """Calcula el factor G de Lea-Catcheside para un tiempo dado."""
        xf = tof / time_val
        xs = tos / time_val
        Gf = 2 * xf * (1 - xf * (1 - np.exp(-1 / xf)))
        Gs = 2 * xs * (1 - xs * (1 - np.exp(-1 / xs)))
        return Gs - (pf * Gf + ps * Gs - ps * Gf)

    def calculateIsoEffectveDose(self, comp_arrays, params, bratio, tir=1.0,
                                   scale=1.0, origin=(0,0,0), spacing=(1,1,1)):
        """
        Dosis Isoefectiva = modelo LQ con Lea-Catcheside.
        Traducido de isodose.m (Gonzalez - Santa Cruz 2012, Radiat. Res. 178)

        comp_arrays: dict con 'boro', 'thn', 'fast', 'gamma' (c/u es array 3D o None)
        params: lista de 14 parámetros [alpha_r, beta_r, alpha_boro, beta_boro, ...]
        bratio: tumor/blood boron concentration ratio
        """
        if comp_arrays['boro'] is None:
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

        # Gr: photon reference time factor (time=30 min)
        Gr = self._leaCatchesideG(tof, tos, pf_g, 1 - pf_g, 30.0)

        # Gij: BNCT time factors
        xf = tof / tir
        xs = tos / tir
        Gf = 2 * xf * (1 - xf * (1 - np.exp(-1 / np.maximum(xf, 1e-10))))
        Gs = 2 * xs * (1 - xs * (1 - np.exp(-1 / np.maximum(xs, 1e-10))))

        # Componentes físicas
        D_boro = comp_arrays['boro'] * bratio
        D_thn = comp_arrays['thn'] if comp_arrays['thn'] is not None else np.zeros_like(D_boro)
        D_fast = comp_arrays['fast'] if comp_arrays['fast'] is not None else np.zeros_like(D_boro)
        D_g = comp_arrays['gamma'] if comp_arrays['gamma'] is not None else np.zeros_like(D_boro)

        Dt = D_boro + D_thn + D_fast + D_g

        with np.errstate(divide='ignore', invalid='ignore'):
            fb = np.divide(D_boro, Dt, out=np.zeros_like(Dt), where=Dt > 0)
            fn = np.divide(D_thn + D_fast, Dt, out=np.zeros_like(Dt), where=Dt > 0)
            fg = np.divide(D_g, Dt, out=np.zeros_like(Dt), where=Dt > 0)

        # Proporciones relativas
        denom_12 = fb + fn
        a_12 = np.divide(fb, denom_12, out=np.zeros_like(fb), where=denom_12 > 0)
        a_21 = np.divide(fn, denom_12, out=np.zeros_like(fn), where=denom_12 > 0)

        denom_13 = fb + fg
        a_13 = np.divide(fb, denom_13, out=np.zeros_like(fb), where=denom_13 > 0)
        a_31 = np.divide(fg, denom_13, out=np.zeros_like(fg), where=denom_13 > 0)

        denom_23 = fn + fg
        a_23 = np.divide(fn, denom_23, out=np.zeros_like(fn), where=denom_23 > 0)
        a_32 = np.divide(fg, denom_23, out=np.zeros_like(fg), where=denom_23 > 0)

        # Gii
        G_11 = pf_bnct * Gf + ps_bnct * Gs
        G_22 = G_11
        G_33 = pf_g * Gf + ps_g * Gs

        # Gij
        G_12 = Gs - (a_12 * pf_bnct + a_21 * pf_bnct) * (Gs - Gf)
        G_23 = Gs - (a_23 * pf_bnct + a_32 * pf_g) * (Gs - Gf)
        G_31 = Gs - (a_31 * pf_g + a_13 * pf_bnct) * (Gs - Gf)

        # Término lineal
        D_lineal = (alpha[0] * D_boro + alpha[1] * D_thn +
                    alpha[2] * D_fast + alpha[3] * D_g) * tir

        # Término cuadrático
        D_cuad = (G_11 * beta_mat[0] * D_boro**2 +
                  G_22 * beta_mat[1] * D_thn**2 +
                  G_22 * beta_mat[2] * D_fast**2 +
                  G_33 * beta_mat[3] * D_g**2 +
                  2 * G_12 * np.sqrt(beta_mat[0] * beta_mat[1]) * D_boro * D_thn +
                  2 * G_12 * np.sqrt(beta_mat[0] * beta_mat[2]) * D_boro * D_fast +
                  2 * G_31 * np.sqrt(beta_mat[0] * beta_mat[3]) * D_boro * D_g +
                  2 * G_11 * np.sqrt(beta_mat[1] * beta_mat[2]) * D_thn * D_fast +
                  2 * G_23 * np.sqrt(beta_mat[1] * beta_mat[3]) * D_thn * D_g +
                  2 * G_23 * np.sqrt(beta_mat[2] * beta_mat[3]) * D_fast * D_g) * tir**2

        # Dosis isoefectiva
        DisoE = 0.5 * (alphabeta_r / Gr) * (
            np.sqrt(1 + 4 * Gr * beta_r / alpha_r**2 * (D_lineal + D_cuad)) - 1
        ) * scale

        # Reemplazar NaN/Inf por 0
        DisoE = np.nan_to_num(DisoE)

        self.createDoseVolume(DisoE, origin, spacing, "Dosis_IsoE_BNCT")
        print(f"Dosis Isoefectiva: max={np.max(DisoE):.4f}, min={np.min(DisoE):.4f}")

    def interpolateDose(self, dose_array, source_data, target_node):
        """Remuestrea la dosis usando la matriz de transformación real de Slicer"""
        try:
            from scipy.interpolate import RegularGridInterpolator
        except ImportError:
            slicer.util.pip_install("scipy")
            from scipy.interpolate import RegularGridInterpolator

        # 1. Malla original de PHITS (en el espacio de la simulación)
        nz, ny, nx = dose_array.shape
        ox, oy, oz = source_data['origin']
        sx, sy, sz = source_data['spacing']
        
        # Coordenadas de los centros de los voxeles de PHITS
        x_src = ox + np.arange(nx) * sx
        y_src = oy + np.arange(ny) * sy
        z_src = oz + np.arange(nz) * sz

        # 2. Obtener la rejilla de la imagen de destino (CT/MRI)
        t_dims = target_node.GetImageData().GetDimensions() # (X, Y, Z)
        t_origin = target_node.GetOrigin()
        t_spacing = target_node.GetSpacing()
        ox_off, oy_off, oz_off = source_data.get('user_offset', (0,0,0))
        
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

        interp = RegularGridInterpolator((z_src, y_src, x_src), dose_array, 
                                        method='linear', bounds_error=False, fill_value=np.nan)
        
        # Procesar slice por slice en Z para no saturar memoria
        i_1d = np.arange(t_dims[0])
        j_1d = np.arange(t_dims[1])
        jj, ii = np.meshgrid(j_1d, i_1d, indexing='ij')
        
        # Pre-calcular matriz IJKToRAS una sola vez
        mat = vtk.vtkMatrix4x4()
        target_node.GetIJKToRASMatrix(mat)
        m_np = slicer.util.arrayFromVTKMatrix(mat)
        
        d_interp = np.zeros((t_dims[2], t_dims[1], t_dims[0]))
        
        for k_idx in range(t_dims[2]):
            pts_ijk = np.vstack([ii.ravel(), jj.ravel(),
                                  np.full(ii.size, k_idx), np.ones(ii.size)])
            pts_ras = m_np @ pts_ijk
            
            x_pts = pts_ras[0, :] - ox_off
            y_pts = pts_ras[1, :] - oy_off
            z_pts = pts_ras[2, :] - oz_off
            
            pts_interp = np.vstack([z_pts, y_pts, x_pts]).T
            d_interp[k_idx, :, :] = interp(pts_interp).reshape((t_dims[1], t_dims[0]))
            
            if (k_idx + 1) % 20 == 0 or k_idx == 0:
                print(f"  Slice {k_idx+1}/{t_dims[2]} interpolated")
        
        d_interp = np.nan_to_num(d_interp)
        
        print(f"Dosis Interp Max: {np.max(d_interp):.6e}, Min: {np.min(d_interp):.6e}")
        return d_interp

    def parsePHITS(self, path):
        """Parser real de PHITS (Traducido de f_lectura_phits_neutron.m)"""
        logging.info(f"Analizando archivo: {path}")
        
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        # Extraer dimensiones (nx, ny, nz) y mallas (xmin, xdel...)
        params = {}
        patterns = {
            'xmin': r'xmin\s*=\s*([\d\.\-\+eE]+)', 'xmax': r'xmax\s*=\s*([\d\.\-\+eE]+)', 'xdel': r'xdel\s*=\s*([\d\.\-\+eE]+)',
            'ymin': r'ymin\s*=\s*([\d\.\-\+eE]+)', 'ymax': r'ymax\s*=\s*([\d\.\-\+eE]+)', 'ydel': r'ydel\s*=\s*([\d\.\-\+eE]+)',
            'zmin': r'zmin\s*=\s*([\d\.\-\+eE]+)', 'zmax': r'zmax\s*=\s*([\d\.\-\+eE]+)', 'zdel': r'zdel\s*=\s*([\d\.\-\+eE]+)',
            'nx': r'nx\s*=\s*(\d+)', 'ny': r'ny\s*=\s*(\d+)', 'nz': r'nz\s*=\s*(\d+)'
        }

        for line in lines:
            for key, pattern in patterns.items():
                if key not in params:
                    match = re.search(pattern, line)
                    if match:
                        params[key] = float(match.group(1))

        nx, ny, nz = int(params['nx']), int(params['ny']), int(params['nz'])
        # --- CONVERSIÓN CM -> MM ---
        spacing = (params['xdel'] * 10.0, params['ydel'] * 10.0, params['zdel'] * 10.0)
        origin = (
            (params['xmin'] + params['xdel']/2.0) * 10.0, 
            (params['ymin'] + params['ydel']/2.0) * 10.0, 
            (params['zmin'] + params['zdel']/2.0) * 10.0
        )

        # Detectar tipo de partícula
        particle = 'unknown'
        for line in lines:
            m = re.search(r'^\s*part\s*=\s*(\w+)', line)
            if m:
                particle = m.group(1).lower()
                break

        # En PHITS, los datos pueden empezar tras 'newpage:' (Neutrones) o 'hc' (Gamma)
        # IMPORTANTE: No buscar ambos a la vez para evitar doble conteo si ambos están presentes
        newpage_indices = [i for i, line in enumerate(lines) if 'newpage:' in line.lower()]
        if not newpage_indices:
            newpage_indices = [i for i, line in enumerate(lines) if 'hc' in line.lower()]
        
        if not newpage_indices:
            raise ValueError("No se encontraron bloques de datos (newpage/hc) en el archivo.")

        # Replicamos exactamente el salto y cálculo de MATLAB
        # salto=25; m=round(nx*ny/10); m1=m-1;
        total_pages = len(newpage_indices)
        num_msets = total_pages // nz
        if num_msets > 3: num_msets = 3
        elif num_msets == 0 and total_pages > 0: num_msets = 1

        import math
        m_matlab = int(math.ceil(nx * ny / 10.0))
        m1 = m_matlab - 1
        
        msets_data = []
        
        def robust_float(s):
            try:
                s_clean = s.replace('f', '.').replace('d', 'e').replace('D', 'e')
                s_clean = re.sub(r'(\d)([\+\-])(\d)', r'\1E\2\3', s_clean)
                return float(s_clean)
            except:
                return 0.0 # MATLAB usa 0 o NaN aquí

        for m in range(num_msets):
            mset_array = np.zeros((nz, ny, nx)) # En Slicer queremos (nz, ny, nx)
            logging.info(f"Cargando mset {m+1} al estilo MATLAB...")
            
            for k in range(nz):
                page_idx = m * nz + k
                if page_idx >= len(newpage_indices): break
                
                idx = newpage_indices[page_idx]
                
                # Determinamos el salto dinámicamente buscando la línea 'hc:'
                salto = 1
                if 'newpage:' in lines[idx].lower():
                    # Buscamos el encabezado de datos 'hc:' en las siguientes líneas
                    for offset in range(1, 50):
                        if (idx + offset) < len(lines) and 'hc:' in lines[idx + offset].lower():
                            salto = offset + 1
                            break
                    else:
                        # Si no encontramos hc:, usamos el salto estándar de 25 como fallback
                        salto = 25
                
                # MATLAB: A=lines(a(k)+salto:a(k)+salto+m1);
                # Esto lee 'm_matlab' líneas
                A_lines = lines[idx + salto : idx + salto + m_matlab]
                
                # MATLAB: B10 = first lines, B11 = last line padded to 10
                all_values = []
                for line in A_lines:
                    vals = [robust_float(v) for v in line.split()]
                    all_values.extend(vals)
                    # Si la línea tiene menos de 10, pad con 0/NaN como B12 de MATLAB
                    if len(vals) < 10:
                        all_values.extend([0.0] * (10 - len(vals)))
                
                # MATLAB: A2 = A1(:) -> Stream de valores
                # B(:,j) = A2(1:nx)'
                B_matlab = np.zeros((nx, ny))
                for j in range(ny):
                    start_ptr = j * nx
                    B_matlab[:, j] = all_values[start_ptr : start_ptr + nx]
                
                # B1 = B (en MATLAB B1 es nx x ny)
                # En Slicer, para que coincida visualmente y con las coordenadas (Ymin en index 0),
                # debemos hacer flip vertical ya que PHITS entrega de Ymax a Ymin.
                mset_array[k, :, :] = np.flipud(B_matlab.T)
            
            msets_data.append(mset_array)

        return {
            'data': msets_data,
            'origin': origin,
            'spacing': spacing,
            'shape': (nz, ny, nx),
            'particle': particle
        }

    def createDoseVolume(self, array, origin, spacing, name):
        volumeNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScalarVolumeNode", name)
        
        # Slicer espera (Z, Y, X)
        volumeNode.SetOrigin(origin)
        volumeNode.SetSpacing(spacing)
        slicer.util.updateVolumeFromArray(volumeNode, array)
        
        # Asegurar que el origen/espaciado se mantengan (updateVolumeFromArray puede resetearlos)
        volumeNode.SetOrigin(origin)
        volumeNode.SetSpacing(spacing)
        
        print(f"--- Volumen Creado: {name} ---")
        print(f"  Array shape: {array.shape}")
        print(f"  Origin: {volumeNode.GetOrigin()}")
        print(f"  Spacing: {volumeNode.GetSpacing()}")
        
        # --- Configuración de Visualización ---
        volumeNode.CreateDefaultDisplayNodes()
        displayNode = volumeNode.GetDisplayNode()
        
        if displayNode:
            # 1. Usar escala de colores tipo Rainbow
            # Intentamos obtener el nodo de color Rainbow de forma robusta
            rainbowColorNode = slicer.mrmlScene.GetNodeByID("vtkMRMLColorTableNodeFileRainbow.txt")
            if not rainbowColorNode:
                # Fallback a búsqueda por nombre
                rainbowColorNode = slicer.util.getNode("Rainbow")
            
            if rainbowColorNode:
                # Intentamos varios métodos según la versión de Slicer
                # SetAndObserveColorNodeID es el estándar en Slicer 5.x
                colorNodeID = rainbowColorNode.GetID()
                if hasattr(displayNode, 'SetAndObserveColorNodeID'):
                    displayNode.SetAndObserveColorNodeID(colorNodeID)
                elif hasattr(displayNode, 'SetAndUseColorNodeID'):
                    displayNode.SetAndUseColorNodeID(colorNodeID)
                elif hasattr(displayNode, 'SetColorNodeID'):
                    displayNode.SetColorNodeID(colorNodeID)
            else:
                logging.warning(f"No se pudo encontrar la escala de colores Rainbow para {name}")
            
            # 2. Ajustar automáticamente el brillo y contraste
            displayNode.AutoWindowLevelOn()
            
            # 3. Hacer que los valores bajos sean transparentes (umbral relativo al máximo)
            displayNode.SetThreshold(np.max(array) * 0.01 if np.max(array) > 0 else 0, np.max(array))
            displayNode.SetApplyThreshold(True)
        
        # Mostrar en los visores y ajustar vista
        slicer.util.setSliceViewerLayers(background=volumeNode, fit=True)
        
        logging.info(f"Volumen '{name}' creado.")
        return volumeNode

    def _cleanupOldVolumes(self):
        prefixes = [
            "Boron_Kerma", "FastNeutron_Kerma", "ThermalNeutron_Kerma",
            "Gamma_Photon", "Raw_Boro_PHITS", "Raw_Fast_PHITS",
            "Raw_Thermal_PHITS", "Raw_Gamma_PHITS", "Dosis_Final_BNCT",
            "Dosis_Fisica_BNCT", "Dosis_RBE_BNCT", "Dosis_IsoE_BNCT"
        ]
        node_col = slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode")
        for i in range(node_col.GetNumberOfItems()):
            node = node_col.GetItemAsObject(i)
            if node and any(node.GetName().startswith(p) for p in prefixes):
                slicer.mrmlScene.RemoveNode(node)

    def _organizeInSubjectHierarchy(self, field_count, has_gamma):
        shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
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
