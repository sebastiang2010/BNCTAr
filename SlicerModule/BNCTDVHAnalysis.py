import os, json
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy as np

_ORGAN_CONFIG_PATH = os.path.join(os.path.dirname(__file__), '..', 'json', 'organ_config.jsonc')
_DEBUG_LOG_PATH = os.path.join(os.path.dirname(__file__), '_dvh_debug.txt')

def _log(msg):
    with open(_DEBUG_LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')

def _load_jsonc(path):
    import re
    with open(path) as f:
        text = f.read()
    text = re.sub(r'//.*?\n|/\*.*?\*/', '', text, flags=re.S)
    return json.loads(text)

class BNCTDVHAnalysis(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "BNCT DVH Analysis"
        self.parent.categories = ["BNCT"]
        self.parent.dependencies = []
        self.parent.contributors = ["Antigravity AI / BNCTAr Team"]
        self.parent.helpText = "Modulo para DVH BNCT."
        self.parent.acknowledgementText = "BNCTAr en 3D Slicer."

class BNCTDVHAnalysisWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        try: os.remove(_DEBUG_LOG_PATH)
        except: pass

        # Limpiar nodos viejos (con proteccion por si las clases VTK aun no se registraron)
        try:
            for cls, prefix in [("vtkMRMLPlotSeriesNode", "DVH_"),
                                 ("vtkMRMLTableNode", "Table_DVH_")]:
                coll = slicer.mrmlScene.GetNodesByClass(cls)
                for i in range(coll.GetNumberOfItems()):
                    node = coll.GetItemAsObject(i)
                    if node and node.GetName().startswith(prefix):
                        slicer.mrmlScene.RemoveNode(node)
            old_chart = slicer.mrmlScene.GetFirstNodeByName("Chart_DVH")
            if old_chart:
                slicer.mrmlScene.RemoveNode(old_chart)
        except:
            pass

        # --- UI ---
        dvhCollapsibleButton = ctk.ctkCollapsibleButton()
        dvhCollapsibleButton.text = "Analisis de DVH"
        dvhCollapsibleButton.collapsed = True
        self.layout.addWidget(dvhCollapsibleButton)
        dvhFormLayout = qt.QFormLayout(dvhCollapsibleButton)

        organ_widget = qt.QWidget()
        organ_wl = qt.QVBoxLayout(organ_widget)
        organ_wl.setContentsMargins(0, 0, 0, 0)
        self.dvhTodosCheck = qt.QCheckBox("Todos")
        self.dvhTodosCheck.checked = True
        organ_wl.addWidget(self.dvhTodosCheck)
        scroll = qt.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(180)
        sw = qt.QWidget()
        self.dvhScrollLayout = qt.QVBoxLayout(sw)
        self.dvhScrollLayout.setSpacing(2)
        scroll.setWidget(sw)
        organ_wl.addWidget(scroll)
        dvhFormLayout.addRow("Organos:", organ_widget)
        self.dvhTodosCheck.connect('toggled(bool)', self._onDvhTodosToggled)

        self.dvhDoseTypeCombo = qt.QComboBox()
        self.dvhDoseTypeCombo.addItems(["Dosis Fisica", "Dosis RBE", "Dosis Isoefectiva"])
        dvhFormLayout.addRow("Tipo de Dosis:", self.dvhDoseTypeCombo)

        self.dvhModeloCombo = qt.QComboBox()
        self.dvhModeloCombo.setMinimumWidth(300)
        self.dvhModeloCombo.visible = False
        dvhFormLayout.addRow("Modelo:", self.dvhModeloCombo)

        self.dvhB10Spin = qt.QDoubleSpinBox()
        self.dvhB10Spin.setRange(0, 1000)
        self.dvhB10Spin.value = 15.0
        dvhFormLayout.addRow("Conc. Boro (ppm):", self.dvhB10Spin)

        self.dvhTirSpin = qt.QDoubleSpinBox()
        self.dvhTirSpin.setRange(0, 1000)
        self.dvhTirSpin.value = 1.0
        dvhFormLayout.addRow("TIR (min):", self.dvhTirSpin)

        self.maskVolumeSelector = slicer.qMRMLNodeComboBox()
        self.maskVolumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.maskVolumeSelector.selectNodeUponCreation = True
        self.maskVolumeSelector.addEnabled = False
        self.maskVolumeSelector.removeEnabled = False
        self.maskVolumeSelector.noneEnabled = False
        self.maskVolumeSelector.showHidden = False
        self.maskVolumeSelector.setMRMLScene(slicer.mrmlScene)
        dvhFormLayout.addRow("Mascara ROIs:", self.maskVolumeSelector)

        self.dvhButton = qt.QPushButton("Calcular DVH")
        self.dvhButton.enabled = True
        dvhFormLayout.addRow(self.dvhButton)

        self.dvhButton.connect('clicked(bool)', self.onDVHButton)
        self.dvhDoseTypeCombo.connect('currentIndexChanged(int)', self._onDvhDoseTypeChanged)

        self._dvh_models = self._loadModels()
        self.dvhModeloCombo.clear()
        for m in self._dvh_models:
            self.dvhModeloCombo.addItem(f"{m['id']}. {m['name']}")

        self._dvh_config = self._loadConfigForDvh()
        self.layout.addStretch(1)

    def cleanup(self):
        ScriptedLoadableModuleWidget.cleanup(self)

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
            if w: w.delete()
        for row in config:
            cb = qt.QCheckBox(f"{row['name']}  [mesh {int(row['nmesh'])}]")
            cb.checked = True
            self.dvhScrollLayout.addWidget(cb)
        self.dvhTodosCheck.blockSignals(True)
        self.dvhTodosCheck.checked = True
        self.dvhTodosCheck.blockSignals(False)

    def _onDvhDoseTypeChanged(self, idx):
        self.dvhModeloCombo.visible = (idx == 2)

    def _loadModels(self):
        path = os.path.join(os.path.dirname(__file__), '..', 'json', 'models.jsonc')
        if os.path.exists(path):
            try:
                return _load_jsonc(path).get('isoeffective', [])
            except Exception as e:
                logging.warning(f"No se pudo cargar models.jsonc: {e}")
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
            model_params = [m['alpha_r'], m['beta_r'], m['alpha_boro'], m['beta_boro'],
                            m['alpha_thn'], m['beta_thn'], m['alpha_fast'], m['beta_fast'],
                            m['tof'], m['tos'], m['pf_g'], m['ps_g'], m['pf_bnct'], m['ps_bnct']]

        mask_node = self.maskVolumeSelector.currentNode()
        if not mask_node:
            slicer.util.errorDisplay("Seleccione un volumen de mascara.")
            return

        logic = BNCTDVHAnalysisLogic()
        try:
            qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
            logic.calculateDVH(organs=selected, dose_type=dose_type,
                               b10_conc=b10, tir=tir,
                               model_params=model_params, mask_volume=mask_node)
        except Exception as e:
            slicer.util.errorDisplay(f"Error en DVH: {str(e)}")
            import traceback
            logging.error(traceback.format_exc())
            _log(f"ERROR: {e}\n{traceback.format_exc()}")
        finally:
            qt.QApplication.restoreOverrideCursor()

    def _loadConfigForDvh(self):
        if os.path.exists(_ORGAN_CONFIG_PATH):
            try:
                data = _load_jsonc(_ORGAN_CONFIG_PATH)
                config = data.get('rois', [])
                if config:
                    self._populateDvhCheckboxes(config)
                    return config
            except Exception as e:
                logging.warning(f"No se pudo cargar organ_config.jsonc: {e}")
        mat_path = os.path.join(os.path.dirname(__file__), 'mat', 'info_rois.mat')
        if os.path.exists(mat_path):
            try:
                logic = BNCTDVHAnalysisLogic()
                config = logic.loadConfigFromMat(mat_path)
                d = os.path.dirname(_ORGAN_CONFIG_PATH)
                if not os.path.exists(d):
                    os.makedirs(d)
                with open(_ORGAN_CONFIG_PATH, 'w', encoding='utf-8') as f:
                    json.dump({'rois': config}, f, indent=2, ensure_ascii=False)
                self._populateDvhCheckboxes(config)
                return config
            except Exception as e:
                logging.warning(f"Fallback MAT fallo: {e}")
        return []

class BNCTDVHAnalysisLogic(ScriptedLoadableModuleLogic):

    def _cleanupOldDVH(self):
        for cls, prefix in [("vtkMRMLPlotSeriesNode", "DVH_"),
                             ("vtkMRMLTableNode", "Table_DVH_")]:
            coll = slicer.mrmlScene.GetNodesByClass(cls)
            for i in range(coll.GetNumberOfItems()):
                node = coll.GetItemAsObject(i)
                if node and node.GetName().startswith(prefix):
                    slicer.mrmlScene.RemoveNode(node)
        old = slicer.mrmlScene.GetFirstNodeByName("Chart_DVH")
        if old:
            slicer.mrmlScene.RemoveNode(old)

    def _getKermaVolumesForMesh(self, mesh_number):
        import re
        volumes = []
        for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
            name = node.GetName()
            if name == 'Gamma_Photon':
                volumes.append((node, name))
                continue
            if 'Kerma' in name:
                match = re.search(r'Kerma(\d+)$', name)
                if match and int(match.group(1)) == mesh_number:
                    volumes.append((node, name))
        return volumes

    def _calcDoseForMesh(self, nmesh, dose_type, b10, tir, ratio,
                          rbe_b, rbe_f, rbe_t, rbe_g, model_params):
        volumes = self._getKermaVolumesForMesh(nmesh)
        if not volumes:
            return None
        comp = {'boro': None, 'fast': None, 'thermal': None, 'gamma': None}
        for node, vol_name in volumes:
            arr = slicer.util.arrayFromVolume(node).astype(np.float64)
            if 'Boron' in vol_name or 'boro' in vol_name.lower():
                comp['boro'] = arr if comp['boro'] is None else comp['boro'] + arr
            elif 'Fast' in vol_name:
                comp['fast'] = arr if comp['fast'] is None else comp['fast'] + arr
            elif 'Thermal' in vol_name or 'thn' in vol_name.lower():
                comp['thermal'] = arr if comp['thermal'] is None else comp['thermal'] + arr
            elif 'Gamma' in vol_name or 'photon' in vol_name.lower():
                comp['gamma'] = arr if comp['gamma'] is None else comp['gamma'] + arr
        if comp['boro'] is None:
            return None
        if dose_type == "Dosis Fisica":
            return ((comp['boro'] * b10 * ratio +
                     (comp['fast'] if comp['fast'] is not None else 0) +
                     (comp['thermal'] if comp['thermal'] is not None else 0) +
                     (comp['gamma'] if comp['gamma'] is not None else 0)) * tir)
        elif dose_type == "Dosis RBE":
            return ((comp['boro'] * b10 * ratio * rbe_b +
                     (comp['fast'] if comp['fast'] is not None else 0) * rbe_f +
                     (comp['thermal'] if comp['thermal'] is not None else 0) * rbe_t +
                     (comp['gamma'] if comp['gamma'] is not None else 0) * rbe_g) * tir)
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
        D_boro = comp['boro'] * ratio
        D_thn = comp['thermal'] if comp['thermal'] is not None else np.zeros_like(D_boro)
        D_fast = comp['fast'] if comp['fast'] is not None else np.zeros_like(D_boro)
        D_g = comp['gamma'] if comp['gamma'] is not None else np.zeros_like(D_boro)
        Dt = D_boro + D_thn + D_fast + D_g
        with np.errstate(divide='ignore', invalid='ignore'):
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
        with np.errstate(divide='ignore', invalid='ignore'):
            a_12 = np.divide(fb, fb+fn, out=np.zeros_like(fb), where=(fb+fn)>0)
            a_21 = np.divide(fn, fb+fn, out=np.zeros_like(fn), where=(fb+fn)>0)
            a_13 = np.divide(fb, fb+fg, out=np.zeros_like(fb), where=(fb+fg)>0)
            a_31 = np.divide(fg, fb+fg, out=np.zeros_like(fg), where=(fb+fg)>0)
            a_23 = np.divide(fn, fn+fg, out=np.zeros_like(fn), where=(fn+fg)>0)
            a_32 = np.divide(fg, fn+fg, out=np.zeros_like(fg), where=(fn+fg)>0)
        G_12 = lc_g(tof, tos, a_12*pf_bnct + a_21*pf_bnct, 0, tir)
        G_23 = lc_g(tof, tos, a_23*pf_bnct + a_32*pf_g, 0, tir)
        G_31 = lc_g(tof, tos, a_31*pf_g + a_13*pf_bnct, 0, tir)
        D_lineal = (alpha_boro * D_boro + alpha_thn * D_thn +
                    alpha_fast * D_fast + alpha_r * D_g) * tir
        D_cuad = (G_11 * beta_boro * D_boro**2 + G_22 * beta_thn * D_thn**2 +
                  G_22 * beta_fast * D_fast**2 + G_33 * beta_r * D_g**2 +
                  2*G_12*np.sqrt(beta_boro*beta_thn)*D_boro*D_thn +
                  2*G_12*np.sqrt(beta_boro*beta_fast)*D_boro*D_fast +
                  2*G_31*np.sqrt(beta_boro*beta_r)*D_boro*D_g +
                  2*G_11*np.sqrt(beta_thn*beta_fast)*D_thn*D_fast +
                  2*G_23*np.sqrt(beta_thn*beta_r)*D_thn*D_g +
                  2*G_23*np.sqrt(beta_fast*beta_r)*D_fast*D_g) * tir**2
        dose = 0.5 * (alphabeta_r / Gr) * (
            np.sqrt(1 + 4 * Gr * beta_r / alpha_r**2 * (D_lineal + D_cuad)) - 1)
        return np.nan_to_num(dose)

    def calculateDVH(self, organs, dose_type, b10_conc, tir, model_params=None, mask_volume=None):
        import re
        _log(f"\n=== DVH START: {len(organs)} organs, {dose_type} ===")

        if mask_volume is None:
            raise ValueError("Debe seleccionar un volumen de mascara.")

        # Encontrar referencia
        ref_node = None
        for node in slicer.mrmlScene.GetNodesByClass("vtkMRMLScalarVolumeNode"):
            name = node.GetName()
            if any(x in name for x in ['Kerma', 'Gamma', 'Dosis', 'DVH', 'Temp']):
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

        # Limpiar nodos viejos
        self._cleanupOldDVH()
        plotChartNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode", "Chart_DVH")

        all_names = []
        any_dose_computed = False
        global_max_dose = 0.0
        spacing = ref_node.GetSpacing()
        voxel_vol_cc = spacing[0] * spacing[1] * spacing[2] / 1000.0

        for org in organs:
            nmesh = int(org['nmesh'])
            ratio = org['ratio']
            rbe_b = org['rbe_boron']
            rbe_f = org['rbe_fast']
            rbe_t = org['rbe_thermal']
            rbe_g = org['rbe_gamma']

            print(f"[DVH] Procesando {org['name']} (mesh {nmesh}, gray_level {org.get('gray_level',0)})...")
            dose = self._calcDoseForMesh(nmesh, dose_type, b10_conc, tir,
                                          ratio, rbe_b, rbe_f, rbe_t, rbe_g,
                                          model_params)
            if dose is None:
                print(f"[DVH]  -> SKIP: no hay kerma volumes para mesh {nmesh}")
                continue
            if dose.shape != ct_shape:
                try:
                    from scipy.ndimage import zoom
                    factors = [max(1, ct_shape[0]) / max(1, dose.shape[0]),
                               max(1, ct_shape[1]) / max(1, dose.shape[1]),
                               max(1, ct_shape[2]) / max(1, dose.shape[2])]
                    dose = zoom(dose, factors, order=1)
                except Exception as e:
                    print(f"[DVH]  -> SKIP: zoom fallo - {e}")
                    continue
            any_dose_computed = True

            gl = org.get('gray_level', 0)
            mask = (mask_array == gl)
            roi_doses = dose[mask]
            if len(roi_doses) == 0:
                print(f"[DVH]  -> SKIP: gray_level={gl} sin voxeles en mascara para '{org['name']}'")
                continue

            roi_doses = roi_doses[~np.isnan(roi_doses)]
            if len(roi_doses) == 0:
                print(f"[DVH]  -> SKIP: todos NaN para '{org['name']}'")
                continue

            sorted_doses = np.sort(roi_doses, axis=None)
            n = len(sorted_doses)
            dvh_vol = np.linspace(100.0, 0.0, n)

            # Agregar punto (0, 100) al inicio para que la curva arranque desde arriba
            sorted_doses = np.insert(sorted_doses, 0, 0.0)
            dvh_vol = np.insert(dvh_vol, 0, 100.0)
            n = len(sorted_doses)

            vol_cc = voxel_vol_cc * (n - 1)  # voxeles reales sin contar el punto (0,100)
            Dmax = float(sorted_doses[-1])
            Dmin = float(sorted_doses[0])
            Dmean = float(np.mean(roi_doses))
            if Dmax > global_max_dose:
                global_max_dose = Dmax

            _log(f"  OK {org['name']}: {n} vox, vol={vol_cc:.4f}cc, "
                 f"Dmin={Dmin:.6e} Dmax={Dmax:.6e} Dmedia={Dmean:.6e}")
            _log(f"    primeras 5 dosis: {sorted_doses[:5].tolist()}")
            _log(f"    ultimas 5 dosis: {sorted_doses[-5:].tolist()}")

            safe_name = re.sub(r'[^a-zA-Z0-9]', '_', org['name'])
            all_names.append(org['name'])

            # --- Crear tabla con formato correcto (array 2D) ---
            tn = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", f"Table_DVH_{safe_name}")
            table_array = np.column_stack((sorted_doses, dvh_vol))
            slicer.util.updateTableFromArray(tn, table_array, ["Dose", "VolumePct"])

            # Forzar nombres de columna y notificar cambio
            tn.GetTable().GetColumn(0).SetName("Dose")
            tn.GetTable().GetColumn(1).SetName("VolumePct")
            tn.Modified()

            print(f"[DVH] Tabla '{tn.GetName()}': {tn.GetTable().GetNumberOfRows()} filas")
            print(f"[DVH] Rango dosis: {sorted_doses[0]:.6e} - {sorted_doses[-1]:.6e} Gy")

            # --- Crear serie del gráfico como DISPERSIÓN (PlotTypeScatter) ---
            ps = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", f"DVH_{safe_name}")
            ps.SetAndObserveTableNodeID(tn.GetID())
            ps.SetXColumnName("Dose")
            ps.SetYColumnName("VolumePct")
            # Usar tipo Scatter (1) para que respete los valores de columna
            ps.SetPlotType(slicer.vtkMRMLPlotSeriesNode.PlotTypeScatter)
            # Opcional: mostrar como línea sólida y puntos
            ps.SetLineStyle(slicer.vtkMRMLPlotSeriesNode.LineStyleSolid)
            ps.SetMarkerStyle(slicer.vtkMRMLPlotSeriesNode.MarkerStyleNone)  # sin puntos
            ps.SetUniqueColor()
            plotChartNode.AddAndObservePlotSeriesNodeID(ps.GetID())

        if not any_dose_computed:
            raise ValueError("No se pudo calcular dosis para ningun organo.")
        if not all_names:
            raise ValueError("No se encontraron voxeles para ningun organo.")

        _log(f"  global_max_dose = {global_max_dose:.6e}")
        _log(f"  Nombres: {all_names}")

        # Configurar propiedades del gráfico
        dose_tag = dose_type.replace("Dosis ", "")
        plotChartNode.SetTitle(f"DVH - {', '.join(all_names[:3])}" +
                               ("..." if len(all_names) > 3 else f" ({dose_tag})"))
        plotChartNode.SetXAxisTitle("Dosis [Gy]")
        plotChartNode.SetYAxisTitle("Volumen [%]")
        plotChartNode.SetLegendVisibility(True)

        # Forzar rangos con los valores reales de dosis (aunque sean pequeños)
        xmin = 0.0
        xmax = global_max_dose * 1.05 if global_max_dose > 0 else 1.0
        plotChartNode.SetXAxisRange(xmin, xmax)
        plotChartNode.SetYAxisRange(0, 105)
        _log(f"  Axis ranges: X=[{xmin}, {xmax:.6e}] Y=[0, 105]")

        # Mostrar en el layout
        layoutManager = slicer.app.layoutManager()
        layoutWithPlot = slicer.modules.plots.logic().GetLayoutWithPlot(layoutManager.layout)
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

        # Mostrar resumen al usuario
        slicer.util.infoDisplay(
            f"DVH calculado: {len(all_names)} organos\n"
            f"Dosis max: {global_max_dose:.4e} Gy\n"
            f"Ver tabla _dvh_debug.txt para datos crudos",
            windowTitle="DVH Completo")

    def loadConfigFromMat(self, path):
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
                    'rbe_boron': float(read_ref(rbe_group['boron'][0,0])),
                    'rbe_fast': float(read_ref(rbe_group['fast_neutron'][0,0])),
                    'rbe_thermal': float(read_ref(rbe_group['nitrogen'][0,0])),
                    'rbe_gamma': float(read_ref(rbe_group['photon'][0,0])),
                }
            return {'rbe_boron': 1.3, 'rbe_fast': 3.2, 'rbe_thermal': 3.2, 'rbe_gamma': 1.0}
        with h5py.File(path, 'r') as f:
            info = f['info_dvh_cnea']
            num_rows = len(info['nroi'])
            for i in range(num_rows):
                rbe_grp = info['RBE'][i,0]
                rbe_vals = read_rbe(rbe_grp) if isinstance(rbe_grp, h5py.Group) else read_rbe(rbe_grp)
                boron_val = 15.0
                if 'boron' in info:
                    boron_val = float(read_ref(info['boron'][i,0]))
                row = {
                    'nroi': int(float(read_ref(info['nroi'][i,0]))),
                    'nmesh': int(float(read_ref(info['nmesh'][i,0]))),
                    'name': str(read_ref(info['roi_name'][i,0])),
                    'rbe_boron': float(rbe_vals['rbe_boron']),
                    'rbe_fast': float(rbe_vals['rbe_fast']),
                    'rbe_thermal': float(rbe_vals['rbe_thermal']),
                    'rbe_gamma': float(rbe_vals['rbe_gamma']),
                    'ratio': float(read_ref(info['ratio'][i,0])),
                    'tir': float(read_ref(info['tir'][i,0])),
                    'tumor': bool(int(float(read_ref(info['tumor'][i,0])))),
                    'gray_level': int(float(read_ref(info['gray_level'][i,0]))),
                }
                config.append(row)
        return config