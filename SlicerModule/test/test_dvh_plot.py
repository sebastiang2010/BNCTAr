"""
Test standalone: verifica que la creacion de tabla y plot series
tengan los mismos datos (simula lo que hace calculateDVH).
Ejecutar en Slicer Python console:
  exec(open('test_dvh_plot.py').read())
"""
import numpy as np
import slicer

print("=" * 60)
print("TEST: Verificacion tabla vs plot en DVH")
print("=" * 60)

# --- 1. Generar datos simulados (como los reales del usuario) ---
np.random.seed(42)
n_voxels = 100
dosis_reales = np.random.uniform(0.5, 4.0, n_voxels)  # Gy, valores realistas
dosis_reales.sort()
vol_pct = 100.0 * (1.0 - np.arange(n_voxels) / n_voxels)

print(f"\nDatos de entrada:")
print(f"  Voxeles: {n_voxels}")
print(f"  Dosis min: {dosis_reales[0]:.4f}  max: {dosis_reales[-1]:.4f}  media: {dosis_reales.mean():.4f}")
print(f"  Vol% primer: {vol_pct[0]:.2f}  ultimo: {vol_pct[-1]:.2f}")

# --- 2. Crear tabla (IGUAL que en calculateDVH) ---
tn = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", "_Test_DVH_Table")
slicer.util.updateTableFromArray(tn, [dosis_reales, vol_pct], ["Dose", "VolumePct"])
print(f"\nTabla creada: {tn.GetName()}")
print(f"  Columnas: {tn.GetNumberOfColumns()}")

# --- 3. Verificar columnas ---
col0 = tn.GetColumn(0)  # Dose
col1 = tn.GetColumn(1)  # VolumePct
print(f"  Col0 nombre: '{col0.GetName()}'  valores={col0.GetNumberOfValues()}")
print(f"  Col1 nombre: '{col1.GetName()}'  valores={col1.GetNumberOfValues()}")

# --- 4. Leer datos de la tabla ---
dosis_tabla = np.array([col0.GetValue(i) for i in range(col0.GetNumberOfValues())])
vol_tabla = np.array([col1.GetValue(i) for i in range(col1.GetNumberOfValues())])

# --- 5. Crear plot series (IGUAL que en calculateDVH) ---
ps = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotSeriesNode", "_Test_DVH_Series")
ps.SetAndObserveTableNodeID(tn.GetID())
ps.SetXColumnName("Dose")
ps.SetYColumnName("VolumePct")
ps.SetPlotType(0)  # PlotTypeLine

# --- 6. VERIFICAR: comparar datos originales vs tabla ---
print(f"\n--- VERIFICACION ---")

# Verificar que dosis coinciden
diff_dosis = np.max(np.abs(dosis_reales - dosis_tabla))
diff_vol = np.max(np.abs(vol_pct - vol_tabla))

print(f"  Diferencia max en Dosis: {diff_dosis:.15f}")
print(f"  Diferencia max en Vol%%: {diff_vol:.15f}")

if diff_dosis < 1e-12 and diff_vol < 1e-12:
    print(f"\n  ✅ PASS: Tabla y datos ORIGINALES coinciden exactamente")
else:
    print(f"\n  ❌ FAIL: Hay diferencias entre tabla y datos originales")

# --- 7. Verificar que el plot series apunta a las columnas correctas ---
xcol = ps.GetXColumnName()
ycol = ps.GetYColumnName()
print(f"\n  PlotSeries X column: '{xcol}'  Y column: '{ycol}'")

if xcol == "Dose" and ycol == "VolumePct":
    print(f"  ✅ PASS: PlotSeries apunta a columnas correctas")
else:
    print(f"  ❌ FAIL: Columnas del plot no coinciden")

# --- 8. Crear chart y mostrar ---
chart = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLPlotChartNode", "_Test_DVH_Chart")
chart.AddAndObservePlotSeriesNodeID(ps.GetID())
chart.SetTitle("TEST DVH - Verificacion")
chart.SetXAxisTitle("Dosis [Gy]")
chart.SetYAxisTitle("Volumen [%]")
chart.SetLegendVisibility(False)

# Mostrar en layout
layoutManager = slicer.app.layoutManager()
layoutWithPlot = slicer.modules.plots.logic().GetLayoutWithPlot(layoutManager.layout)
layoutManager.setLayout(layoutWithPlot)
plotWidget = layoutManager.plotWidget(0)
if plotWidget:
    pvn = plotWidget.mrmlPlotViewNode()
    if pvn:
        pvn.SetPlotChartNodeID(chart.GetID())

print(f"\n  ✅ Chart _Test_DVH_Chart creado y mostrado en layout")
print(f"  📊 Miralo en Slicer - deberia ver curva descendiendo de 100% a 0%")
print(f"  🎯 Dosis en eje X: 0 a ~4.5 Gy  (valores REALES, no 1e-16)")
print("=" * 60)

# Limpiar al final (descomentar si se quiere)
# slicer.mrmlScene.RemoveNode(tn)
# slicer.mrmlScene.RemoveNode(ps)
# slicer.mrmlScene.RemoveNode(chart)
