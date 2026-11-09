# AGENTS.md — BNCTAr Dosimetría

Reglas del proyecto. Aplican a todo análisis, código y reporte generado en este repo
(`C:\programas\BNCTAr\BNCTAr-dosimetria`).

> Nota: `SlicerModule/AGENTS.md` contiene reglas ajenas a este proyecto (bot de trading
> Binance/HFT). **No aplican** a este repo. Este archivo es la fuente de reglas raíz.

---

## 0. Fuente externa de SOLO LECTURA — UTCP

`C:\Seba\BNCT\BNCT-clinica\mat\UTCP` es una **fuente externa de solo lectura**.

**PROHIBIDO** dentro de esa carpeta:

- Crear, editar, renombrar, mover o borrar cualquier archivo, **incluido su `agents.md`**.
- Correr scripts que escriban ahí: `save`, `writetable`, `exportgraphics`, `print`,
  `fopen(...,'w')`, `delete`, `movefile`, `copyfile` con destino dentro de UTCP.
- Ejecutar scripts *desde* UTCP: `B2_v30.m` y `B3_v15.m` usan `addpath(fullfile(pwd,'utils'))`
  y rutas relativas; correrlos ahí genera archivos de salida en `resultados/` y `reporte/`
  de UTCP. **No se ejecutan en su ubicación original.**

**Permitido:** lectura (`load`, `readmatrix`, `imread`, `dir`, `fopen(...,'r')`)
para inspeccionar datos o contrastar resultados.

**Regla de oro:** si algo de UTCP hace falta, se **copia** a este proyecto y se trabaja
sobre la copia. El original queda intacto como referencia auditable.

### 0.1 Dónde va lo copiado

| Contenido copiado de UTCP | Destino en este repo |
|---|---|
| Código radiobiológico (`ntcp_isoE.m`, `isoE_params.m`, `recalc_isoE_ROI.m`, `dvh_isoE.m`, `isoE_tabla_organos.m`, `ntcp_skin_riback_v1.m`) | `utcp/` (adaptado, no copia literal) |
| Parámetros y catálogo de modelos | `model_list.txt` (extendido) + `utcp/utcp_ntcp_registry.m` |
| Fixtures de contraste (`resultados/dvh_B2_v30.mat`) | `utcp/test/golden/` — **solo si se necesita test de regresión** |
| Mallas (`B2/c1/mesh_*.mat`, `B3/…`) | **NO se copian** (4,3 GB). Este proyecto lee PHITS/MCNP directamente. |
| Reportes, póster, papers, `old/`, `imagenes/`, `prims/` | **NO se copian** |

---

## 1. Reglas radiobiológicas heredadas de UTCP (obligatorias)

1. **Tumor = `GTV70`, no `GTV 3mm`.** Priorizar por nombre `GTV70`; fallback a
   `patient.roi_tumor == 1`. Nunca por índice fijo (`idx_tumor = 2` está prohibido:
   en la tabla real el índice 2 es `RT LENS`).
2. **La marca de tumor es manual y explícita** (columna `Tumor` de `tbl_roi`).
   Debe haber **exactamente un** ROI marcado como tumor para calcular TCP/UTCP.
   Si hay más de uno, gana `GTV70`; si hay cero, error con diálogo.
3. **Restricción principal en B2 = PIEL (SKIN)**, no cerebro. El haz de B2 es muy
   térmico; la dosis limitante es la piel. Reflejarlo en Métodos y Resultados.
4. **Segmentación de piel:** erodar `BODY` con `GROSOR_PIEL = 1` vóxel (isótropo,
   6 vecinos). `SKIN` es un ROI sintético, **no** figura entre los 20 ROIs del DVH.
5. **Relleno de máscara con `for` por corte (NO `imfill` 3D).** Origen literal del
   usuario: *"NO puedes sacar el for, esta mal ya lo dije muchas veces, imfill no
   funciona 3D"*.
6. **Discrepancias RBE al reportar** (pipeline RBE vs scripts históricos):
   boro `1.30 / 3.80 / 1.50` vs `2.5`; fast/N `3.20` vs `3.0`; $\gamma$ `1.0` en ambos.

---

## 2. Contrato de datos del módulo UTCP (`utcp/`)

- **Entrada canónica:** un elemento de `R` devuelto por `f_dvh.m`:
  `R(i).name`, `R(i).dvh_dose` (ascendente), `R(i).dvh_vol` (% acumulado),
  `R(i).vol_cc`, `R(i).tumor`.
- **Unidades:** los modelos NTCP/TCP de UTCP esperan **Gy isoE**, no Gy físicos ni
  Gy (RBE). Todo cálculo radiobiológico debe partir del DVH de `Diso`
  (`patient.show == 3`). Si el DVH guardado viene de `Df` o `Dw`, abortar con error
  explícito.
- **Conversión:** `utcp_dvh_differential(R)` produce `[dose, dvol]` normalizado
  ($\sum \Delta v = 1$), deduplicando dosis repetidas (si no, `interp1` falla).
- **Nada de estado global nuevo:** persistir la selección de modelo en
  `patient.ntcp_model` / `patient.utcp_result` (compatible con `btn_save_case`).

$$\mathrm{UTCP} = \mathrm{TCP}_{\mathrm{tumor}} \cdot \prod_{j \in \text{OAR}} \left(1 - \mathrm{NTCP}_j\right)$$

---

## 3. Reglas de código MATLAB

1. **Versionado obligatorio:** al modificar código, crear archivo nuevo con versión
   incrementada (`f_dvh_v2.m`, `dosimetria_hn_v2.m`) o commit descriptivo en git.
   No sobrescribir sin historial.
2. **Gate R-CHK:** antes de entregar, `checkcode('archivo.m')` con **0 errores**.
   Cada warning se justifica individualmente. `checkcode` es solo lectura.
3. **Verificación con MATLAB MCP:** `matlab_check_matlab_code`,
   `matlab_evaluate_matlab_code`, `matlab_run_matlab_file`,
   `matlab_run_matlab_test_file`, `matlab_detect_matlab_toolboxes`.
4. **Nada suelto:** tests → `test/`; utilidades → `utils/`; resultados `.mat` →
   `resultados/`; reportes → `reporte/`. No dejar `.m`, `.mat`, `.png`, `.log` en la raíz.
5. **Código muerto conocido — no usar como base:** `ntcp2.m` (usa `a_12…a_13`
   indefinidas), `f_evaluacion.m`, `f_error.m`, `f_analisis_mesh.m`
   (`struc`, `niveles`, `nu2srt` indefinidas), `dosis_fisica.m` (stub).
6. **`.mlapp` vs `.m` divergen.** Son ramas paralelas con lógica distinta
   (`dosimetria_hn_App.mlapp` tiene el DVH viejo sin `f_dvh`). Elegir una y declararlo;
   no asumir paridad.

---

## 4. Fórmulas matemáticas (LaTeX)

En archivos `.md` de este proyecto, toda fórmula va en LaTeX:

- `$$...$$` bloque, `$...$` inline. Nunca Unicode plano.
- Identificadores multi-letra con `\mathrm{}`: `\mathrm{NTCP}`, `D_{\mathrm{boro}}`,
  `\mathrm{RBE}_{\mathrm{boro}}`. Solo subíndices de una letra o griega quedan en cursiva.
- `\left(` `\right(` autoajustables, `\cdot` para multiplicar, `\,` antes de unidades.
- Unidad con calificador separada: `Gy (RBE)`, `Gy (IsoE)`, `Gy (EQD2)` —
  nunca `Gy(RBE)`.
- Fórmulas largas con `aligned` dentro de `$$...$$`.

---

## 5. Seguridad y datos

- Nunca subir: `.env`, API keys, credenciales, logs, `__pycache__`, archivos temporales.
- `.mat` de pacientes no se commitean si superan unos pocos MB; usar `resultados/` ignorado.
- Ignorados por git: `graphify-out/`, `ruvector.db`, `.foxcode/`, salidas PHITS.
