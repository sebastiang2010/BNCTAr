# Reporte — Implementación de TCP, NTCP y UTCP en el módulo BNCT de 3D Slicer

Fecha: 2026-11-09. Implementa la spec completa (fases 0–5). Motor radiobiológico
puro Python + numpy, sin Qt ni slicer (testeable fuera de 3D Slicer).

## 1. Archivos creados y modificados

### Nuevos — `SlicerModule/`
| Archivo | Contenido |
|---|---|
| `utcp_models.py` | Registro de modelos: TCP Provenzano de referencia (c1=22, c2=0.33, $\alpha=0.34\,\mathrm{Gy}^{-1}$, $\beta=8.3\times10^{-9}\,\mathrm{Gy}^{-2}$, fracciones=1), LKB Mucosa/Cerebro, LKB configurable, Riback (bloqueado), defaults isoE por tejido, validadores |
| `utcp_dvh.py` | DVH acumulado→diferencial (dedupe + normalización, fiel a `ntcp_isoE.m`), EUD, Dx %, VDy, `voxel_volume_cc` |
| `utcp_tcp.py` | Provenzano: $H_i=\exp(-A_r D_i - G_r B_r D_i^2)^{n_f/c_2}$, $\mathrm{TCP}=\exp\!\left(-c_1\left(\sum_i H_i\,dv_i\right)^{c_2}\right)$, PEUD, $G_r$ de Catcheside |
| `utcp_ntcp.py` | LKB: $a=0 \Rightarrow D_{\mathrm{ref}}=D_{\max}\frac{\alpha/\beta+D_{\max}}{\alpha/\beta+2}$ (EQD2); $a>0 \Rightarrow \mathrm{EUD}=\left(\sum \Delta v_i D_i^{a}\right)^{1/a}$; $\mathrm{NTCP}=\frac12\left[1+\mathrm{erf}\!\left(t/\sqrt2\right)\right]$ en ambas ramas (autoritativo `ntcp_isoE.m`, no `ntcp.m` que omite $\sqrt2$) |
| `utcp_engine.py` | Puerto fiel de `isodose.m`: $G_{11}$ (suma, no lc), $G_{22}$, $G_{33}$, $G_{12}$, $G_{23}$, $G_{31}$, dosis lineal y cuadrática, $\mathrm{DisoE}=0.5\frac{\alpha_r/\beta_r}{G_r}\left(\sqrt{1+\frac{4G_r}{\alpha_r(\alpha_r/\beta_r)}(D_{\mathrm{lin}}+D_{\mathrm{cuad}})}-1\right)$; NaN donde no hay dosis + cobertura; `compute_utcp_point` multi-órgano y multi-mesh |
| `utcp_time_optimizer.py` | Barrido temporal (isoE **recalculada** por t, nunca escalada), curvas TCP(t)/NTCP$_j$(t)/UTCP(t), restricciones Dmax/Dmean/Dx % en Gy / Gy(RBE) / Gy(IsoE), máximo sin restricciones y máximo factible |
| `test/test_utcp_engine.py` | 30 tests unitarios sin slicer |
| `test/test_golden_matlab.py` | Validación contra MATLAB con entradas sintéticas idénticas |
| `test/matlab_ref/gen_golden.m` + `ntcp_isoE.m` | Generador de golden (copia autorizada; UTCP original intacto) + referencia |

### Modificado
- `SlicerModule/BNCT.py` (+432/−1):
  - imports UTCP (con fallback `_UTCP_IMPORT_ERROR`), `import sys`.
  - **FIX**: `_readConfigFromTable` ahora preserva claves del jsonc que la tabla no edita (`universo`, `tir`, `utcp_*`) — antes las borraba en cada edición de celda.
  - **FIX**: `_syncOrganConfigJson` conserva el top-level (`aire`) del jsonc.
  - Pestaña DVH Analysis: tabla ROI | Tipo | Modelo IsoE | Modelo NTCP | Incluir UTCP; botones *Calcular UTCP* y *Optimizar tiempo*; gráfico `Chart_UTCP` (TCP azul, UTCP rojo) y tabla `Table_UTCP`; persistencia de la selección por ROI en `organ_config.jsonc`.
  - Tumor: lee el flag existente de la tabla Dosimetry (único selector, no duplicado); 0 tumores → error; >1 → resolución explícita por **nombre** (nunca por índice).
- `AGENTS.md` (raíz): reglas del proyecto (UTCP solo lectura, tumor único manual, unidades Gy(IsoE), LaTeX, checkcode).

## 2. Funciones reutilizadas
`_getKermaVolumesForMesh` (descubrimiento de kerma por mesh), `load_jsonc`,
`_normName`, `dvhTirSpin` / `dvhB10Spin`, nodo `Mascaras_Combinadas_Gris`,
`models.jsonc` (secciones `tcp`/`ntcp` ahora sí consumidas), flag *Tumor* de la
tabla Dosimetry, `isodose.m`/`TCPc.m` como referencia de puerto (no copiados).

## 3. Modelos radiobiológicos implementados
- **TCP Provenzano** (referencia B2/B3 + 3 de `models.jsonc`).
- **NTCP LKB**: Mucositis ≥3 (TD50 = 39.8 Gy(IsoE), m = 0.17, a = 0, α/β = 10 Gy — marcado PROVISIONAL si se usa en piel), SS ≥2 Dattoli 2024 (TD50 = 13.9049 Gy(IsoE), m = 0.1766, a = 10, α/β = 3 Gy), LKB configurable.
- **RiBack 2026 (piel)**: arquitectura lista, cálculo bloqueado (`NotImplementedError`) hasta validar sus parámetros.
- $$\mathrm{UTCP} = \mathrm{TCP}_{\mathrm{tumor}} \cdot \prod_{j \in \mathrm{OAR}} \left(1 - \mathrm{NTCP}_j\right)$$ con supuesto de **independencia condicional** documentado en la salida.

## 4. Validaciones realizadas y resultados
- **34/34 tests OK**: `python -m unittest SlicerModule.test.test_golden_matlab SlicerModule.test.test_utcp_engine`.
- Golden MATLAB (mismas entradas): DisoE coincide a ~$10^{-9}$ relativo (NaN idénticos); NTCP a $10^{-9}$ (mucosa y cerebro, $t = 30/45/60$ min); TCP/PEUD a $10^{-7}$ (amplificación exponencial de $\exp(\cdot)$, documentada).
- Cubre: conversión DVH, normalización, EUD uniforme/heterogéneo, NTCP(TD50)=0.5, límites [0,1], UTCP 1 y varios órganos, sin tumor, parámetros incompletos, NaN/cobertura, distinción Gy / Gy(RBE) / Gy(IsoE), restricciones factibles.

## 5. Limitaciones y parámetros pendientes
- RiBack sin parámetros validados; LKB configurable requiere TD50/m/ab.
- Divergencias isoE preexistentes en `BNCT.py` (`calculateIsoEffectveDose` y `_calcIsoDose` difieren de `isodose.m`) documentadas y **no tocadas**; el motor UTCP usa su puerto fiel propio.
- DVH legacy de `BNCT.py` no aplica `scale` ni pesos; el motor lee kerma crudo directamente (boro × ppm × ratio).
- Dos campos comparten un único TIR (pesos independientes); restricciones dosimétricas sin UI (API lista en `utcp_time_optimizer`).
- `organ_config.jsonc` actual marca GTV 3mm como tumor — revisar el checkbox antes de calcular.

## 6. Procedimiento de prueba en 3D Slicer
1. Recargar el módulo BNCT (o reiniciar Slicer).
2. Pestaña *Dosimetry*: marcar exactamente un tumor (columna Tumor).
3. Generar máscaras → `Mascaras_Combinadas_Gris`.
4. Pestaña *DVH Analysis*: B-10 y TIR, seleccionar OARs, modelo NTCP por órgano, "Incluir en UTCP".
5. *Calcular UTCP* → valores + advertencias; *Optimizar tiempo* → curvas y $t$ óptimo (`Chart_UTCP`).

La selección persiste en `organ_config.jsonc` (claves `utcp_isoe_model`,
`utcp_ntcp_model`, `utcp_include` por ROI). El proyecto UTCP original
(`C:\Seba\BNCT\BNCT-clinica\mat\UTCP`) quedó intacto (solo lectura).
