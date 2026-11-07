# BNCT_PositionOptimizer

Optimizador de posicionamiento del paciente para terapia por captura neutrónica de boro (BNCT).
Busca la pose rígida (3 traslaciones + 3 rotaciones) que **maximiza la UTCP** —*Utopía TCP*,
Provenzano et al. 2019— manteniendo el paciente fuera de colisión con el reactor y por debajo
de los límites clínicos de NTCP en mucosa y piel.

Todas las magnitudes del proyecto están en **mm**, **grados** y **Gy** (Gy-Eq después de aplicar
los pesos RBE). El paquete es autocontenido: `run_optimizer.m` corre de punta a punta con datos
sintéticos, sin archivos externos y sin display.

---

## 1. Propósito

El problema que resuelve el paquete es: dada una anatomía del paciente, una geometría del puerto
del reactor y cuatro mapas de dosis PHITS fijos en el marco del haz, encontrar la pose
$p=[x,y,z,\alpha,\beta,\gamma]$ que mejor cubre el tumor sin pasar los límites de toxicidad ni
chocar con la estructura.

La función de mérito es la UTCP:

$$
\mathrm{UTCP} = \mathrm{TCP} \cdot \left(1-\mathrm{NTCP}_{m}\right)^{w_m} \cdot \left(1-\mathrm{NTCP}_{p}\right)^{w_p}
$$

donde $\mathrm{NTCP}_m$ es la mucositis oral, $\mathrm{NTCP}_p$ la toxicidad de piel y $w_m, w_p$
son exponentes de penalización configurables. Si un NTCP supera su límite $N_{\max}$, la pose se
rechaza y la UTCP se fuerza a $-\infty$.

---

## 2. Estado verificado

| Aspecto | Estado |
|---|---|
| Demo sintético end-to-end | Corre completo en **MATLAB R2023b** (headless, sin display) |
| Optimización por defecto | `surrogateopt` con 150 evaluaciones ≈ **10 s** |
| Validación de finalistas | 3 finalistas con **colisión exacta** (`pointToMeshDistance`) y **robustez ±2 mm / ±2°** |
| Dependencia de datos reales | Ninguna para el demo: `makeDemoData` genera paciente, reactor y dosis |

> El demo **no** es una validación física del transporte neutrónico; es una verificación de
> software del flujo completo (carga → SDF → optimización → validación → bundle).

---

## 3. Requisitos

- **MATLAB R2023b** o superior.
- **Global Optimization Toolbox** — `surrogateopt`, `ga`, `patternsearch`, `bayesopt`.
- **Image Processing Toolbox** — `niftiread`, `niftiinfo`, `bwdist` (y `imdilate` solo en el
  relleno `shell`), `niftiwrite` para exportar máscaras.
- **Statistics and Machine Learning Toolbox** — requerido por `bayesopt`
  (`optimizableVariable`).

Funciones de terceros: ninguna. Todas las dependencias son de MathWorks.

Verificación rápida del entorno:

```matlab
ver('optim')          % Global Optimization Toolbox
ver('images')         % Image Processing Toolbox
ver('stats')          % Statistics and Machine Learning Toolbox
```

---

## 4. Estructura de archivos

```
BNCT_PositionOptimizer/
├── run_optimizer.m      % demo end-to-end (entry point)
├── README.md
├── src/                 % 15 funciones, todas públicas
└── data/                % reservado para datos del usuario (vacío en el repo)
```

### 4.1 Módulos

| Archivo | Rol |
|---|---|
| `run_optimizer.m` | Demo: arma datos sintéticos, precomputa SDF, corre la optimización, valida 3 finalistas e imprime el resumen. |
| `src/makeDemoData.m` | Genera paciente elipsoidal + tumor + OARs, reactor no convexo, dosis sintética y parámetros radiobiológicos de ejemplo. |
| `src/loadPatient.m` | Lee una máscara NIfTI, la binariza, extrae la superficie con `isosurface` y arma la grilla en mm. |
| `src/loadReactorSTL.m` | Carga la geometría del reactor desde `.stl` o `.vtk` ASCII y normaliza a `points`/`faces`. |
| `src/loadPHITSDose.m` | Carga y valida las 4 componentes de dosis PHITS (`.mat`, directorio de `.mat` o struct en memoria). |
| `src/signedDistanceMap3D.m` | Rasteriza la malla y precomputa el SDF firmado en una grilla de voxels (relleno por paridad). |
| `src/eul2rot.m` | Matriz de rotación $R = R_z R_y R_x$ a partir de ángulos de Euler en grados. |
| `src/rigidTransform.m` | Arma la homogénea 4×4 paciente → haz a partir de `p`. |
| `src/transformPatient.m` | Transforma la superficie del paciente y muestrea la dosis en tumor y OARs (Nivel A). |
| `src/sampleDoseMaps.m` | Interpola las 4 componentes de dosis en puntos ya expresados en el marco del haz. |
| `src/collisionCheck.m` | Chequeo de colisión: SDF voxelizado si está disponible, distancia exacta como respaldo. |
| `src/pointToMeshDistance.m` | Distancia exacta punto–triángulo (Ericson), sin signo; chequeo determinista de finalistas. |
| `src/calculateUTCP.m` | Calcula TCP (Perotti 2022), NTCP (LKB) y UTCP con límites clínicos. |
| `src/objectiveFunction.m` | Función de costo escalonada que combina colisión, máscaras vacías, NTCP y UTCP. |
| `src/optimizePosition.m` | Driver de optimización: `surrogateopt`, `ga`, `patternsearch`, `bayesopt` o `grid`, con historial. |
| `src/validatePHITS.m` | Valida finalistas (colisión exacta + robustez) y exporta el bundle para el recálculo PHITS. |

### 4.2 Flujo de datos

```
loadPatient ─┐
loadReactorSTL ─┼─→ data ─→ optimizePosition ─→ history ─→ seleccionarFinalistas *
loadPHITSDose ─┘                │                              │
                                │                              ▼
                                └── objectiveFunction ──→ validatePHITS ─→ bundle (Nivel B)
                                        │
                        signedDistanceMap3D (precomputado una vez)
```

\* `seleccionarFinalistas` es una función **local de `run_optimizer.m`**, no forma parte de `src/`:
toma el top-3 de poses aceptadas separadas entre sí por más de 3 mm o 1.5°, y si no hay ninguna
aceptada usa la mejor pose del historial como reserva.

---

## 5. Inicio rápido

```matlab
cd BNCT_PositionOptimizer
run_optimizer
```

`run_optimizer.m` hace todo por su cuenta: agrega `src/` al path, cambia el directorio de trabajo
al del proyecto (con `onCleanup` para restaurarlo) y ejecuta el flujo completo.

### 5.1 Qué imprime

1. Encabezado y tamaño de los datos sintéticos (voxels del paciente, triángulos del reactor).
2. Grilla y método de relleno del SDF del reactor.
3. **Línea de base**: métricas de la pose de setup desalineada `data.pInicial`
   (`UTCP`, `TCP`, `NTCP` mucosa y piel, `dmin` y motivo).
4. **Optimización**: método, número de evaluaciones, tiempo, mejor pose `p*`, matriz `T*` y
   conteo de evaluaciones por motivo (`ok`, `colision`, `ntcp_mucosa`, `ntcp_piel`, …).
5. **Validación** de los 3 finalistas: UTCP, TCP, NTCP, `dmin` exacta, rango de UTCP perturbada y
   cantidad de perturbaciones que violan límites (sobre 12).
6. Ruta del bundle para el recálculo PHITS.

Salida (formato, valores abreviados):

```
=== OPTIMIZADOR DE POSICION BNCT (demo sintetica) ===
Datos sinteticos listos: paciente 31x23x26 voxels, reactor 288 triangulos.
SDF del reactor: grilla 126x121x121 @ 2 mm (fill=parity).

--- Pose de setup inicial p0 = [10 -8 6 8 -6 4] ---
  UTCP=<...> | TCP=<...> | NTCP mucosa=<...> | NTCP piel=<...> | dmin=<...> mm | [<estado>: <motivo>]

--- Optimizacion (surrogateopt, 150 evaluaciones, <t> s) ---
Mejor pose p* = [<x> <y> <z> <alpha> <beta> <gamma>]  (mm, grados)
Matriz T* (paciente -> haz):
   <4x4>
Evaluaciones por motivo: ok=<n>, colision=<n>, ntcp_mucosa=<n>, ntcp_piel=<n>, ...
```

Los conteos y las métricas dependen de la corrida; la estructura del texto es la que se muestra.

### 5.2 Dónde queda el bundle

`makeDemoData` fija `options.bundleDir = fullfile(tempdir, 'BNCT_validacion_demo')`, así que en el
demo el bundle se escribe en:

```
tempdir/BNCT_validacion_demo/
```

Si **no** se define `options.bundleDir`, `validatePHITS` usa el valor por defecto
`tempdir/BNCT_validacion_<yyyyMMdd_HHmmss>` y devuelve la ruta como segunda salida:

```matlab
[report, bundleDir] = validatePHITS(finalistas, data);
```

---

## 6. Uso con datos reales

### 6.1 Anatomía del paciente

`loadPatient` devuelve una estructura con la máscara binaria, la superficie y la grilla:

```matlab
p = loadPatient('cuerpo.nii.gz');          % cuerpo / referencia
t = loadPatient('tumor.nii.gz');
p.tumorMask = t.mask;
p.tumorGrid = struct('gridX', t.gridX, 'gridY', t.gridY, ...
                     'gridZ', t.gridZ, 'voxelSize', t.voxelSize);

m = loadPatient('mucosa.nii.gz');
p.oars.mucosa = struct('mask', m.mask);

s = loadPatient('piel.nii.gz');
p.oars.piel = struct('mask', s.mask);
```

Campos de `patientStruct` que consume el pipeline:

| Campo | Tipo | Descripción |
|---|---|---|
| `.mask` | `[nx ny nz]` double 0/1 | Volumen de referencia (cuerpo). |
| `.tumorMask` | `[nx ny nz]` double 0/1 | **Obligatorio**. ROI tumoral. |
| `.oars.<nombre>.mask` | `[nx ny nz]` double 0/1 | OARs. El pipeline usa `mucosa` y `piel` por nombre. |
| `.surface.points` | `Nx3` mm | **Obligatorio** para colisiones; lo genera `loadPatient`. |
| `.surface.faces` | `Mx3` | Triángulos de la superficie. |
| `.voxelSize` | `[1x3]` mm | Espaciado; si falta `niftiinfo`, es `[1 1 1]`. |
| `.origin` | `[1x3]` mm | Origen de la grilla; default `[0 0 0]`. |
| `.gridX/.gridY/.gridZ` | `1xn` mm | Ejes `ndgrid` (`gridX` ↔ dimensión 1). |
| `.tumorGrid` | struct | Grilla propia del tumor (opcional). |

Cada ROI puede traer su propia grilla dentro de su entrada de `oars`; si no la trae, se usa la del
cuerpo. La convención es **`ndgrid`**: `mask(i,j,k)` con `i↔x`, `j↔y`, `k↔z`, consistente con
`sampleDoseMaps` y con los mapas PHITS.

> `loadPatient` usa `niftiread`/`niftiinfo` (Image Processing Toolbox). Si `niftiinfo` no devuelve
> metadata útil, emite un warning y cae a `voxelSize = [1 1 1] mm`, `origin = [0 0 0]`.

### 6.2 Geometría del reactor

```matlab
[wallStruct, bodyStruct] = loadReactorSTL('puerto.stl', 'cuerpo_reactor.stl');
```

Formatos soportados: `.stl` (vía `stlread`) y `.vtk` **ASCII** legacy (`POINTS` +
`POLYGONS`/`CELLS`, triángulos o quads). VTK binario y `.vtp` producen un error explícito.
`bodyPath` es opcional: si se omite, `bodyStruct` vuelve vacío y las colisiones se evalúan solo
contra `wallStruct`. Las mallas se asumen en **mm** y en el **marco del haz**.

### 6.3 Dosis PHITS

`loadPHITSDose` acepta tres formas de entrada:

```matlab
doseData = loadPHITSDose('dosis_phits.mat');     % (a) un .mat con el struct
doseData = loadPHITSDose('carpeta_dosis/');      % (b) un directorio de .mat
doseData = loadPHITSDose(miStruct);              % (c) struct ya armado
```

Formato canónico que exige (y devuelve normalizado):

| Campo | Tamaño | Unidad | Origen PHITS |
|---|---|---|---|
| `.gridX/.gridY/.gridZ` | `1xn` | mm | Ejes `ndgrid` |
| `.dose_B` | `[nx ny nz]` | Gy | boro |
| `.dose_H` | `[nx ny nz]` | Gy | neutrones rápidos |
| `.dose_N` | `[nx ny nz]` | Gy | neutrones térmicos / epitermales |
| `.dose_gamma` | `[nx ny nz]` | Gy | fotones |

Se reconocen nombres alternativos (`xc`/`ejeX` → `gridX`, `boro` → `dose_B`, `fast` → `dose_H`,
`thn`/`termico` → `dose_N`, `Dg`/`gamma` → `dose_gamma`). Los `NaN` (regiones de aire sin tally)
se reemplazan por 0 Gy y se cuentan en `.nNaN`. También se completan `.origin`, `.spacing` y
`.orientation`.

> Si en el directorio solo hay `.out` crudos de PHITS, la función **no** los parsea: pide
> convertirlos antes con los parsers del proyecto real (`f_lectura_phits_neutron.m`,
> `f_lectura_phits_gamma.m`) y guardar `.mat`.

### 6.4 SDF del reactor (precómputo obligatorio para velocidad)

```matlab
wallStruct = signedDistanceMap3D(wallStruct, ...
    struct('resolution', 2, 'margin', [120 80 80], 'fillMethod', 'parity'));
```

El resultado **reemplaza** al `wallStruct` original: conserva `points`/`faces` y agrega
`.sdf`, `.gridX`, `.gridY`, `.gridZ`, `.voxelSize`, `.occupancy`, `.fillMethod` y `.resolution`.
`collisionCheck` detecta esos campos y pasa a evaluar por interpolación trilineal (O(N)) en lugar
de distancia exacta.

Opciones:

| Opción | Default | Descripción |
|---|---|---|
| `.resolution` | `2` | Resolución isotrópica de la grilla, mm. |
| `.margin` | `50` | Margen alrededor del bbox; escalar o `1x3` en mm. **Debe cubrir toda la anatomía y su rango de movimiento.** |
| `.gridX/.gridY/.gridZ` | — | Ejes explícitos (override de `resolution`/`margin`). |
| `.fillMethod` | `'parity'` | `'parity'` (relleno sólido por paridad) o `'shell'` (pared fina dilatada 1 voxel). |

### 6.5 Parámetros radiobiológicos

`radioBio` alimenta a `calculateUTCP`:

```matlab
radioBio = struct();
radioBio.rbeTumor  = [6.0 3.0 2.0 1.0];   % [B H N gamma] Gy -> Gy-Eq
radioBio.rbeNormal = [3.0 3.0 2.0 1.0];
radioBio.tcp = struct('Ar', 0.0218, 'Br', 0.00218, 'c1', 2.435, 'c2', 0.120, ...
    'tof', 34.6320, 'tos', 1.2121e+03, ...
    'pf_g', 0.53, 'ps_g', 0.47, 'pf_bnct', 0.2, 'ps_bnct', 0.8);
radioBio.ntcp_mucosa = struct('TD50', 39.8, 'm', 0.17, 'ab', 10, 'a', 0);
radioBio.ntcp_piel   = struct('TD50', 55.0, 'm', 0.50, 'ab', 10, 'a', 1);
radioBio.w_m = 0.4;
radioBio.w_p = 0.2;
radioBio.N_max_mucosa = 0.10;
radioBio.N_max_piel   = 0.20;
```

| Campo | Descripción |
|---|---|
| `.rbeTumor`, `.rbeNormal` | Pesos `[B H N gamma]` aplicados como suma ponderada para obtener Gy-Eq. |
| `.tcp` | Struct del modelo Perotti 2022 (`TCPc.m`): `Ar, Br, c1, c2, tof, tos, pf_g, ps_g, pf_bnct, ps_bnct`. |
| `.tcpValue` | Alternativa: TCP fijo. Tiene precedencia sobre `.tcp` si existe. |
| `.ntcp_mucosa`, `.ntcp_piel` | Struct LKB (`TD50, m, ab, a`) **o** un valor numérico directo de NTCP. |
| `.w_m`, `.w_p` | Exponentes de penalización en la UTCP (default 0). |
| `.N_max_mucosa`, `.N_max_piel` | Límites clínicos 0–1 (default 1, es decir sin restricción efectiva). |

Los NTCP se calculan con dos variantes según `a`:

- `a == 0` → modelo de **dosis máxima** (Gonzalez 2017), con
  $\mathrm{EQD2} = D\,(\mathit{ab}+D)/(\mathit{ab}+2)$.
- `a ~= 0` → modelo **EUD** (Dattoli Viegas 2024), con
  $\mathrm{EUD} = \left(\sum_i D_i^{a} v_i\right)^{1/a}$.

Los valores del ejemplo son los del demo; `ntcp_piel` y los pesos RBE son *placeholders*
didácticos.

### 6.6 Ensamblado del struct `data`

```matlab
data = struct();
data.patientStruct  = p;              % anatomía (6.1)
data.wallStruct     = wallStruct;     % reactor con SDF (6.2 + 6.4)
data.bodyStruct     = bodyStruct;     % opcional
data.doseData       = doseData;       % 4 componentes PHITS (6.3)
data.radioBio       = radioBio;       % radiobiología (6.5)
data.safetyDistance = 1;              % mm (default 1)

data.options = struct();
data.options.robustDx     = 2;                          % mm (validatePHITS)
data.options.robustDtheta = 2;                          % grados (validatePHITS)
data.options.bundleDir    = fullfile(pwd, 'bundle_nivelB');
data.options.log          = containers.Map;             % registro de evaluaciones
data.options.onEvaluation = @(pp, c, info) fprintf('%s %.2f\n', info.reason, c);
```

Campos de `data.options` reconocidos por el código:

| Campo | Consumido por | Efecto |
|---|---|---|
| `.log` (`containers.Map`) | `objectiveFunction` | Guarda `eval_000001`, `eval_000002`, … con `p`, `cost`, `reason`, `dmin`. |
| `.onEvaluation` (`function_handle`) | `objectiveFunction` | Callback `@(p, cost, info)` en cada evaluación. |
| `.robustDx` | `validatePHITS` | Perturbación de traslación, mm (default 2). |
| `.robustDtheta` | `validatePHITS` | Perturbación de rotación, grados (default 2). |
| `.bundleDir` | `validatePHITS` | Carpeta del bundle (default `tempdir/BNCT_validacion_<timestamp>`). |

> `data.pInicial` solo lo usa el demo como pose de setup; el optimizador no lo necesita.

---

## 7. Algoritmos disponibles

```matlab
[result, history] = optimizePosition(data, opts);
```

### 7.1 Opciones de `optimizePosition`

| Opción | Default | Descripción |
|---|---|---|
| `.method` | `'surrogateopt'` | `'surrogateopt'` \| `'ga'` \| `'patternsearch'` \| `'bayesopt'` \| `'grid'`. |
| `.lb` | `[-30 -30 -30 -15 -15 -15]` | Límites inferiores (mm, grados). |
| `.ub` | `[30 30 30 15 15 15]` | Límites superiores (mm, grados). |
| `.maxFunEvals` | `500` | Presupuesto de evaluaciones. |
| `.x0` | centro de `[lb, ub]` | Semilla para `patternsearch` y `grid`. |
| `.verbose` | `false` | `true` → `Display = 'iter'`; `false` → `'off'` (headless-safe). |
| `.plot` | `false` | Gráfico de convergencia si hay display (si falla, solo emite warning). |

Se valida `lb < ub` componente a componente.

### 7.2 Comparación

| Método | Mapeo de `maxFunEvals` | Notas |
|---|---|---|
| `surrogateopt` | `MaxFunctionEvaluations` | **Default**. Adecuado para objetivos con mesetas y rechazos (el costo es escalonado). |
| `ga` | `MaxGenerations = max(5, ceil(maxFunEvals/50))`, `PopulationSize = 50` | Población fija de 50; el presupuesto real es `población × generaciones`. |
| `patternsearch` | `MaxFunctionEvaluations` | Arranca en `x0`; bueno para refinamiento local. |
| `bayesopt` | `MaxObjectiveEvaluations` | Variables `tx, ty, tz, alpha, beta, gamma`; objetivo determinista. Modela colisión y NTCP como **restricciones** (`bayesWrap`). |
| `grid` | No usa presupuesto | Barrido coordenado determinista: 2 pasadas, 5 puntos por dimensión (la segunda pasada sobre una ventana reducida). |

Los métodos basados en toolboxes verifican la existencia de la función antes de correr y fallan con
un mensaje claro si falta la toolbox; `grid` no requiere ninguna.

> `bayesopt` accede a `data.safetyDistance` y a `data.radioBio.N_max_mucosa` / `.N_max_piel` para
> armar las restricciones: esos campos deben existir cuando se elige este método.

### 7.3 Ejemplos

```matlab
% Default: surrogateopt
[result, history] = optimizePosition(data, struct('method', 'surrogateopt', 'maxFunEvals', 150));

% Búsqueda global con población
[result, history] = optimizePosition(data, struct('method', 'ga', 'maxFunEvals', 300));

% Refinamiento local desde una pose conocida
[result, history] = optimizePosition(data, ...
    struct('method', 'patternsearch', 'x0', result.pBest, 'maxFunEvals', 80));

% Barrido determinista (sin toolbox de optimización global)
[result, history] = optimizePosition(data, struct('method', 'grid'));
```

### 7.4 Salidas

`result`:

| Campo | Descripción |
|---|---|
| `.pBest` | `[1x6]` mejor pose (mm, grados). |
| `.TBest` | `4x4` de `rigidTransform(pBest)`. |
| `.utcp`, `.tcp`, `.ntcp`, `.dmin`, `.cost`, `.reason` | Métricas recalculadas en `pBest`. |
| `.counts` | Conteo por motivo: `ok`, `colision`, `ntcp_mucosa`, `ntcp_piel`, `tumor_vacio`, `oar_vacio`, `otros`. |
| `.nEvals` | Evaluaciones registradas por el wrapper del solver (no cuenta el recálculo final de métricas en `pBest`). |
| `.method`, `.elapsed`, `.lb`, `.ub` | Método usado, tiempo y límites efectivos. |

`history`: `struct` con `.p` (`Nx6`), `.cost` (`Nx1`), `.reason` (`Nx1` cell) y `.dmin` (`Nx1`),
una fila por evaluación.

### 7.5 Función de costo

`objectiveFunction(p, data)` devuelve un escalar **finito** siempre (nunca `-Inf` ni `NaN`), para
que los solvers no se rompan:

| Situación | Costo | `info.reason` |
|---|---|---|
| `dmin < safetyDistance` | `1e6 + 1e3*(safety - dmin)` | `'colision'` |
| Máscara tumoral u OAR vacía | `1e5` | `'tumor_vacio'` / `'oar_vacio'` |
| Límite NTCP violado | `1e5 + 1e3*exceso` | `'ntcp_mucosa'` / `'ntcp_piel'` |
| Caso aceptable | `-UTCP` | `'ok'` |

El orden refleja la gravedad: ninguna pose en colisión puede ganarle a una pose válida, y ninguna
violación de NTCP puede ganarle a una aceptable.

---

## 8. Estrategia de dos niveles

### Nivel A — ranking sobre dosis fijas (lo que hace este paquete)

Las 4 matrices PHITS quedan **fijas en el marco del haz**. Para cada pose candidata `T(p)` se
transforman las máscaras del paciente al marco del haz y se muestrea la dosis en cada voxel:

$$
D_k^{p}(r) = D_k^{\text{PHITS}}\!\left(T(p)\,r\right), \qquad
r_{\text{haz}} = R\,r_{\text{pac}} + t
$$

Con esas muestras se calculan TCP, NTCP y UTCP, y las poses se clasifican. Este nivel es barato
(una interpolación por voxel y por componente) y permite explorar cientos de poses.

```matlab
[patientT, sampledDose] = transformPatient(data.patientStruct, T, data.doseData);
[tcp, ntcp, utcp, info] = calculateUTCP(sampledDose, patientT, data.radioBio);
```

`sampledDose` devuelve `.tumor` (`Nt×4`, `[D_B D_H D_N D_gamma]`), `.oars.<nombre>.dose`,
`.oars.<nombre>.mean`, los volúmenes de voxel en cm³ y `.voxelsFueraDominio` por región.

### Nivel B — recálculo PHITS de los finalistas

`validatePHITS` **no ejecuta PHITS**. Deja en el bundle todo lo necesario (pose, `T`, anatomía
transformada en NIfTI) para que el usuario lance el transporte real con la geometría ya
posicionada. Después se comparan las métricas del Nivel B contra el ranking del Nivel A.

### ⚠️ Advertencia sobre el Nivel A

**La dosis del cilindro trasladada no reproduce el transporte real.** Mover las máscaras sobre un
mapa de dosis fijo ignora, entre otras cosas:

- el cambio real de moderación y atenuación al modificar la cantidad de tejido en el camino;
- la absorción y el scattering que dependen de la geometría efectiva;
- la distribución de boro (que depende del tejido que ocupa cada voxel, no de la posición del
  paciente);
- los contornos y la fuga por las paredes del puerto.

Por eso el Nivel A es **solo un ranking de candidatas**. Ninguna conclusión clínica debe sacarse
sin el recálculo del Nivel B.

---

## 9. Colisiones y SDF

### 9.1 Dos niveles de chequeo

| Contexto | Método | Costo | Precisión |
|---|---|---|---|
| Bucle del optimizador | SDF voxelizado (`signedDistanceMap3D` + interpolación trilineal) | O(N) por pose | ~1 voxel (`resolution` mm) |
| Validación de finalistas | `pointToMeshDistance` (punto–triángulo exacto) | O(N×M) con prefiltro de esfera envolvente | Exacta, sin signo |

`collisionCheck` elige automáticamente:

1. Si la malla trae `.sdf` + `.gridX/.gridY/.gridZ` → interpola el SDF (rápido).
2. Si no hay SDF y la malla tiene menos de 5000 triángulos → distancia exacta.
3. Si no hay SDF y la malla es grande → usa igual la distancia exacta y avisa **una vez por
   sesión** que conviene precomputar el SDF.

```matlab
[collision, dmin] = collisionCheck(patientStruct, wallStruct, bodyStruct, T, safetyDistance);
% collision = dmin < safetyDistance
```

`dmin` es la distancia mínima del paciente al reactor: **positiva fuera del sólido, negativa
adentro** cuando se usa el SDF. Con la verificación exacta es sin signo. `safetyDistance` es
opcional y vale 1 mm por defecto.

### 9.2 Por qué no `alphaShape` ni convex hull

El puerto del reactor tiene una **abertura cóncava** (el orificio por donde entra el haz). Un
`alphaShape` o un *convex hull* cerrarían esa concavidad y convertirían el hueco en sólido, con lo
cual el paciente chocaría siempre contra el reactor. El paquete reconstruye el sólido **por paridad
de cruces** (`fillMethod = 'parity'`, ray casting a lo largo de z columna por columna), que
preserva los túneles perpendiculares al eje de relleno — como la abertura del haz.

- `'parity'` requiere una malla cerrada (2-manifold).
- Túneles **paralelos** al eje de relleno (z) se cerrarían: en ese caso usar
  `fillMethod = 'shell'` o reorientar el relleno.

El SDF se calcula como:

$$
\text{SDF} = \Delta \cdot \left(\text{bwdist}(occ) - \text{bwdist}(\neg occ)\right)
$$

con $\Delta$ = `resolution` en mm (grilla isotrópica).

### 9.3 Fuera de dominio

Los puntos que caen fuera de la caja del SDF se evalúan como la **distancia a la caja** de la
grilla:

```matlab
gap = max(max(lo - P, P - hi), [], 2);
gap = max(gap, 0);
d(out) = vecnorm(gap, 2, 2);
```

Es una **cota inferior conservadora**: nunca subestima el riesgo. Puede generar colisiones
espurias solo en poses extremas, por eso el `.margin` del SDF debe cubrir toda la anatomía y su
rango de movimiento.

---

## 10. Validación y bundle (`validatePHITS`)

```matlab
[report, bundleDir] = validatePHITS(finalistas, data);
```

`finalistas` acepta `K×6`, un cell array de `1×6` o un struct array con campo `.p`.

### 10.1 Qué hace por cada finalista

1. **Colisión exacta**: `pointToMeshDistance` de la superficie transformada contra `wallStruct` y
   `bodyStruct` (independiente del SDF).
2. **Métricas**: `TCP`, `NTCP` (mucosa y piel) y `UTCP` por el mismo camino que el optimizador.
3. **Robustez**: 12 perturbaciones axiales — ±`robustDx` (2 mm) en cada traslación y
   ±`robustDtheta` (2°) en cada rotación. Se reporta el spread de UTCP y cuántas perturbaciones
   violan límites (NTCP o colisión por SDF). La pose nominal se reporta aparte.
4. **Exportación** de las máscaras transformadas y de la pose.

### 10.2 Archivos generados

| Archivo | Contenido |
|---|---|
| `transform.json` | Por finalista: `p`, `T` (4×4), `tcp`, `ntcpMucosa`, `ntcpPiel`, `utcp`, `dminSDF_mm`, `dminExacta_mm`, `colisionExacta`, `utcpPertMin/Max/Media`, `nPertViolan`. Encabezado con `unidadTraslacion: "mm"`, `unidadRotacion: "grados"` y la convención. |
| `resumen.csv` | Una fila por finalista con `id, x_mm, y_mm, z_mm, alpha_deg, beta_deg, gamma_deg, tcp, ntcp_mucosa, ntcp_piel, utcp, dmin_exacta_mm, utcp_pert_min, utcp_pert_max, n_pert_violan`. |
| `finalistaK_tumor.nii` | Máscara tumoral transformada (`niftiwrite`). |
| `finalistaK_<oar>.nii` | Máscara de cada OAR transformada (una por OAR presente). |
| `finalistaK_pose.txt` | `p` y `T` en texto plano, listos para copiar a otro script. |

Si `niftiwrite` no está disponible, cada máscara se guarda como `finalistaK_<region>.mat` con la
variable `maskT`.

### 10.3 Cómo se construyen las máscaras exportadas

Se hace un *warp* inverso: para cada punto `q` de la grilla original se consulta la máscara
original en `T⁻¹·q` con interpolación **nearest**, y se umbraliza en 0.5. El volumen resultante,
leído con la misma grilla que la entrada, representa la anatomía ya desplazada por `T`
(marco del haz).

> **Cuidado:** `niftiwrite` se llama solo con la matriz, sin `voxelSize` ni `origin`, así que el
> header del NIfTI **no hereda** el espaciado ni el origen originales. Al cargar el archivo hay
> que reasignar el espaciado y el origen de la grilla de entrada para que las coordenadas sean
> correctas.

### 10.4 Uso de `T` en 3D Slicer / PHITS

`T` es la homogénea 4×4 que mapea **paciente → haz**; su bloque `T(1:3,1:3)` es la rotación y
`T(1:3,4)` la traslación. En `pose.txt` y `transform.json` se escribe completa, con la aclaración
de la convención.

- En **PHITS**, usar `T` para ubicar la geometría del paciente respecto del eje del haz: la
  geometría de entrada está en el marco del paciente y `T` la lleva al marco donde viven los
  tallies.
- En **3D Slicer**, `T` se carga como transformación lineal 4×4. Slicer trabaja internamente en
  RAS y en archivos en LPS: hay que verificar la convención de ejes antes de confiar en la
  superposición, porque una conversión omitida invierte dos ejes y da una pose visualmente
  plausible pero incorrecta.

### 10.5 Campos de `report`

`report` es un struct array con un elemento por finalista:
`id, p, T, tcp, ntcp, utcp, dminSDF, dminExacta, colisionExacta, utcpPertMin, utcpPertMax,
utcpPertMedia, utcpPertStd, nPertViolan, archivos` y, al final, `bundleDir`.

---

## 11. Convención de transformación

Toda la geometría se mueve con una única convención, documentada en `rigidTransform.m`,
`eul2rot.m` y en el JSON del bundle:

$$
r_{\text{haz}} = R(\alpha,\beta,\gamma)\, r_{\text{pac}} + t
$$

$$
R = R_z(\gamma)\, R_y(\beta)\, R_x(\alpha)
$$

con `t = [x y z]` en mm y `[alpha beta gamma]` en **grados** (se convierten a radianes dentro de
`eul2rot`). `p = 0` es la **pose de referencia PHITS**: el paciente tal como está definido en su
marco original, sin desplazamiento ni rotación.

```matlab
T = rigidTransform([10 -8 6 8 -6 4]);   % paciente -> haz
rh = T(1:3,1:3) * rp + T(1:3,4);        % un punto
```

> La combinación `Rz*Ry*Rx` es intencional: se define una sola vez en `eul2rot.m` y `rigidTransform.m`,
> y todos los demás módulos consumen `T` ya armada. No cambiarla sin actualizar también este README:
> una inversión de orden produce poses que parecen razonables pero apuntan a otro lado.

---

## 12. Referencia de la API pública

### `patientStruct = loadPatient(maskPath)`

| Argumento | Tipo | Descripción |
|---|---|---|
| `maskPath` | char/string | Ruta a un NIfTI (`.nii` / `.nii.gz`). |

Devuelve `.rawData`, `.mask` (double 0/1), `.surface` (`.points` `Nx3` mm, `.faces` `Mx3`),
`.voxelSize`, `.origin`, `.gridX/.gridY/.gridZ`.

```matlab
p = loadPatient('cuerpo.nii.gz');
```

### `[wallStruct, bodyStruct] = loadReactorSTL(wallPath, bodyPath)`

`bodyPath` es opcional. Devuelve `points`, `faces` (1-based) y `fuente`.

```matlab
[wall, body] = loadReactorSTL('puerto.stl');
```

### `doseData = loadPHITSDose(source)`

`source`: ruta a `.mat`, ruta a directorio o struct en memoria. Ver §6.3.

### `sdfStruct = signedDistanceMap3D(meshStruct, options)`

`options`: `.resolution` (2), `.margin` (50), `.gridX/.gridY/.gridZ`, `.fillMethod` (`'parity'`).
Devuelve `.sdf` (mm, positivo afuera), `.gridX/.gridY/.gridZ`, `.voxelSize`, `.points`, `.faces`,
`.occupancy`, `.fillMethod`, `.resolution`.

```matlab
wall = signedDistanceMap3D(wall, struct('resolution', 2, 'margin', [120 80 80]));
```

### `[collision, dmin] = collisionCheck(patientStruct, wallStruct, bodyStruct, T, safetyDistance)`

`safetyDistance` opcional (default 1 mm). `dmin` con signo si hay SDF.

### `[dmin, idxTri] = pointToMeshDistance(points, meshPoints, meshFaces)`

Distancia **sin signo** en mm de cada punto (`Nx3`) a la malla; devuelve también el índice del
triángulo más cercano. Prefiltro por esfera envolvente por triángulo.

```matlab
[d, tri] = pointToMeshDistance(ptsHaz, wall.points, wall.faces);
dmin = min(d);
```

### `R = eul2rot(angles)` · `T = rigidTransform(p)`

`angles = [alpha beta gamma]` en grados → `R` ortonormal 3×3. `p = [x y z alpha beta gamma]` →
`T` 4×4.

### `[patientT, sampledDose] = transformPatient(patientStruct, T, doseData)`

Transforma la superficie, muestrea las 4 componentes en tumor y OARs y adjunta `.T`. Los voxels
fuera del dominio de dosis reciben 0 Gy y se cuentan en `.voxelsFueraDominio`.

### `dose = sampleDoseMaps(doseData, points)`

`points` ya en el marco del haz (`Nx3`, mm). Devuelve `.B`, `.H`, `.N`, `.gamma` (`Nx1`, Gy) y
`.nFueraDominio`. Usa `interpn` con `extrapval = 0` y formato `ndgrid` (usar `meshgrid` invertiría
los ejes).

### `[tcp, ntcp, utcp, info] = calculateUTCP(dose, patientStruct, params)`

Ver §6.5. `info` trae `.ok`, `.reason`, `.limitsViolated`, `.deqTumorMean`, `.tumorVoxels`.

```matlab
[tcp, ntcp, utcp, info] = calculateUTCP(sampledDose, patientT, data.radioBio);
if ~info.ok
    fprintf('Rechazada: %s\n', info.reason);
end
```

### `[cost, info] = objectiveFunction(p, data)`

Costo escalonado (ver §7.5). `info` incluye `.p`, `.T`, `.dmin`, `.collision`, `.tcp`, `.ntcp`,
`.utcp`, `.reason`, `.limitsViolated`, `.detail`.

### `[result, history] = optimizePosition(data, opts)`

Ver §7.

### `[report, bundleDir] = validatePHITS(finalistas, data)`

Ver §10.

### `data = makeDemoData()`

Genera el set sintético completo: paciente elipsoidal (semiejes 30×22×25 mm, voxel 2 mm), tumor
esférico de radio 9 mm en `(8,9,3)`, OARs en capas (mucosa `ρ ∈ (0.80, 0.92]`, piel
`ρ ∈ (0.92, 1.00]`), reactor tubular no convexo (x ∈ [40, 90] mm, radio externo 60 mm, abertura de
radio 25 mm), dosis sintética con atenuación exponencial y perfil gaussiano, radiobiología de
ejemplo, `safetyDistance = 1` y `pInicial = [10, -8, 6, 8, -6, 4]`.

---

## 13. Limitaciones conocidas

1. **Precisión del SDF ≈ 1 voxel.** Con `resolution = 2` mm, el SDF puede errar en el orden de un
   voxel. Sirve para cribar; la decisión final se toma con `pointToMeshDistance`.
2. **Máscaras exportadas con `interpn` nearest.** La discretización introduce escalones de un
   voxel en los bordes de las ROIs del bundle.
3. **Dosis estática en el marco del haz.** El Nivel A no re-resuelve el transporte: mover al
   paciente no cambia la física (ver la advertencia de §8).
4. **No se modela la camilla** ni otros accesorios: si están en el camino, hay que agregarlos como
   `bodyStruct` o como malla adicional.
5. **Robustez solo axial.** `validatePHITS` perturba un eje por vez (±dx, ±dθ); no evalúa
   perturbaciones combinadas ni correlacionadas, ni errores de registro.
6. **El demo sintético no es una validación física.** Verifica el software, no el transporte
   neutrónico ni la dosimetría.
7. **Fuera de dominio de dosis = 0 Gy.** Puede subestimar la dosis en bordes si la grilla PHITS no
   cubre toda la anatomía en todo el rango de poses.
8. **Composición de Deq simplificada.** `calculateUTCP` usa una suma ponderada simple
   `Deq = [D_B D_H D_N D_gamma] · rbe`. El modelo completo del proyecto (`ntcp2.m`) combina las
   componentes con factores de reparación `G` y términos `A, B`. La sección marcada
   `COMPOSICION DEQ (INTEGRACION)` está pensada para reemplazarse por esa conversión o para llamar
   directamente a `ntcp2.m` / `TCPc.m`.
9. **TCP genérico de reserva.** Si `radioBio` no trae `.tcp` ni `.tcpValue`, se usa un logístico
   sobre la dosis media y se emite un warning.
10. **Valores del demo son de ejemplo.** `makeDemoData` marca explícitamente `ntcp_piel` y los pesos
    RBE como *placeholders* didácticos; `w_m`, `w_p` y los `N_max` también son valores de demo, no
    criterios clínicos.
11. **`interpn` para dosis es lineal; para máscaras es nearest.** No confundir ambos caminos:
    dosis se interpola (suave), etiquetas se remuestrean con vecino (discreto).
12. **`grid` no es un optimizador global.** Es un barrido coordenado de 2 pasadas: útil como
    línea base reproducible, no como búsqueda exhaustiva.

---

## 14. TODO abierto — definir antes del Nivel B

> **Pregunta pendiente:** ¿las 4 matrices PHITS (`dose_B`, `dose_H`, `dose_N`, `dose_gamma`)
> corresponden al campo del **cilindro completo** (geometría del puerto con el fantoma de
> referencia) o **ya incluyen el boro específico de la esfera tumoral**?
>
> La respuesta cambia por completo la interpretación del Nivel A:
>
> - Si son del cilindro completo, la componente de boro del tumor está **subrepresentada** y el
>   ranking favorece poses que no reflejan la captura real en la lesión.
> - Si ya incluyen el boro tumoral, la dosis de boro está **atada a la posición original del
>   tumor** y trasladar la máscara la vuelve inconsistente.
>
> Definir esto **antes** de lanzar el recálculo del Nivel B, y dejar registrada la decisión junto
> con los archivos de dosis usados.

---

## 15. Referencias

- Provenzano et al. (2019) — formulación de la UTCP (*Utopía TCP*).
- Perotti Bernardini et al. (2022), *Radiat. Res.* 198 — modelo TCP implementado en `TCPc.m`
  (port en `calculateUTCP`).
- Gonzalez et al. (2017), *Phys. Med. Biol.* 62, 7938 — NTCP LKB con `a = 0` (dosis máxima),
  port de `ntcp.m`.
- Dattoli Viegas et al. (2024), *Phys. Med.* 128, 104840 — NTCP con `a ≠ 0` (EUD).
- Ericson, *Real-Time Collision Detection* — algoritmo de punto más cercano en triángulo usado en
  `pointToMeshDistance`.

---

## 16. Notas de mantenimiento

- El paquete no tiene dependencias fuera de MATLAB y sus toolboxes: se puede copiar la carpeta
  `BNCT_PositionOptimizer/` a cualquier instalación R2023b+ y correr `run_optimizer`.
- `src/` contiene solo funciones públicas; no hay helpers privados fuera de los archivos listados.
- Si se cambia la convención de rotación, hay que tocar `eul2rot.m` y `rigidTransform.m` (donde vive
  la definición) y actualizar la sección §11 de este README. El resto de los módulos consume `T` ya
  armada y no necesita cambios.
- Para depurar una pose puntual sin correr la optimización:

```matlab
[cost, info] = objectiveFunction([0 0 0 0 0 0], data);
fprintf('UTCP=%.4f motivo=%s dmin=%.2f mm\n', info.utcp, info.reason, info.dmin);
```
