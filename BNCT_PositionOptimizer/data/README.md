# data/ — archivos de entrada del optimizador (NO versionados, salvo este README)

Esta carpeta contiene los datos reales (pesados) que consume
`run_optimizer.m` con `modoDatos = 'real'`. Los archivos binarios/grandes
están excluidos del repositorio vía `.gitignore`; este README sí se
versiona para documentar qué debe haber acá.

## Archivos esperados

| Archivo | Tamaño aprox. | Contenido | Quién lo lee |
|---|---|---|---|
| `mesh-neutron.out` | ~150 MB | Salida PHITS `[ T-Track ]` de neutrones, KERMA en Gy/source. Contiene **3 msets** en este orden: 1 = Boron (kerma del boro), 2 = FastNeutron (neutrones rápidos), 3 = ThermalNeutron (neutrones térmicos). Grilla 100×200×200 celdas de 0.3 cm (→ 3 mm), eje X a lo largo del haz. | `loadPHITSDose` → `src/leerPHITSNeutron.m` |
| `mesh-photon.out` | ~50 MB | Salida PHITS T-Track de fotones, 1 mset = Gamma (kerma de fotones), misma grilla que neutron. | `loadPHITSDose` → `src/leerPHITSGamma.m` |
| `BODY.vtk` | ~2.6 MB | Superficie externa del PACIENTE (POLYDATA ASCII, ~40 574 puntos / 81 139 triángulos), exportada de 3D Slicer (espacio LPS, mismo marco PHITS: haz entra por x=0, fantoma x∈[0,300] mm). Se mueve con la pose T(p) y colisiona contra el puerto. | `loadReactorSTL` (parser VTK) |
| `cone_beam_con_pared_v4_1.vtk` | ~52 KB | Puerto del REACTOR (POLYDATA ASCII): cilindro r=150 mm (x∈[−400,−200]) + cono 150→85 (x∈[−200,0]) + anillo de boca r=85 en x≈0. Incluye una envolvente lejana (±2600 mm) que `run_optimizer` recorta con `clipMeshBox`. Estático (no se mueve). | `loadReactorSTL` |

Opcional: NIfTI de máscaras anatómicas reales (tumor/mucosa/piel) cuando
estén disponibles; hoy el modo real usa máscaras sintéticas de
`makeDemoData` para UTCP y geometría real solo para colisiones.

## Cómo se cargan

- Dosis PHITS: `loadPHITSDose(struct('neutronFile', 'data/mesh-neutron.out', 'photonFile', 'data/mesh-photon.out'))`.
  Los `.out` traen Gy/source (doseScale=1); la normalización clínica se
  aplica aguas abajo (ver README principal §6).
- Geometría VTK: `loadReactorSTL('data/BODY.vtk')` y
  `loadReactorSTL('data/cone_beam_con_pared_v4_1.vtk')` (el parser VTK
  ASCII maneja secciones POINTS/POLYGONS/CELLS y exponentes Fortran `D`).

## Referencia de formato (fuente de verdad)

El convenio de parseo (orden de msets, salto de cabeceras `hc:`,
dims nx×ny×nz, flips por página) es el mismo de `SlicerModule/BNCT.py`
(función `parsePHITS`, que a su vez dice "Traducido de
f_lectura_phits_neutron.m"). Componentes de neutrones esperadas:

1. **Boron** — kerma por reacción (n,α) en ¹⁰B; domina en tejido con boro.
2. **FastNeutron** — neutrones rápidos; es la componente más penetrante.
3. **ThermalNeutron** — flujo térmico; decae rápido con la profundidad.

Verificación rápida de carga correcta: los máximos de las tres
componentes deben quedar cerca de la entrada (x pequeño), el centro de
masa axial del FastNeutron debe ser más profundo que el del
ThermalNeutron, y el kerma de Boron debe estar fuertemente correlacionado
con el flujo térmico.
