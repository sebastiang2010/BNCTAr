# -*- coding: utf-8 -*-
"""utcp_dvh — Conversion y validacion de DVH para TCP/NTCP/UTCP.

Capa pura Python + numpy: NO importa Qt ni slicer (testeable fuera de 3D Slicer).

Unidades:
- dosis en Gy (IsoE) para todo calculo radiobiologico.
- volumenes fraccionarios normalizados (sum = 1); volumen total en cm^3.

Puerto fiel de la logica diferencial de utils/ntcp_isoE.m (UTCP) y del
adaptador acumulado->diferencial requerido por f_dvh.m (el DVH acumulado
de f_dvh nunca llega a 0 y tiene dosis repetidas -> dedupe + normalizacion).
"""

import numpy as np


def voxel_volume_cc(spacing_mm):
    """Volumen de voxel en cm^3 a partir de spacing (mm).

    spacing_mm: secuencia (sx, sy, sz) en milimetros.
    """
    sx, sy, sz = [float(s) for s in spacing_mm]
    if sx <= 0 or sy <= 0 or sz <= 0:
        raise ValueError("spacing debe ser positivo (mm), got %r" % (spacing_mm,))
    return (sx * sy * sz) / 1000.0


def dvh_from_voxels(dose_iso_voxels, total_volume_cc=None):
    """Construye un DVH diferencial desde la muestra 1D de voxels de una ROI.

    dose_iso_voxels: array 1D de dosis en Gy (IsoE) por voxel (NaNs permitidos;
        se descartan los NaN con aviso en el resultado).
    total_volume_cc: volumen de la ROI en cm^3 (por defecto N_valid * voxel vol
        no se conoce aqui; si es None se usa 1.0 y el caller corrige con vol_cc).

    Devuelve dict: {"dose_centers", "dvol" (fraccion, sum=1), "n_voxels_valid",
    "n_voxels_nan", "coverage" (validos/total), "volume_cc"}.
    """
    d = np.asarray(dose_iso_voxels, dtype=float).ravel()
    finite = np.isfinite(d)
    n_nan = int((~finite).sum())
    dv = d[finite]
    if dv.size == 0:
        raise ValueError("ROI sin voxels con dosis finita (Gy IsoE)")
    vals, counts = np.unique(dv, return_counts=True)
    dvol = counts.astype(float) / counts.sum()
    vol = float(total_volume_cc) if total_volume_cc else 1.0
    return {
        "dose_centers": vals,
        "dvol": dvol,
        "n_voxels_valid": int(dv.size),
        "n_voxels_nan": n_nan,
        "coverage": float(dv.size) / float(d.size),
        "volume_cc": vol,
    }


def cumulative_to_differential(dvh_dose_asc, dvh_vol_pct):
    """Convierte DVH acumulado (dose asc, volumen % acumulado) a diferencial.

    Puerto de cumulative2differential (B2_v30/B3_v15) + ntcp_isoE.m con las
    correcciones necesarias para el DVH de f_dvh.m:
    1. dedupe de dosis repetidas (interp1 de MATLAB falla con repetidos).
    2. dV = -diff(V/100); el ultimo bin hereda el decrecimiento final;
       clamp >= 0 (DVH monotono decreciente puede tener mesetas).
    3. normalizacion sum(dV) = 1; si la suma es 0 -> distribucion uniforme
       sobre los bins (caso DVH plano, comportamiento de ntcp_isoE.m).
    """
    D = np.asarray(dvh_dose_asc, dtype=float).ravel()
    Y = np.asarray(dvh_vol_pct, dtype=float).ravel()
    if D.size != Y.size:
        raise ValueError("dvh_dose y dvh_vol deben tener igual longitud")
    if D.size < 2:
        raise ValueError("DVH acumulado requiere >= 2 puntos")
    order = np.argsort(D)
    D = D[order]
    Y = Y[order]
    # dedupe: promedia Y entre dosis identicas
    uD, inv = np.unique(D, return_inverse=True)
    uY = np.zeros(uD.size)
    cnt = np.zeros(uD.size)
    np.add.at(uY, inv, Y)
    np.add.at(cnt, inv, 1.0)
    uY = uY / np.maximum(cnt, 1.0)
    # diferencial (fiel a ntcp_isoE.m): dV = -diff(V), append 0, clamp, normalizar
    V = uY / 100.0
    if np.any(np.diff(uD) <= 0):
        raise ValueError("dosis repetidas o no crecientes tras dedupe")
    dV = np.append(-np.diff(V), 0.0)
    dV = np.clip(dV, 0.0, None)
    total = dV.sum()
    if total <= 0:
        dV = np.full(dV.size, 1.0 / dV.size)
    else:
        dV = dV / total
    return uD, dV


def validate_differential(dose_centers, dvol):
    """Valida un DVH diferencial: finito, no negativo, normalizado."""
    D = np.asarray(dose_centers, dtype=float)
    w = np.asarray(dvol, dtype=float)
    if D.size != w.size:
        raise ValueError("dose_centers y dvol de distinta longitud")
    if not np.all(np.isfinite(D)) or not np.all(np.isfinite(w)):
        raise ValueError("DVH diferencial con valores no finitos")
    if np.any(w < 0):
        raise ValueError("dvol con fracciones negativas")
    if np.any(D < 0):
        raise ValueError("dosis negativa en DVH")
    s = w.sum()
    if not np.isclose(s, 1.0, atol=1e-6):
        raise ValueError("dvol no normalizado (sum=%g)" % s)
    return True


def eud(dose_centers, dvol, a):
    """EUD = (sum dV_i * D_i^a)^(1/a) para a != 0 (Gy)."""
    D = np.asarray(dose_centers, dtype=float)
    w = np.asarray(dvol, dtype=float)
    a = float(a)
    if a == 0:
        raise ValueError("EUD indefinido para a=0; usar dosis de referencia EQD2")
    if np.any(D < 0):
        raise ValueError("EUD con dosis negativa")
    s = float(np.sum(w * np.power(D, a)))
    if s <= 0:
        return 0.0
    return float(np.power(s, 1.0 / a))


def dose_at_volume(dvh_dose_asc, dvh_vol_pct, pct):
    """Dx%: dosis que cubre pct % del volumen (interp lineal sobre acumulado)."""
    D = np.asarray(dvh_dose_asc, dtype=float).ravel()
    Y = np.asarray(dvh_vol_pct, dtype=float).ravel()
    if D.size == 0:
        raise ValueError("DVH vacio")
    return float(np.interp(float(pct), Y[::-1], D[::-1]))


def volume_at_dose(dvh_dose_asc, dvh_vol_pct, dose):
    """VDy: % de volumen que recibe >= dose Gy (interp lineal sobre acumulado)."""
    D = np.asarray(dvh_dose_asc, dtype=float).ravel()
    Y = np.asarray(dvh_vol_pct, dtype=float).ravel()
    if D.size == 0:
        raise ValueError("DVH vacio")
    return float(np.interp(float(dose), D, Y))
