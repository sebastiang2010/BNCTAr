# -*- coding: utf-8 -*-
"""utcp_ntcp — Modelos NTCP (LKB; puerto fiel de utils/ntcp_isoE.m).

Capa pura Python + numpy: NO importa Qt ni slicer.

LKB sobre DVH diferencial en Gy (IsoE):
    a > 0:  Dref = EUD = (sum dV_i * D_i^a)^(1/a)
    a == 0: Dref = EQD2(Dmax) = Dmax * (ab + Dmax) / (ab + 2)
    t = (Dref - TD50) / (m * TD50)
    NTCP = 0.5 * (1 + erf(t / sqrt(2)))        <- en AMBAS ramas

Nota: ntcp.m de la raiz MATLAB omite /sqrt(2) en la rama a=0 (bug documentado);
ntcp_isoE.m es la referencia autoritativa del proyecto UTCP y es la que se porta.

Modelo Riback 2026 (piel): arquitectura preparada, calculo bloqueado hasta
contar parametros validados (status="pending" en utcp_models).
"""

import math

import numpy as np


def validate_ntcp_params(params):
    """Valida parametros LKB; lanza ValueError si faltan o son invalidos."""
    kind = params.get("kind", "lkb")
    if kind == "riback":
        raise NotImplementedError(
            "Modelo NTCP Riback 2026 no disponible: parametros sin validar "
            "(status=%r). Configurar N0/k_EUD/k_NTCP/D_presc antes de activarlo."
            % params.get("status", "pending"))
    required = ["TD50", "m", "a"]
    missing = [k for k in required if params.get(k) is None]
    if missing:
        raise ValueError("Parametros NTCP incompletos: faltan %s" % ", ".join(missing))
    TD50, m, a = float(params["TD50"]), float(params["m"]), float(params["a"])
    if not math.isfinite(TD50) or TD50 <= 0:
        raise ValueError("NTCP: TD50 debe ser > 0 Gy(IsoE)")
    if not math.isfinite(m) or m <= 0:
        raise ValueError("NTCP: m debe ser > 0")
    if a < 0:
        raise ValueError("NTCP: a debe ser >= 0 (a=0 => EQD2 sobre Dmax)")
    if a == 0:
        if params.get("ab") is None:
            raise ValueError("NTCP: ab (alpha/beta) requerido cuando a=0")
        ab = float(params["ab"])
        if not math.isfinite(ab) or ab <= 0:
            raise ValueError("NTCP: ab debe ser > 0 Gy")
    return dict(params)


def ntcp_lkb(dose_centers, dvol, params):
    """NTCP LKB sobre DVH diferencial (Gy IsoE).

    dose_centers: centros de dosis [Gy IsoE].
    dvol: fracciones de volumen (sum = 1; se re-normaliza con aviso interno).
    params: dict {TD50, m, a, ab?, kind?}.

    Devuelve dict {"ntcp", "dref_gy", "t"}.
    """
    p = validate_ntcp_params(params)
    D = np.asarray(dose_centers, dtype=float).ravel()
    w = np.asarray(dvol, dtype=float).ravel()
    if D.size != w.size or D.size == 0:
        raise ValueError("DVH OAR vacio o desalineado")
    if not np.all(np.isfinite(D)) or not np.all(np.isfinite(w)):
        raise ValueError("DVH OAR con valores no finitos (NaN/Inf)")
    if np.any(D < 0) or np.any(w < 0):
        raise ValueError("DVH OAR con dosis o volumenes negativos")
    s = w.sum()
    if s <= 0:
        raise ValueError("DVH OAR sin volumen efectivo")
    if not np.isclose(s, 1.0, atol=1e-6):
        w = w / s

    a = float(p["a"])
    if a == 0:
        dmax = float(D.max())
        ab = float(p["ab"])
        dref = dmax * (ab + dmax) / (ab + 2.0)
    else:
        dref = float(np.power(np.sum(w * np.power(D, a)), 1.0 / a))

    t = (dref - float(p["TD50"])) / (float(p["m"]) * float(p["TD50"]))
    ntcp = 0.5 * (1.0 + math.erf(t / math.sqrt(2.0)))
    ntcp = min(max(ntcp, 0.0), 1.0)
    return {"ntcp": ntcp, "dref_gy": dref, "t": t}


def ntcp_eval(dose_centers, dvol, params):
    """Dispatcher de modelos NTCP (extension futura: riback, etc.)."""
    kind = params.get("kind", "lkb")
    if kind == "riback":
        raise NotImplementedError(
            "NTCP Riback 2026 bloqueado hasta validar sus parametros")
    if kind != "lkb":
        raise ValueError("Modelo NTCP desconocido: kind=%r" % kind)
    return ntcp_lkb(dose_centers, dvol, params)
