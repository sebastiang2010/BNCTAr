# -*- coding: utf-8 -*-
"""utcp_tcp — Modelo TCP de Provenzano (puerto de TCPc.m y B2_v30/B3_v15).

Capa pura Python + numpy: NO importa Qt ni slicer.

Formulacion (TCPc.m):
    H_i = exp(-Ar*D_i - Gr*Br*D_i^2)^(nf/c2)
    TCP = exp(-c1 * (sum_i H_i * dv_i)^c2)
    dv_i = fraccion de volumen * V_tumor [cm^3]

Gr (factor de reparacion sublethal, doble exonencion de Catcheside, t=30 min,
identico a isodose.m):
    xf = tof/t, xs = tos/t;  Gf = 2xf(1-xf(1-e^{-1/xf})); Gs = 2xs(1-xs(1-e^{-1/xs}))
    Gr = Gs - pf_g*(Gs - Gf)      (f_g=1, f_bnct=0 en isodose.m)

Dosis en Gy (IsoE). TCP en [0, 1].
"""

import math

import numpy as np


def catcheside_g(tof, tos, tir, pf_g):
    """Gr para dosis continua de t minutos (isodose.m con f_g=1, f_bnct=0)."""
    t = float(tir)
    if t <= 0:
        raise ValueError("tir debe ser > 0 minutos")
    xf = float(tof) / t
    xs = float(tos) / t
    Gf = 2.0 * xf * (1.0 - xf * (1.0 - math.exp(-1.0 / xf)))
    Gs = 2.0 * xs * (1.0 - xs * (1.0 - math.exp(-1.0 / xs)))
    pf = float(pf_g)
    return Gs - pf * (Gs - Gf)


def validate_tcp_params(params):
    """Valida parametros TCP Provenzano; lanza ValueError si faltan/son invalidos."""
    required = ["c1", "c2", "Ar", "Br", "tof", "tos", "pf_g"]
    missing = [k for k in required if params.get(k) is None]
    if missing:
        raise ValueError("Parametros TCP incompletos: faltan %s" % ", ".join(missing))
    if params.get("fracciones") is None:
        params = dict(params, fracciones=1)
    for k in ["c1", "c2", "Ar", "Br", "fracciones"]:
        v = float(params[k])
        if not math.isfinite(v) or v <= 0:
            raise ValueError("TCP: %s debe ser > 0 finito (got %r)" % (k, v))
    if not (0.0 <= float(params["pf_g"]) <= 1.0):
        raise ValueError("TCP: pf_g debe estar en [0, 1]")
    if float(params["tof"]) <= 0 or float(params["tos"]) <= 0:
        raise ValueError("TCP: tof/tos deben ser > 0")
    return dict(params)


def tcp_provenzano(dose_gy_isoe, dvol, tumor_volume_cc, params, tir=30.0):
    """TCP de Provenzano sobre DVH diferencial.

    dose_gy_isoe: centros de dosis [Gy IsoE].
    dvol: fracciones de volumen normalizadas (sum = 1).
    tumor_volume_cc: volumen tumoral [cm^3].
    params: dict {c1, c2, Ar, Br, fracciones, tof, tos, pf_g}.
    tir: duracion de irradiacion [min] (por defecto 30, como TCPc.m).

    Devuelve dict {"tcp", "peud", "g_factor"}.
    """
    p = validate_tcp_params(params)
    D = np.asarray(dose_gy_isoe, dtype=float).ravel()
    w = np.asarray(dvol, dtype=float).ravel()
    if D.size != w.size:
        raise ValueError("dose y dvol de distinta longitud")
    if D.size == 0:
        raise ValueError("DVH tumoral vacio")
    if not np.all(np.isfinite(D)) or not np.all(np.isfinite(w)):
        raise ValueError("DVH tumoral con valores no finitos (NaN/Inf)")
    if np.any(D < 0) or np.any(w < 0):
        raise ValueError("DVH tumoral con dosis o volumenes negativos")
    V = float(tumor_volume_cc)
    if not math.isfinite(V) or V <= 0:
        raise ValueError("volumen tumoral debe ser > 0 cm^3 (got %r)" % V)
    s = w.sum()
    if not np.isclose(s, 1.0, atol=1e-6):
        raise ValueError("dvol no normalizado (sum=%g)" % s)

    gr = catcheside_g(p["tof"], p["tos"], tir, p["pf_g"])
    nf = float(p["fracciones"])
    exponent = nf / float(p["c2"])
    # exp(-Ar*D - Gr*Br*D^2)^exponent = exp(exponent * (-Ar*D - Gr*Br*D^2))
    h = np.exp(exponent * (-p["Ar"] * D - gr * p["Br"] * np.square(D)))
    dv_cm3 = w * V
    inner = float(np.sum(h * dv_cm3))
    if inner < 0:
        raise ValueError("sumatorio TCP negativo (entradas invalidas)")
    if inner == 0.0:
        tcp = 1.0
    else:
        tcp = math.exp(-float(p["c1"]) * math.pow(inner, float(p["c2"])))
    tcp = min(max(tcp, 0.0), 1.0)

    # PEUD: dosis uniforme equivalente (TCPc.m)
    peud = _peud_from_tcp(tcp, p, V, gr)
    return {"tcp": tcp, "peud": peud, "g_factor": gr}


def _peud_from_tcp(tcp, p, V_cc, gr):
    """Puerto de PEUD de TCPc.m; None si no es resoluble."""
    c1, c2 = float(p["c1"]), float(p["c2"])
    Ar, Br = float(p["Ar"]), float(p["Br"])
    if tcp >= 1.0 or tcp <= 0.0:
        return None
    try:
        c = math.log(-math.log(tcp) / (c1 * math.pow(V_cc, c2)))
        disc = Ar * Ar - 4.0 * gr * Br * c
        if disc < 0:
            return None
        return (-Ar + math.sqrt(disc)) / (2.0 * gr * Br)
    except (ValueError, ZeroDivisionError, OverflowError):
        return None
