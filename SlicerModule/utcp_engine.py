# -*- coding: utf-8 -*-
"""utcp_engine — Dosis fotonica isoefectiva y UTCP multi-organo.

Capa pura Python + numpy: NO importa Qt ni slicer.

Puerto FIEL de isodose.m (BNCTAr):

    parameters (14) = [alpha_r, beta_r, alpha_boro, beta_boro, alpha_thn,
                       beta_thn, alpha_fast, beta_fast, tof, tos,
                       pf_g, ps_g, pf_bnct, ps_bnct]

    Gr (t=30 min): Gr = Gs - pf_g*(Gs - Gf)     [f_g=1, f_bnct=0]
    A tir: G_11 = pf_bnct*Gf + ps_bnct*Gs (SUMA); G_22 = G_11;
           G_33 = pf_g*Gf + ps_g*Gs
           G_12 = Gs - (a_12*pf_bnct + a_21*pf_bnct)*(Gs-Gf)
           G_23 = Gs - (a_23*pf_bnct + a_32*pf_g)*(Gs-Gf)
           G_31 = Gs - (a_31*pf_g   + a_13*pf_bnct)*(Gs-Gf)
    Dlin  = (ab*(D.boro*bratio) + at*D.thn + af*D.fast + ar*D.g)*tir
    Dquad = (G_11*bb*(D.boro*bratio)^2 + G_22*bt*D.thn^2 + G_22*bf*D.fast^2
           + G_33*br*D.g^2
           + 2G_12*sqrt(bb*bt)*(D.boro*bratio)*D.thn
           + 2G_12*sqrt(bb*bf)*(D.boro*bratio)*D.fast
           + 2G_31*sqrt(bb*br)*(D.boro*bratio)*D.g
           + 2G_11*sqrt(bt*bf)*D.thn*D.fast
           + 2G_23*sqrt(bt*br)*D.thn*D.g
           + 2G_23*sqrt(bf*br)*D.fast*D.g) * tir^2
    DisoE = 0.5*(ar/br_r/Gr)*(sqrt(1 + (4*Gr/ar/(ar/br_r))*(Dlin+Dquad)) - 1)

Contrato de componentes: los arrays de entrada son el tally crudo PHITS
(semantica de tasa); el boro se multiplica por b10_conc [ppm] AQUI (a diferencia
del path DVH de BNCT.py, que no lo hace: divergencia documentada).
"""

import math

import numpy as np

from utcp_tcp import catcheside_g, tcp_provenzano, validate_tcp_params
from utcp_ntcp import ntcp_eval, validate_ntcp_params
from utcp_dvh import dvh_from_voxels


def _g_pair(tof, tos, tir):
    """Gf, Gs a un tiempo de irradiacion dado (isodose.m, xf=tof/tir)."""
    t = float(tir)
    if t <= 0:
        raise ValueError("tir debe ser > 0 minutos")
    xf = float(tof) / t
    xs = float(tos) / t
    Gf = 2.0 * xf * (1.0 - xf * (1.0 - math.exp(-1.0 / xf)))
    Gs = 2.0 * xs * (1.0 - xs * (1.0 - math.exp(-1.0 / xs)))
    return Gf, Gs


def isoeffective_dose(comps, params14, bratio, tir, b10_conc):
    """Dosis fotonica isoefectiva [Gy IsoE] desde las 4 componentes.

    comps: dict con claves 'boro', 'thn', 'fast', 'gamma' (arrays del tally
        crudo; boro SIN multiplicar por ppm).
    params14: secuencia de 14 parametros (orden isodose.m).
    bratio: concentracion boro tejido/sangre de la ROI.
    tir: minutos de irradiacion.
    b10_conc: concentracion de B-10 [ppm].

    Devuelve (diso_array, info) donde info reporta cobertura: los voxels con
    suma de componentes nula (Dt <= 0) quedan como NaN (comportamiento
    MATLAB: las fracciones a_ij dan NaN) y se cuentan en info.
    """
    p = [float(x) for x in params14]
    if len(p) != 14:
        raise ValueError("params14 requiere 14 valores (orden isodose.m)")
    (ar, br, ab, bb, at, bt, af, bf, tof, tos,
     pf_g, ps_g, pf_bnct, ps_bnct) = p
    if ar <= 0 or br <= 0:
        raise ValueError("alpha_r/beta_r deben ser > 0")
    if not math.isfinite(tir) or tir <= 0:
        raise ValueError("tir debe ser > 0 minutos")

    keys = ("boro", "thn", "fast", "gamma")
    missing = [k for k in keys if comps.get(k) is None]
    if missing:
        raise ValueError(
            "Componentes dosimetricas faltantes: %s. No se sustituyen por "
            "otras distribuciones." % ", ".join(missing))

    boro = np.asarray(comps["boro"], dtype=float) * float(b10_conc) * float(bratio)
    thn = np.asarray(comps["thn"], dtype=float)
    fast = np.asarray(comps["fast"], dtype=float)
    gam = np.asarray(comps["gamma"], dtype=float)
    shapes = {arr.shape for arr in (boro, thn, fast, gam)}
    if len(shapes) != 1:
        raise ValueError("Componentes con shapes distintas: %r" % (shapes,))

    dt = boro + thn + fast + gam
    with np.errstate(divide="ignore", invalid="ignore"):
        fb = np.where(dt > 0, boro / dt, np.nan)
        fn = np.where(dt > 0, (thn + fast) / dt, np.nan)
        fg = np.where(dt > 0, gam / dt, np.nan)
        a12 = fb / (fb + fn)
        a21 = fn / (fb + fn)
        a13 = fb / (fb + fg)
        a31 = fg / (fb + fg)
        a23 = fn / (fn + fg)
        a32 = fg / (fn + fg)

    Gf30, Gs30 = _g_pair(tof, tos, 30.0)
    gr = Gs30 - pf_g * (Gs30 - Gf30)
    Gf, Gs = _g_pair(tof, tos, tir)
    span = Gs - Gf

    g11 = pf_bnct * Gf + ps_bnct * Gs   # suma, NO lc
    g22 = g11
    g33 = pf_g * Gf + ps_g * Gs
    g12 = Gs - (a12 * pf_bnct + a21 * pf_bnct) * span
    g23 = Gs - (a23 * pf_bnct + a32 * pf_g) * span
    g31 = Gs - (a31 * pf_g + a13 * pf_bnct) * span

    dlin = (ab * boro + at * thn + af * fast + ar * gam) * tir
    dquad = (
        g11 * bb * np.square(boro)
        + g22 * bt * np.square(thn)
        + g22 * bf * np.square(fast)
        + g33 * br * np.square(gam)
        + 2.0 * g12 * math.sqrt(bb * bt) * boro * thn
        + 2.0 * g12 * math.sqrt(bb * bf) * boro * fast
        + 2.0 * g31 * math.sqrt(bb * br) * boro * gam
        + 2.0 * g11 * math.sqrt(bt * bf) * thn * fast
        + 2.0 * g23 * math.sqrt(bt * br) * thn * gam
        + 2.0 * g23 * math.sqrt(bf * br) * fast * gam
    ) * tir * tir

    abr = ar / br
    inside = 1.0 + (4.0 * gr / ar / abr) * (dlin + dquad)
    inside = np.where(dt > 0, inside, np.nan)
    with np.errstate(invalid="ignore"):
        diso = 0.5 * (abr / gr) * (np.sqrt(np.maximum(inside, 0.0)) - 1.0)
    diso = np.where(dt > 0, diso, np.nan)

    n_total = int(diso.size)
    n_nan = int(np.isnan(diso).sum())
    info = {
        "n_voxels": n_total,
        "n_voxels_nan": n_nan,
        "coverage": (n_total - n_nan) / n_total if n_total else 0.0,
        "g_factor": gr,
    }
    return diso, info


def roi_sample(dose_array, mask):
    """Muestra 1D de voxels de una ROI: arr[mask] (NaNs conservados)."""
    arr = np.asarray(dose_array, dtype=float)
    m = np.asarray(mask, dtype=bool)
    if arr.shape != m.shape:
        raise ValueError(
            "dose %r y mask %r con shapes distintas" % (arr.shape, m.shape))
    return arr[m]


def _comps_for(comps, roi):
    """Selecciona el set de componentes de un ROI.

    comps puede ser un unico dict {'boro','thn','fast','gamma'} (un solo mesh)
    o un mapa {nmesh: comps} cuando ROIs de distinto mesh tienen Kerma propios.
    """
    if "boro" in comps:
        return comps
    key = int(roi.get("mesh", 0) or 0)
    if key not in comps:
        raise ValueError("Sin componentes dosimetricas para mesh %d (ROI %r)"
                         % (key, roi.get("name", "?")))
    return comps[key]


def compute_utcp_point(comps, tumor, oars, tir, b10_conc):
    """TCP + NTCP_j + UTCP a un tiempo de irradiacion dado.

    comps: dict {'boro','thn','fast','gamma'} tally crudo (tasa), o mapa
        {nmesh: comps} con componentes por mesh.
    tumor: dict {name, mask, volume_cc, bratio, isoe_params14, tcp_params,
        mesh?}.
    oars: lista de dicts {name, mask, volume_cc?, bratio, isoe_params14,
           ntcp_params, include(bool), provisional?, mesh?}.
    tir: minutos; b10_conc: ppm.

    UTCP = TCP * prod(1 - NTCP_j) sobre los OARs con include=True.
    Supuesto documentado: independencia condicional entre complicaciones.

    Devuelve dict {tir, tcp, ntcp: {name: {...}}, utcp, warnings, coverage}.
    """
    warnings = []
    if tumor is None:
        raise ValueError("No hay ROI marcada como tumor: TCP/UTCP indefinidos")
    tcp_params = tumor.get("tcp_params")
    if tcp_params is None:
        raise ValueError("Sin parametros TCP para la ROI tumoral %r"
                         % tumor.get("name", "?"))
    validate_tcp_params(tcp_params)

    diso_t, info_t = isoeffective_dose(
        _comps_for(comps, tumor), tumor["isoe_params14"],
        tumor["bratio"], tir, b10_conc)
    sample_t = roi_sample(diso_t, tumor["mask"])
    if sample_t.size == 0:
        raise ValueError("Mascara tumoral vacia (%r)" % tumor["name"])
    try:
        dvh_t = dvh_from_voxels(sample_t, tumor.get("volume_cc"))
    except ValueError as exc:
        raise ValueError(
            "Tumor %r sin dosis isoefectiva valida: %s" % (tumor["name"], exc))
    if dvh_t["coverage"] < 1.0:
        warnings.append(
            "Tumor %s: cobertura dosimetrica %.1f%% (%d voxels NaN descartados)"
            % (tumor["name"], 100.0 * dvh_t["coverage"], dvh_t["n_voxels_nan"]))
    res_t = tcp_provenzano(dvh_t["dose_centers"], dvh_t["dvol"],
                           tumor["volume_cc"], tcp_params, tir=tir)
    tcp = res_t["tcp"]

    ntcp_by_name = {}
    utcp = tcp
    coverage = {"tumor": info_t["coverage"]}
    for oar in oars or []:
        name = oar.get("name", "?")
        if not oar.get("include", False):
            ntcp_by_name[name] = {"included": False, "ntcp": None}
            continue
        try:
            validate_ntcp_params(oar["ntcp_params"])
        except (ValueError, NotImplementedError) as exc:
            ntcp_by_name[name] = {"included": True, "ntcp": None,
                                  "error": str(exc)}
            warnings.append("OAR %s excluido del producto UTCP: %s"
                            % (name, exc))
            continue
        diso_o, info_o = isoeffective_dose(
            _comps_for(comps, oar), oar["isoe_params14"], oar["bratio"],
            tir, b10_conc)
        coverage[name] = info_o["coverage"]
        sample_o = roi_sample(diso_o, oar["mask"])
        if sample_o.size == 0:
            warnings.append("OAR %s con mascara vacia: excluido" % name)
            ntcp_by_name[name] = {"included": True, "ntcp": None,
                                  "error": "mascara vacia"}
            continue
        try:
            dvh_o = dvh_from_voxels(sample_o, oar.get("volume_cc"))
        except ValueError as exc:
            warnings.append("OAR %s sin dosis isoefectiva valida: %s"
                            % (name, exc))
            ntcp_by_name[name] = {"included": True, "ntcp": None,
                                  "error": str(exc)}
            continue
        if dvh_o["coverage"] < 1.0:
            warnings.append(
                "OAR %s: cobertura %.1f%% (%d voxels NaN)"
                % (name, 100.0 * dvh_o["coverage"], dvh_o["n_voxels_nan"]))
        res_o = ntcp_eval(dvh_o["dose_centers"], dvh_o["dvol"],
                          oar["ntcp_params"])
        ntcp_by_name[name] = {
            "included": True,
            "ntcp": res_o["ntcp"],
            "dref_gy": res_o["dref_gy"],
            "model": oar["ntcp_params"].get("name", "LKB"),
        }
        if oar["ntcp_params"].get("provisional"):
            warnings.append(
                "OAR %s usa modelo PROVISIONAL (%s): NTCP orientativa"
                % (name, ntcp_by_name[name]["model"]))
        utcp *= (1.0 - res_o["ntcp"])

    utcp = min(max(utcp, 0.0), 1.0)
    return {
        "tir": float(tir),
        "tcp": tcp,
        "tcp_peud": res_t["peud"],
        "ntcp": ntcp_by_name,
        "utcp": utcp,
        "warnings": warnings,
        "coverage": coverage,
        "assumption": ("UTCP = TCP * prod(1-NTCP_j); supone independencia "
                       "condicional entre complicaciones"),
    }
