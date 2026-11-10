# -*- coding: utf-8 -*-
"""utcp_time_optimizer — Barrido temporal y optimizacion de UTCP.

Capa pura Python + numpy: NO importa Qt ni slicer.

Para cada tiempo de irradiacion t:
    1. las componentes del tally son TASAS -> la dependencia temporal la
       aplica isoeffective_dose via tir (Dlin ~ t, Dquad ~ t^2). Las comps
       NO se re-escalan manualmente (escalado doble = bug documentado).
    2. se recalcula la dosis fotonica isoefectiva con el modelo de cada ROI
       (dependencia temporal explicita: NUNCA se escala una DisoE fija).
    3. DVH isoefectivos del tumor y de cada OAR incluido.
    4. TCP(t), NTCP_j(t), UTCP(t).

Rendimiento: la evaluacion se restringe a los voxels de la mascara de cada
ROI (arrays 1D); el barrido completo no procesa el volumen entero en cada
tiempo. La matematica es la misma de compute_utcp_point (mismas primitivas
del engine); test de regression: barrido == compute_utcp_point por t.

Restricciones dosimetricas: metricas (Dmax/Dmean/Dx%) evaluadas en la
magnitud indicada (fisica, RBE o isoE), independientes de los modelos NTCP
(Gy IsoE). Se reportan el maximo de UTCP sin restricciones y el maximo
factible sujeto a restricciones.
"""

import numpy as np

from utcp_engine import isoeffective_dose, _comps_for
from utcp_dvh import dvh_from_voxels
from utcp_tcp import tcp_provenzano, validate_tcp_params
from utcp_ntcp import ntcp_eval, validate_ntcp_params


def _mask_indices(comps, roi):
    """Indices planos de la mascara de una ROI (validando shapes)."""
    comps_roi = _comps_for(comps, roi)
    mask = np.asarray(roi["mask"], dtype=bool)
    first = np.asarray(comps_roi["boro"])
    if first.shape != mask.shape:
        raise ValueError(
            "Shape de componentes %r != shape de mascara %r para ROI %r"
            % (first.shape, mask.shape, roi.get("name", "?"))
        )
    return comps_roi, np.flatnonzero(mask.ravel())


def _flat_comps(comps_roi, idx):
    """Componentes restringidas a los voxels de la ROI (arrays 1D)."""
    return {k: np.asarray(v, dtype=float).ravel()[idx] for k, v in comps_roi.items()}


def _scaled(comps, tir, tir_ref):
    """DEPRECIADO: las comps son tasas; el escalado temporal lo hace el
    motor via tir. Se conserva solo por compatibilidad de API."""
    f = float(tir) / float(tir_ref)
    if "boro" in comps:
        return {k: np.asarray(v, dtype=float) * f for k, v in comps.items()}
    return {
        mesh: {k: np.asarray(v, dtype=float) * f for k, v in c.items()}
        for mesh, c in comps.items()
    }


def _metric_value(vals, metric, pct=None):
    """Dmax/Dmean/Dx% sobre un sample ya enmascarado (1D, NaNs ok)."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    m = (metric or "Dmax").lower()
    if m == "dmax":
        return float(vals.max())
    if m == "dmean":
        return float(vals.mean())
    if m == "dmin":
        return float(vals.min())
    if m == "dx":
        p = float(pct if pct is not None else 2.0)
        return float(np.percentile(vals, 100.0 - p))
    raise ValueError("Metrica de restriccion desconocida: %r" % metric)


def _roi_dose(comps, roi, tir, b10_conc, quantity):
    """Dosis (1D, voxels de la ROI) en la magnitud pedida a tiempo tir.

    comps son TASAS: fisica/RBE = suma de componentes * tir; isoE delega en
    isoeffective_dose (que aplica tir / tir^2 internamente).
    """
    comps_roi, idx = _mask_indices(comps, roi)
    flat = _flat_comps(comps_roi, idx)
    kind = (quantity or "isoE").lower()
    if kind == "isoe":
        dose, _ = isoeffective_dose(
            flat, roi["isoe_params14"], roi["bratio"], tir, b10_conc
        )
        return np.asarray(dose, dtype=float)
    boro = flat["boro"] * float(b10_conc) * float(roi["bratio"])
    thn = flat["thn"]
    fast = flat["fast"]
    gam = flat["gamma"]
    if kind == "rbe":
        rbe = roi.get("rbe")
        if rbe is None:
            raise ValueError(
                "Restriccion RBE sin valores RBE en %r" % roi.get("name", "?")
            )
        dose = (
            rbe["boro"] * boro
            + rbe["thn"] * thn
            + rbe["fast"] * fast
            + rbe["gamma"] * gam
        )
    elif kind == "fisica":
        dose = boro + thn + fast + gam
    else:
        raise ValueError("Cantidad de restriccion desconocida: %r" % kind)
    return np.asarray(dose, dtype=float) * float(tir)


def _roi_by_name(roi_spec, tumor, oars):
    if roi_spec in (None, "tumor"):
        return tumor
    for o in oars or []:
        if o["name"] == roi_spec:
            return o
    raise ValueError("ROI de restriccion inexistente: %r" % roi_spec)


def check_constraints(comps, tumor, oars, tir, b10_conc, constraints):
    """Evalua restricciones a un tiempo dado. Devuelve (ok, detalles)."""
    detalles = []
    ok = True
    for c in constraints or []:
        roi = _roi_by_name(c.get("roi"), tumor, oars)
        vals = _roi_dose(comps, roi, tir, b10_conc, c.get("quantity", "isoE"))
        val = _metric_value(vals, c.get("metric", "Dmax"), c.get("pct"))
        limit = float(c["limit"])
        met = bool(np.isfinite(val) and val <= limit)
        ok = ok and met
        detalles.append(
            {
                "roi": c.get("roi", "tumor"),
                "metric": c.get("metric", "Dmax"),
                "value": float(val) if np.isfinite(val) else None,
                "limit": limit,
                "quantity": c.get("quantity", "isoE"),
                "ok": met,
            }
        )
    return ok, detalles


def _prep_point(comps, tumor, oars, b10_conc):
    """Precomputa por ROI los samples 1D de componentes (invariantes con t)."""
    pre = {"tumor": None, "oars": {}}
    comps_roi, idx = _mask_indices(comps, tumor)
    pre["tumor"] = {
        "comps": _flat_comps(comps_roi, idx),
        "p14": tumor["isoe_params14"],
        "bratio": tumor["bratio"],
        "volume_cc": tumor.get("volume_cc"),
        "tcp_params": tumor.get("tcp_params"),
        "name": tumor.get("name", "?"),
    }
    for oar in oars or []:
        if not oar.get("include", False):
            continue
        comps_roi, idx = _mask_indices(comps, oar)
        pre["oars"][oar.get("name", "?")] = {
            "comps": _flat_comps(comps_roi, idx),
            "p14": oar["isoe_params14"],
            "bratio": oar["bratio"],
            "volume_cc": oar.get("volume_cc"),
            "ntcp_params": oar.get("ntcp_params"),
            "name": oar.get("name", "?"),
        }
    return pre


def _eval_point(pre, t, b10_conc):
    """TCP/NTCP/UTCP a tiempo t sobre los samples precomputados.

    Misma matematica que compute_utcp_point (mismas primitivas), pero sobre
    arrays 1D ya enmascarados: O(voxels de ROI) por t en vez de O(volumen).
    """
    warnings = []
    tu = pre["tumor"]
    if tu["tcp_params"] is None:
        raise ValueError("Sin parametros TCP para la ROI tumoral %r" % tu["name"])
    validate_tcp_params(tu["tcp_params"])
    diso_t, info_t = isoeffective_dose(
        tu["comps"], tu["p14"], tu["bratio"], t, b10_conc
    )
    try:
        dvh_t = dvh_from_voxels(np.asarray(diso_t, dtype=float), tu["volume_cc"])
    except ValueError as exc:
        raise ValueError(
            "Tumor %r sin dosis isoefectiva valida: %s" % (tu["name"], exc)
        )
    if dvh_t["coverage"] < 1.0:
        warnings.append(
            "Tumor %s: cobertura dosimetrica %.1f%% (%d voxels "
            "NaN descartados)"
            % (tu["name"], 100.0 * dvh_t["coverage"], dvh_t["n_voxels_nan"])
        )
    res_t = tcp_provenzano(
        dvh_t["dose_centers"],
        dvh_t["dvol"],
        tu["volume_cc"],
        tu["tcp_params"],
        tir=float(t),
    )
    tcp = res_t["tcp"]
    ntcp_by_name = {}
    utcp = tcp
    coverage = {"tumor": info_t["coverage"]}
    for name, oar in pre["oars"].items():
        try:
            validate_ntcp_params(oar["ntcp_params"])
        except (ValueError, NotImplementedError) as exc:
            ntcp_by_name[name] = {"included": True, "ntcp": None, "error": str(exc)}
            warnings.append("OAR %s excluido del producto UTCP: %s" % (name, exc))
            continue
        diso_o, info_o = isoeffective_dose(
            oar["comps"], oar["p14"], oar["bratio"], t, b10_conc
        )
        coverage[name] = info_o["coverage"]
        try:
            dvh_o = dvh_from_voxels(np.asarray(diso_o, dtype=float), oar["volume_cc"])
        except ValueError as exc:
            warnings.append("OAR %s sin dosis isoefectiva valida: %s" % (name, exc))
            ntcp_by_name[name] = {"included": True, "ntcp": None, "error": str(exc)}
            continue
        if dvh_o["coverage"] < 1.0:
            warnings.append(
                "OAR %s: cobertura %.1f%% (%d voxels NaN)"
                % (name, 100.0 * dvh_o["coverage"], dvh_o["n_voxels_nan"])
            )
        res_o = ntcp_eval(dvh_o["dose_centers"], dvh_o["dvol"], oar["ntcp_params"])
        ntcp_by_name[name] = {
            "included": True,
            "ntcp": res_o["ntcp"],
            "dref_gy": res_o["dref_gy"],
            "model": oar["ntcp_params"].get("name", "LKB"),
        }
        utcp *= 1.0 - res_o["ntcp"]
    utcp = min(max(utcp, 0.0), 1.0)
    return {
        "tir": float(t),
        "tcp": tcp,
        "ntcp": ntcp_by_name,
        "utcp": utcp,
        "warnings": warnings,
        "coverage": coverage,
    }


def sweep_utcp(
    comps, tir_ref, tumor, oars, tir_values, b10_conc, constraints=None, progress=None
):
    """Barrido temporal de TCP/NTCP/UTCP.

    comps: componentes crudas (TASAS) evaluadas en la escena; la dependencia
        temporal la aplica el motor via tir (tir_ref se ignora y se conserva
        por compatibilidad de firma).
    tir_values: lista/array de tiempos [min].
    constraints: lista de dicts {roi, metric, pct?, limit, quantity?}.
    progress: callable progreso(i, total) llamado en cada paso (opcional).

    Devuelve dict {tir_values, tcp, utcp, ntcp: {name: [...]},
                   t_opt, utcp_max, t_feasible, utcp_feasible_max,
                   warnings, coverage, constraint_report: [{tir, detalles}]}.
    """
    tvals = np.asarray(tir_values, dtype=float).ravel()
    if tvals.size < 2:
        raise ValueError("El barrido requiere >= 2 tiempos")
    if np.any(tvals <= 0):
        raise ValueError("Los tiempos del barrido deben ser > 0")

    tcp = np.zeros(tvals.size)
    utcp = np.zeros(tvals.size)
    oar_names = [o["name"] for o in (oars or []) if o.get("include", False)]
    ntcp_curves = {n: np.full(tvals.size, np.nan) for n in oar_names}
    constraint_report = []
    warnings_all = []
    coverage = {}
    pre = _prep_point(comps, tumor, oars, b10_conc)

    for i, t in enumerate(tvals):
        try:
            res = _eval_point(pre, t, b10_conc)
        except ValueError as exc:
            constraint_report.append(
                {"tir": float(t), "error": str(exc), "detalles": []}
            )
            if progress is not None:
                progress(i, tvals.size)
            continue
        tcp[i] = res["tcp"]
        utcp[i] = res["utcp"]
        if i == 0:
            warnings_all = list(res["warnings"])
            coverage = dict(res["coverage"])
        for n in oar_names:
            info = res["ntcp"].get(n) or {}
            if info.get("ntcp") is not None:
                ntcp_curves[n][i] = info["ntcp"]
        if constraints:
            ok, detalles = check_constraints(
                comps, tumor, oars, t, b10_conc, constraints
            )
            constraint_report.append({"tir": float(t), "ok": ok, "detalles": detalles})
        if progress is not None:
            progress(i, tvals.size)

    i_max = int(np.nanargmax(utcp))
    t_opt = float(tvals[i_max])
    utcp_max = float(utcp[i_max])

    t_feasible, utcp_feas = None, None
    if constraints:
        feas = np.array([d.get("ok", False) for d in constraint_report])
        if feas.any():
            idx = np.where(feas, utcp, -np.inf)
            j = int(np.nanargmax(np.where(np.isfinite(idx), idx, -np.inf)))
            if np.isfinite(idx[j]):
                t_feasible = float(tvals[j])
                utcp_feas = float(utcp[j])
    return {
        "tir_values": tvals,
        "tcp": tcp,
        "utcp": utcp,
        "ntcp": ntcp_curves,
        "t_opt": t_opt,
        "utcp_max": utcp_max,
        "t_feasible": t_feasible,
        "utcp_feasible_max": utcp_feas,
        "warnings": warnings_all,
        "coverage": coverage,
        "constraint_report": constraint_report,
    }
