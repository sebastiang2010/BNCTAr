# -*- coding: utf-8 -*-
"""utcp_time_optimizer — Barrido temporal y optimizacion de UTCP.

Capa pura Python + numpy: NO importa Qt ni slicer.

Para cada tiempo de irradiacion t:
    1. las componentes del tally son TASAS -> escalan linealmente con t
       (el escalado lo hace isoeffective_dose via tir).
    2. se recalcula la dosis fotonica isoefectiva con el modelo de cada ROI
       (dependencia temporal explicita: NUNCA se escala una DisoE fija).
    3. DVH isoefectivos del tumor y de cada OAR incluido.
    4. TCP(t), NTCP_j(t), UTCP(t).

Restricciones dosimetricas: metricas (Dmax/Dmean/Dx%) evaluadas en la
magnitud indicada (fisica, RBE o isoE), independientes de los modelos NTCP
(Gy IsoE). Se reportan el maximo de UTCP sin restricciones y el maximo
factible sujeto a restricciones.
"""

import numpy as np

from utcp_engine import isoeffective_dose, compute_utcp_point, _comps_for


def _scaled(comps, tir, tir_ref):
    """Escalado lineal de las componentes (tasa) de tir_ref a tir minutos.

    Soporta comps plano {k: arr} o mapa por mesh {nmesh: {k: arr}}.
    """
    f = float(tir) / float(tir_ref)
    if "boro" in comps:
        return {k: np.asarray(v, dtype=float) * f for k, v in comps.items()}
    return {mesh: {k: np.asarray(v, dtype=float) * f for k, v in c.items()}
            for mesh, c in comps.items()}


def _metric_value(dose_array, mask, metric, pct=None):
    """Dmax/Dmean/Dx% sobre la mascara de una ROI."""
    vals = dose_array[np.asarray(mask, dtype=bool)]
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


def _dose_for_constraint(comps, tumor, oars, tir, b10_conc, spec):
    """Dosis (fisica/RBE/isoE) escalada a tir para la ROI de la restriccion."""
    kind = (spec.get("quantity") or "isoE").lower()
    roi = spec.get("roi")
    if roi == "tumor" or roi is None:
        roi_def, bratio, p14 = tumor, tumor["bratio"], tumor["isoe_params14"]
        mask = tumor["mask"]
    else:
        roi_def = next((o for o in oars if o["name"] == roi), None)
        if roi_def is None:
            raise ValueError("Restriccion sobre ROI inexistente: %r" % roi)
        bratio, p14, mask = roi_def["bratio"], roi_def["isoe_params14"], roi_def["mask"]
    comps_roi = _comps_for(comps, roi_def)
    if kind == "isoe":
        dose, _ = isoeffective_dose(comps_roi, p14, bratio, tir, b10_conc)
        return dose
    # fisica o RBE: suma lineal de componentes con las tasas escaladas
    boro = np.asarray(comps_roi["boro"], dtype=float) * float(b10_conc) * float(bratio)
    total = boro + np.asarray(comps_roi["thn"], dtype=float) \
        + np.asarray(comps_roi["fast"], dtype=float) \
        + np.asarray(comps_roi["gamma"], dtype=float)
    if kind == "rbe":
        rbe = roi_def.get("rbe")
        if rbe is None:
            raise ValueError("Restriccion RBE sin valores RBE en %r" % roi)
        dose = (rbe["boro"] * boro + rbe["thn"] * comps_roi["thn"]
                + rbe["fast"] * comps_roi["fast"]
                + rbe["gamma"] * comps_roi["gamma"])
    else:
        dose = total
    return np.asarray(dose, dtype=float) * (float(tir) / float(spec.get("tir_ref", 1.0)))


def check_constraints(comps, tumor, oars, tir, b10_conc, constraints):
    """Evalua restricciones a un tiempo dado. Devuelve (ok, detalles)."""
    detalles = []
    ok = True
    for c in constraints or []:
        dose = _dose_for_constraint(comps, tumor, oars, tir, b10_conc, c)
        val = _metric_value(dose, _mask_of(c, tumor, oars),
                            c.get("metric", "Dmax"), c.get("pct"))
        limit = float(c["limit"])
        met = bool(np.isfinite(val) and val <= limit)
        ok = ok and met
        detalles.append({"roi": c.get("roi", "tumor"), "metric": c.get("metric", "Dmax"),
                         "value": float(val) if np.isfinite(val) else None,
                         "limit": limit, "quantity": c.get("quantity", "isoE"),
                         "ok": met})
    return ok, detalles


def _mask_of(spec, tumor, oars):
    if spec.get("roi") in (None, "tumor"):
        return tumor["mask"]
    for o in oars:
        if o["name"] == spec["roi"]:
            return o["mask"]
    raise ValueError("ROI de restriccion inexistente: %r" % spec["roi"])


def sweep_utcp(comps, tir_ref, tumor, oars, tir_values, b10_conc,
               constraints=None):
    """Barrido temporal de TCP/NTCP/UTCP.

    comps: componentes crudas evaluadas a tir_ref (tasa).
    tir_values: lista/array de tiempos [min].
    constraints: lista de dicts {roi, metric, pct?, limit, quantity?, tir_ref?}.

    Devuelve dict {tir_values, tcp, utcp, ntcp: {name: [...]},
                   t_opt, utcp_max, t_feasible, utcp_feasible_max,
                   constraint_report: [{tir, detalles}]}.
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

    for i, t in enumerate(tvals):
        c_scaled = _scaled(comps, t, tir_ref)
        try:
            res = compute_utcp_point(c_scaled, tumor, oars, t, b10_conc)
        except ValueError as exc:
            constraint_report.append({"tir": float(t), "error": str(exc),
                                      "detalles": []})
            continue
        tcp[i] = res["tcp"]
        utcp[i] = res["utcp"]
        for n in oar_names:
            info = res["ntcp"].get(n) or {}
            if info.get("ntcp") is not None:
                ntcp_curves[n][i] = info["ntcp"]
        if constraints:
            ok, detalles = check_constraints(c_scaled, tumor, oars, t,
                                             b10_conc, constraints)
            constraint_report.append({"tir": float(t), "ok": ok,
                                      "detalles": detalles})

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
        "constraint_report": constraint_report,
    }
