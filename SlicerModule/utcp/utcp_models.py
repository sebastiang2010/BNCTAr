# -*- coding: utf-8 -*-
"""utcp_models — Registro (registry) de modelos radiobiologicos TCP/NTCP.

Capa pura de Python + numpy: NO importa Qt ni slicer (testeable fuera de 3D Slicer).

Fuentes:
  * TCP Provenzano: TCPc.m (raiz del repo, Perotti 2022 Radiat. Res. 198) y las
    implementaciones locales de B2_v30.m / B3_v15.m del proyecto UTCP.
    Parametros de referencia B2/B3: c1=22, c2=0.33, alpha=0.34 Gy^-1,
    beta=8.3e-9 Gy^-2, fracciones=1. Son valores de REFERENCIA del modelo
    original, no universales.
  * TCP adicionales: seccion "tcp" de json/models.jsonc (3 modelos).
  * NTCP LKB: ntcp_isoE.m (UTCP, solo lectura) — forma estandar
    NTCP = 0.5*(1+erf(t/sqrt(2))) en AMBAS ramas (a=0 y a>0).
    Endpoints registrados:
      - Mucosa PROVISIONAL: TD50=39.8 Gy(IsoE), m=0.17, a=0, ab=10
        (mucositis >= 3, Gonzalez 2017 / Strigari 2012). NO validado para
        piel/hueso/tronco/medula.
      - Cerebro: TD50=13.9049 Gy(IsoE), m=0.1766, a=10, ab=3
        (somnolencia SS >= 2, Dattoli 2024).
  * NTCP Riback 2026 (piel): arquitectura preparada, calculo BLOQUEADO hasta
    contar parametros validados (evaluar() levanta NotImplementedError).
"""

import copy

# --------------------------------------------------------------------------
# TCP
# --------------------------------------------------------------------------

#: Parametros de referencia del modelo Provenzano usado en B2_v30.m/B3_v15.m.
TCP_PROVENZANO_REF = {
    "name": "TCP Provenzano H&N (ref B2/B3)",
    "kind": "provenzano",
    "Ar": 0.34,          # alpha fotonica de referencia [1/Gy]
    "Br": 8.3e-9,        # beta fotonica de referencia [1/Gy^2]
    "c1": 22.0,
    "c2": 0.33,
    "fracciones": 1,     # B2/B3 usan 1 fraccion
    # Reparacion (Schmid 2010, mismas que el modelo isoE de tumor H&N)
    "tof": 34.632,
    "tos": 1212.1,
    "pf_g": 0.53,
    "ps_g": 0.47,
    "pf_bnct": 0.2,
    "ps_bnct": 0.8,
    "source": "B2_v30.m / B3_v15.m (UTCP) + TCPc.m",
}


def _tcp_from_models_jsonc(entries):
    """Adapta los modelos TCP de models.jsonc (claves c_one/c_two) al registro."""
    out = []
    for m in entries:
        out.append({
            "name": m["name"],
            "kind": "provenzano",
            "Ar": float(m["Ar"]),
            "Br": float(m["Br"]),
            "c1": float(m["c_one"]),
            "c2": float(m["c_two"]),
            "fracciones": 1,
            "tof": float(m["tof"]),
            "tos": float(m["tos"]),
            "pf_g": float(m["pf_g"]),
            "ps_g": float(m["ps_g"]),
            "pf_bnct": float(m["pf_bnct"]),
            "ps_bnct": float(m["ps_bnct"]),
            "source": "models.jsonc",
        })
    return out


# --------------------------------------------------------------------------
# NTCP
# --------------------------------------------------------------------------

NTCP_MUCOSA = {
    "name": "LKB mucositis >= 3 (Gonzalez 2017)",
    "kind": "lkb",
    "TD50": 39.8,        # Gy(IsoE)
    "m": 0.17,
    "a": 0,              # 0 -> Dmax (rama EQD2)
    "ab": 10.0,          # Gy; solo rama a=0
    "endpoint": "mucositis grado >= 3",
    "provisional_para": ["piel", "hueso", "tronco", "medula"],
    "source": "ntcp_isoE.m (UTCP) / models.jsonc ntcp id 1",
}

NTCP_BRAIN = {
    "name": "LKB SS >= 2 (Dattoli 2024)",
    "kind": "lkb",
    "TD50": 13.9049,     # Gy(IsoE)
    "m": 0.1766,
    "a": 10,             # EUD con a=10 (Sud 2021)
    "ab": 3.0,           # Gy; no usado en rama a>0 (se conserva por contrato)
    "endpoint": "sindrome de somnolencia grado >= 2",
    "source": "ntcp_isoE.m (UTCP) / models.jsonc ntcp id 2",
}

#: Plantilla configurable: el usuario debe completar TD50/m/a/ab.
NTCP_GENERIC_TEMPLATE = {
    "name": "LKB configurable",
    "kind": "lkb",
    "TD50": None,
    "m": None,
    "a": 1.0,
    "ab": None,
    "endpoint": "",
    "source": "definido por el usuario",
}

#: Modelo Riback 2026 (piel): NO activo hasta tener parametros validados.
NTCP_RIBACK_SKIN = {
    "name": "Riback 2026 (piel)",
    "kind": "riback",
    "status": "pending",
    "endpoint": "toxicidad cutanea (Riback 2026)",
    "note": ("Parametros pendientes de validacion (N0, k_EUD, k_NTCP, "
             "D_presc, umbral de irradiacion, base de areas). Calculo "
             "bloqueado; no se reemplaza automaticamente por el modelo de "
             "mucosa sin advertencia explicita."),
    "source": "B2_v36.m / config/skin_riback_model3.json (UTCP)",
}


def builtin_ntcp_models():
    """Lista de modelos NTCP incluidos por defecto (copias defensivas)."""
    return [copy.deepcopy(NTCP_MUCOSA),
            copy.deepcopy(NTCP_BRAIN),
            copy.deepcopy(NTCP_GENERIC_TEMPLATE),
            copy.deepcopy(NTCP_RIBACK_SKIN)]


def builtin_tcp_models(tcp_jsonc_entries=None):
    """Lista de modelos TCP: referencia B2/B3 + los de models.jsonc."""
    out = [copy.deepcopy(TCP_PROVENZANO_REF)]
    if tcp_jsonc_entries:
        out.extend(_tcp_from_models_jsonc(tcp_jsonc_entries))
    return out


def get_model_by_name(models, name):
    for m in models:
        if m.get("name") == name:
            return copy.deepcopy(m)
    return None


# --------------------------------------------------------------------------
# Perfiles isoE por tejido (valores por defecto del pipeline UTCP)
# --------------------------------------------------------------------------

#: Modelo isoE por defecto segun tejido (B2_v30/B3_v15 -> isoE_tabla_organos):
#: tumor -> H&N in-vivo; BRAIN -> Brain; resto -> Mucosa (PROVISIONAL).
ISOE_DEFAULT_BY_TISSUE = {
    "tumor": "H&N tumor (in-vivo - Perotti 2022)",
    "brain": "Brain (Dattoli 2024)",
    "default": "Mucosa (Gonzalez 2017)",  # PROVISIONAL
}


def default_isoe_model_for_roi(roi_name, is_tumor, available_names):
    """Devuelve el nombre del modelo isoE por defecto para una ROI.

    Orden de decision (sin indices fijos ni nombres hardcodeados obligatorios):
      1. ROI marcada como tumor -> modelo de tumor in-vivo.
      2. Nombre normalizado que contenga 'BRAIN' -> modelo Brain.
      3. Resto -> Mucosa (PROVISIONAL para tejidos no validados).
    Si el nombre no esta disponible en `available_names`, devuelve el primero.
    """
    if is_tumor:
        want = ISOE_DEFAULT_BY_TISSUE["tumor"]
    else:
        norm = (roi_name or "").upper()
        want = ISOE_DEFAULT_BY_TISSUE["brain"] if "BRAIN" in norm \
            else ISOE_DEFAULT_BY_TISSUE["default"]
    if want in available_names:
        return want
    return available_names[0] if available_names else want


def validate_tcp_params(model):
    """Valida un modelo TCP; devuelve lista de errores (vacía si es válido)."""
    errs = []
    for k in ("Ar", "Br", "c1", "c2"):
        v = model.get(k)
        if not isinstance(v, (int, float)) or v != v or v <= 0:
            errs.append(f"TCP: parametro '{k}' debe ser numerico positivo (={v}).")
    fr = model.get("fracciones", 1)
    if not isinstance(fr, (int, float)) or fr <= 0:
        errs.append(f"TCP: 'fracciones' debe ser positivo (={fr}).")
    return errs


def validate_ntcp_params(model):
    """Valida un modelo NTCP; devuelve lista de errores (vacía si es válido)."""
    kind = model.get("kind")
    if kind == "riback":
        return ["NTCP Riback: no disponible (parametros pendientes de validacion)."]
    errs = []
    td50, m = model.get("TD50"), model.get("m")
    a = model.get("a", 1.0)
    ab = model.get("ab")
    if not isinstance(td50, (int, float)) or td50 <= 0:
        errs.append(f"NTCP '{model.get('name')}': TD50 debe ser > 0 en Gy(IsoE).")
    if not isinstance(m, (int, float)) or m <= 0:
        errs.append(f"NTCP '{model.get('name')}': m debe ser > 0.")
    if a is None or a < 0:
        errs.append(f"NTCP '{model.get('name')}': 'a' debe ser >= 0 (0=Dmax, EUD si a>0).")
    if (a == 0) and (not isinstance(ab, (int, float)) or ab <= 0):
        errs.append(f"NTCP '{model.get('name')}': ab (alpha/beta) debe ser > 0 cuando a=0.")
    return errs
