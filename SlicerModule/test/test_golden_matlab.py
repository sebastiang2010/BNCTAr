# -*- coding: utf-8 -*-
"""Validacion del puerto Python contra valores golden de MATLAB.

Golden generado con SlicerModule/test/matlab_ref/gen_golden.m, que usa
isodose.m y TCPc.m de la RAIZ del repo y ntcp_isoE.m (copia de UTCP/utils).
Requisito del spec (seccion 11.12): reproducir MATLAB con las mismas entradas.

Ejecutar:  python SlicerModule/test/test_golden_matlab.py
"""

import json
import math
import os
import sys
import unittest

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))  # SlicerModule/

from utcp_engine import isoeffective_dose
from utcp_tcp import tcp_provenzano
from utcp_ntcp import ntcp_lkb
from utcp_dvh import cumulative_to_differential

GOLDEN = os.path.join(_HERE, "matlab_ref", "golden_utcp.json")


@unittest.skipUnless(os.path.exists(GOLDEN),
                     "golden_utcp.json no existe (ejecutar gen_golden.m)")
class TestGoldenMatlab(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(GOLDEN, encoding="utf-8") as f:
            cls.g = json.load(f)
        dims = tuple(int(d) for d in cls.g["dims"])
        cls.dims = dims
        cr = cls.g["comps_raw"]
        cls.comps = {
            "boro": np.array(cr[0][" data_"] if isinstance(cr[0], dict) else cr[0],
                             dtype=float).reshape(dims, order="F"),
            "thn": np.array(cr[1][" data_"] if isinstance(cr[1], dict) else cr[1],
                            dtype=float).reshape(dims, order="F"),
            "fast": np.array(cr[2][" data_"] if isinstance(cr[2], dict) else cr[2],
                             dtype=float).reshape(dims, order="F"),
            "gamma": np.array(cr[3][" data_"] if isinstance(cr[3], dict) else cr[3],
                              dtype=float).reshape(dims, order="F"),
        }
        cls.p14 = [float(x) for x in cls.g["p14"]]
        tp = [float(x) for x in cls.g["tcp_params"]]
        cls.tcp_params = {"Ar": tp[0], "Br": tp[1], "c1": tp[2], "c2": tp[3],
                          "tof": tp[4], "tos": tp[5], "pf_g": tp[6],
                          "ps_g": tp[7], "pf_bnct": tp[8], "ps_bnct": tp[9],
                          "fracciones": 1}
        cls.b10 = float(cls.g["b10"])
        cls.bratio = float(cls.g["bratio"])
        cls.dv = float(cls.g["dv_voxel"])

    @staticmethod
    def _vec(x):
        """jsonencode envuelve arrays en {'data_': ...} en algunas versiones."""
        if isinstance(x, dict):
            x = x.get("data_") or x.get("data") or next(iter(x.values()))
        return np.asarray(x, dtype=float).ravel()

    def test_isodose_matches_matlab(self):
        for k, t in enumerate(self.g["tir"]):
            with self.subTest(tir=t):
                gold = self._vec(self.g["diso"][k]).reshape(self.dims, order="F")
                got, _ = isoeffective_dose(self.comps, self.p14, self.bratio,
                                           float(t), self.b10)
                self.assertEqual(np.isnan(got).sum(), np.isnan(gold).sum())
                ok = np.isfinite(gold)
                self.assertTrue(np.allclose(got[ok], gold[ok], rtol=1e-9,
                                            atol=1e-12))

    def test_tcp_matches_matlab(self):
        for k, t in enumerate(self.g["tir"]):
            with self.subTest(tir=t):
                diso = self._vec(self.g["diso"][k])
                deq = diso[np.isfinite(diso)]
                n = deq.size
                dvol = np.full(n, 1.0 / n)
                V = n * self.dv
                res = tcp_provenzano(np.sort(deq), dvol, V,
                                     dict(self.tcp_params), tir=float(t))
                # tolerancia 1e-7: DisoE coincide con MATLAB a ~1e-9 relativo,
                # pero TCP = exp(-c1*(...)^c2) amplifica exponencialmente esa
                # diferencia (misma entrada, dif solo en decimos exactos).
                self.assertAlmostEqual(res["tcp"], float(self.g["tcp"][k]),
                                       places=7)
                peud_g = self.g["peud"][k]
                if peud_g is None or (isinstance(peud_g, float) and math.isnan(peud_g)):
                    self.assertIsNone(res["peud"])
                else:
                    # PEUD amplifica un orden mas (ln(-ln(TCP)) + raiz)
                    self.assertAlmostEqual(res["peud"], float(peud_g), places=6)

    def _ntcp_from_dvh(self, Xiso, params):
        centers, dvol = cumulative_to_differential(Xiso,
                                                   100.0 * (1.0 - (np.arange(Xiso.size)) / Xiso.size))
        return ntcp_lkb(centers, dvol, params)

    def test_ntcp_matches_matlab(self):
        cases = [("ntcp_skin", {"TD50": 39.8, "m": 0.17, "a": 0, "ab": 10}),
                 ("ntcp_brain", {"TD50": 13.9049, "m": 0.1766, "a": 10, "ab": 3})]
        for k, t in enumerate(self.g["tir"]):
            diso = self._vec(self.g["diso"][k])
            Xiso = np.sort(diso[np.isfinite(diso)])
            for field, params in cases:
                with self.subTest(tir=t, field=field):
                    got = self._ntcp_from_dvh(Xiso, params)["ntcp"]
                    self.assertAlmostEqual(got, float(self.g[field][k]),
                                           places=9)

    def test_disoe_not_rbe_scaled(self):
        # Gy(IsoE) difiere de cualquier escalar x dosis fisica (spec 11.11)
        t = float(self.g["tir"][0])
        gold = self._vec(self.g["diso"][0]).reshape(self.dims, order="F")
        phys = (self.comps["boro"] * self.b10 * self.bratio
                + self.comps["thn"] + self.comps["fast"] + self.comps["gamma"]) * t
        ok = np.isfinite(gold) & (phys > 0)
        ratios = gold[ok] / phys[ok]
        self.assertGreater(float(ratios.max() / ratios.min()), 1.01)


if __name__ == "__main__":
    unittest.main(verbosity=2)
