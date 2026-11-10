# -*- coding: utf-8 -*-
"""Tests del motor radiobiologico UTCP (sin slicer, sin Qt).

Ejecutar:  python -m unittest discover -s SlicerModule/test -v
o:         python SlicerModule/test/test_utcp_engine.py
(sys.path apunta a SlicerModule para importar utcp_*.)
"""

import math
import os
import sys
import unittest

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "utcp"))  # SlicerModule/utcp/

from utcp_models import NTCP_BRAIN, NTCP_MUCOSA, TCP_PROVENZANO_REF
from utcp_dvh import (
    cumulative_to_differential,
    dose_at_volume,
    eud,
    validate_differential,
    voxel_volume_cc,
    dvh_from_voxels,
)
from utcp_ntcp import ntcp_lkb
from utcp_tcp import catcheside_g, tcp_provenzano
from utcp_engine import isoeffective_dose, compute_utcp_point


def uniform_dvh(d_gy, n=1000):
    return np.full(n, float(d_gy)), np.full(n, 1.0 / n)


class TestDVH(unittest.TestCase):
    def test_voxel_volume_cc(self):
        self.assertAlmostEqual(voxel_volume_cc((0.3 * 10,) * 3), 0.027)

    def test_cumulative_to_differential_normalizes(self):
        D = np.array([0.0, 1.0, 2.0, 3.0])
        Y = np.array([100.0, 75.0, 50.0, 25.0])
        centers, dvol = cumulative_to_differential(D, Y)
        self.assertAlmostEqual(float(dvol.sum()), 1.0, places=9)
        validate_differential(centers, dvol)

    def test_cumulative_handles_repeated_doses(self):
        D = np.array([1.0, 1.0, 2.0])
        Y = np.array([100.0, 80.0, 40.0])
        centers, dvol = cumulative_to_differential(D, Y)
        self.assertAlmostEqual(float(dvol.sum()), 1.0, places=9)

    def test_cumulative_flat_becomes_uniform(self):
        D = np.array([1.0, 2.0, 3.0])
        Y = np.array([50.0, 50.0, 50.0])
        centers, dvol = cumulative_to_differential(D, Y)
        self.assertAlmostEqual(float(dvol.sum()), 1.0, places=9)
        self.assertTrue(np.allclose(dvol, dvol[0]))

    def test_dvh_from_voxels_drops_nan(self):
        vox = np.array([1.0, 2.0, np.nan, 2.0])
        d = dvh_from_voxels(vox)
        self.assertEqual(d["n_voxels_valid"], 3)
        self.assertAlmostEqual(d["coverage"], 0.75)

    def test_eud_uniform(self):
        D, w = uniform_dvh(10.0)
        self.assertAlmostEqual(eud(D, w, 10.0), 10.0, places=6)

    def test_eud_heterogeneous_hot_spot_dominates(self):
        D = np.array([2.0, 10.0])
        w = np.array([0.5, 0.5])
        self.assertGreater(eud(D, w, 10.0), 6.0)
        self.assertLess(eud(D, w, 0.5), 6.0)

    def test_dose_at_volume(self):
        D = np.array([0.0, 10.0, 20.0])
        Y = np.array([100.0, 50.0, 0.0])
        self.assertAlmostEqual(dose_at_volume(D, Y, 50.0), 10.0)


class TestNTCP(unittest.TestCase):
    def test_ntcp_at_td50_is_half(self):
        D, w = uniform_dvh(13.9049)
        res = ntcp_lkb(D, w, NTCP_BRAIN)
        self.assertAlmostEqual(res["ntcp"], 0.5, places=9)

    def test_ntcp_monotone_and_bounded(self):
        prev = -1.0
        for d in (5.0, 10.0, 14.0, 20.0, 40.0):
            D, w = uniform_dvh(d)
            n = ntcp_lkb(D, w, NTCP_BRAIN)["ntcp"]
            self.assertGreater(n, prev)
            self.assertTrue(0.0 <= n <= 1.0)
            prev = n

    def test_ntcp_a0_uses_eqd2_on_dmax(self):
        D = np.array([10.0, 40.0])
        w = np.array([0.9, 0.1])
        res = ntcp_lkb(D, w, NTCP_MUCOSA)
        dmax = 40.0
        eqd2 = dmax * (10.0 + dmax) / (10.0 + 2.0)
        t = (eqd2 - 39.8) / (0.17 * 39.8)
        self.assertAlmostEqual(res["dref_gy"], eqd2, places=9)
        self.assertAlmostEqual(
            res["ntcp"], 0.5 * (1 + math.erf(t / math.sqrt(2))), places=9
        )

    def test_ntcp_incomplete_params_raise(self):
        D, w = uniform_dvh(10.0)
        with self.assertRaises(ValueError):
            ntcp_lkb(D, w, {"TD50": None, "m": 0.17, "a": 0, "ab": 10})

    def test_riback_blocked(self):
        D, w = uniform_dvh(10.0)
        with self.assertRaises(NotImplementedError):
            ntcp_lkb(D, w, {"kind": "riback", "status": "pending"})


class TestTCP(unittest.TestCase):
    def test_catcheside_g_range(self):
        gr = catcheside_g(34.632, 1212.1, 30.0, 0.53)
        self.assertTrue(0.0 < gr < 1.0)

    def test_tcp_monotone_increasing_with_dose(self):
        p = dict(TCP_PROVENZANO_REF)
        tcps = [
            tcp_provenzano(*uniform_dvh(d), 27.3, p) for d in (2.0, 5.0, 10.0, 20.0)
        ]
        vals = [r["tcp"] for r in tcps]
        self.assertLess(vals[0], vals[-1])
        for v in vals:
            self.assertTrue(0.0 <= v <= 1.0)

    def test_tcp_larger_tumor_lower_tcp(self):
        p = dict(TCP_PROVENZANO_REF)
        D, w = uniform_dvh(10.0)
        small = tcp_provenzano(D, w, 5.0, p)["tcp"]
        big = tcp_provenzano(D, w, 50.0, p)["tcp"]
        self.assertGreater(small, big)

    def test_tcp_invalid_volume_raises(self):
        p = dict(TCP_PROVENZANO_REF)
        D, w = uniform_dvh(10.0)
        with self.assertRaises(ValueError):
            tcp_provenzano(D, w, 0.0, p)

    def test_tcp_incomplete_params_raise(self):
        D, w = uniform_dvh(10.0)
        with self.assertRaises(ValueError):
            tcp_provenzano(D, w, 27.3, {"c1": 22})


class TestIsoE(unittest.TestCase):
    def _comps(self):
        g = np.linspace(0.001, 0.01, 8).reshape(2, 2, 2)
        return {"boro": g * 0.5, "thn": g * 0.2, "fast": g * 0.3, "gamma": g * 0.1}

    def _p14(self):
        return [
            0.35,
            0.035,
            3.090888558919538,
            7.034e-17,
            2.468388470867692,
            1.0205e-16,
            2.468388470867692,
            1.0205e-16,
            38.961,
            216.4502,
            0.8,
            0.2,
            0.2,
            0.8,
        ]

    def test_missing_component_raises(self):
        comps = self._comps()
        comps.pop("fast")
        with self.assertRaises(ValueError):
            isoeffective_dose(comps, self._p14(), 3.5, 30.0, 20.0)

    def test_disoe_scales_nonlinearly_with_tir(self):
        comps = self._comps()
        d30, _ = isoeffective_dose(comps, self._p14(), 3.5, 30.0, 20.0)
        d60, _ = isoeffective_dose(comps, self._p14(), 3.5, 60.0, 20.0)
        d90, _ = isoeffective_dose(comps, self._p14(), 3.5, 90.0, 20.0)
        ok = np.isfinite(d30)
        self.assertGreater(np.nanmax(d60), np.nanmax(d30))
        self.assertGreater(np.nanmax(d90), np.nanmax(d60))
        ratio1 = np.nanmax(d60) / np.nanmax(d30)
        ratio2 = np.nanmax(d90) / np.nanmax(d60)
        self.assertFalse(np.isclose(ratio1, ratio2, rtol=1e-3))

    def test_nan_where_no_dose(self):
        comps = self._comps()
        comps["boro"][0, 0, 0] = 0.0
        comps["thn"][0, 0, 0] = 0.0
        comps["fast"][0, 0, 0] = 0.0
        comps["gamma"][0, 0, 0] = 0.0
        diso, info = isoeffective_dose(comps, self._p14(), 3.5, 30.0, 20.0)
        self.assertTrue(np.isnan(diso[0, 0, 0]))
        self.assertEqual(info["n_voxels_nan"], 1)

    def test_units_not_rbe(self):
        # Gy(IsoE) NO es dosis fisica x RBE constante: la relacion DisoE/Dfis
        # cambia con tir (reparacion + termino cuadratico).
        comps = self._comps()
        diso30, _ = isoeffective_dose(comps, self._p14(), 3.5, 30.0, 20.0)
        diso90, _ = isoeffective_dose(comps, self._p14(), 3.5, 90.0, 20.0)
        phys = (
            comps["boro"] * 20.0 * 3.5 + comps["thn"] + comps["fast"] + comps["gamma"]
        )
        r30 = diso30 / (phys * 30.0)
        r90 = diso90 / (phys * 90.0)
        ok = np.isfinite(r30) & np.isfinite(r90) & (phys > 0)
        self.assertTrue(ok.any())
        self.assertFalse(np.allclose(r30[ok], r90[ok], rtol=1e-4))


class TestUTCP(unittest.TestCase):
    def _setup(self):
        g = np.linspace(0.001, 0.008, 27).reshape(3, 3, 3)
        comps = {"boro": g * 0.5, "thn": g * 0.2, "fast": g * 0.3, "gamma": g * 0.1}
        p14_muc = [
            0.35,
            0.035,
            3.090888558919538,
            7.034e-17,
            2.468388470867692,
            1.0205e-16,
            2.468388470867692,
            1.0205e-16,
            38.961,
            216.4502,
            0.8,
            0.2,
            0.2,
            0.8,
        ]
        p14_br = [
            1.0,
            0.5,
            21.718,
            0.8274,
            36.359,
            6.579e-30,
            36.359,
            6.579e-30,
            42.0,
            228.0,
            0.38,
            0.62,
            0.8,
            0.2,
        ]
        tumor = {
            "name": "GTVX",
            "mask": g > 0.004,
            "volume_cc": 1.0,
            "bratio": 3.5,
            "isoe_params14": p14_muc,
            "tcp_params": dict(TCP_PROVENZANO_REF),
        }
        brain = {
            "name": "BRAIN",
            "mask": g <= 0.004,
            "volume_cc": 2.0,
            "bratio": 1.0,
            "isoe_params14": p14_br,
            "ntcp_params": dict(NTCP_BRAIN),
            "include": True,
        }
        return comps, tumor, brain

    def test_utcp_single_oar_equals_tcp_times_factor(self):
        comps, tumor, brain = self._setup()
        res = compute_utcp_point(comps, tumor, [brain], 45.0, 20.0)
        self.assertTrue(0 <= res["utcp"] <= 1)
        self.assertAlmostEqual(
            res["utcp"], res["tcp"] * (1 - res["ntcp"]["BRAIN"]["ntcp"])
        )

    def test_utcp_multi_oar_multiplicative(self):
        comps, tumor, brain = self._setup()
        mucosa = dict(brain, name="MUCOSA", ntcp_params=dict(NTCP_MUCOSA))
        res = compute_utcp_point(comps, tumor, [brain, mucosa], 45.0, 20.0)
        expected = (
            res["tcp"]
            * (1 - res["ntcp"]["BRAIN"]["ntcp"])
            * (1 - res["ntcp"]["MUCOSA"]["ntcp"])
        )
        self.assertAlmostEqual(res["utcp"], expected, places=12)

    def test_excluded_oar_not_in_product(self):
        comps, tumor, brain = self._setup()
        off = dict(brain, include=False)
        res = compute_utcp_point(comps, tumor, [off], 45.0, 20.0)
        self.assertIsNone(res["ntcp"]["BRAIN"]["ntcp"])
        self.assertAlmostEqual(res["utcp"], res["tcp"])

    def test_no_tumor_raises(self):
        comps, tumor, brain = self._setup()
        with self.assertRaises(ValueError):
            compute_utcp_point(comps, None, [brain], 45.0, 20.0)

    def test_empty_tumor_mask_raises(self):
        comps, tumor, brain = self._setup()
        tumor["mask"] = np.zeros_like(tumor["mask"])
        with self.assertRaises(ValueError):
            compute_utcp_point(comps, tumor, [brain], 45.0, 20.0)

    def test_nan_voxels_warned_and_dropped(self):
        comps, tumor, brain = self._setup()
        for k in comps:
            comps[k][0, 0, 0] = 0.0
        res = compute_utcp_point(comps, tumor, [brain], 45.0, 20.0)
        self.assertTrue(any("cobertura" in w.lower() for w in res["warnings"]))


class TestTimeOptimizer(unittest.TestCase):
    def test_sweep_returns_finite_curves(self):
        from utcp_time_optimizer import sweep_utcp

        g = np.linspace(0.001, 0.008, 27).reshape(3, 3, 3)
        comps = {"boro": g * 0.5, "thn": g * 0.2, "fast": g * 0.3, "gamma": g * 0.1}
        p14 = [
            0.35,
            0.035,
            3.090888558919538,
            7.034e-17,
            2.468388470867692,
            1.0205e-16,
            2.468388470867692,
            1.0205e-16,
            38.961,
            216.4502,
            0.8,
            0.2,
            0.2,
            0.8,
        ]
        tumor = {
            "name": "GTVX",
            "mask": g > 0.005,
            "volume_cc": 1.0,
            "bratio": 3.5,
            "isoe_params14": p14,
            "tcp_params": dict(TCP_PROVENZANO_REF),
        }
        brain = {
            "name": "BRAIN",
            "mask": g <= 0.005,
            "volume_cc": 2.0,
            "bratio": 1.0,
            "isoe_params14": p14,
            "ntcp_params": dict(NTCP_BRAIN),
            "include": True,
        }
        res = sweep_utcp(comps, 1.0, tumor, [brain], np.linspace(5.0, 90.0, 18), 20.0)
        self.assertEqual(res["tcp"].size, 18)
        self.assertTrue(np.all(np.isfinite(res["utcp"])))
        self.assertTrue(np.all((res["utcp"] >= 0) & (res["utcp"] <= 1)))
        self.assertGreater(res["t_opt"], 0)

    def test_constraint_feasible_max(self):
        from utcp_time_optimizer import sweep_utcp

        g = np.linspace(0.001, 0.008, 27).reshape(3, 3, 3)
        comps = {"boro": g * 0.5, "thn": g * 0.2, "fast": g * 0.3, "gamma": g * 0.1}
        p14 = [
            0.35,
            0.035,
            3.090888558919538,
            7.034e-17,
            2.468388470867692,
            1.0205e-16,
            2.468388470867692,
            1.0205e-16,
            38.961,
            216.4502,
            0.8,
            0.2,
            0.2,
            0.8,
        ]
        tumor = {
            "name": "GTVX",
            "mask": g > 0.005,
            "volume_cc": 1.0,
            "bratio": 3.5,
            "isoe_params14": p14,
            "tcp_params": dict(TCP_PROVENZANO_REF),
        }
        brain = {
            "name": "BRAIN",
            "mask": g <= 0.005,
            "volume_cc": 2.0,
            "bratio": 1.0,
            "isoe_params14": p14,
            "ntcp_params": dict(NTCP_BRAIN),
            "include": True,
        }
        cons = [{"roi": "BRAIN", "metric": "Dmax", "limit": 6.0, "quantity": "isoE"}]
        res = sweep_utcp(
            comps,
            1.0,
            tumor,
            [brain],
            np.linspace(5.0, 90.0, 18),
            20.0,
            constraints=cons,
        )
        if res["t_feasible"] is not None:
            self.assertLessEqual(res["t_feasible"], res["t_opt"] + 1e-9)

    def test_sweep_matches_compute_utcp_point(self):
        # Regression: el barrido rapido (samples 1D por ROI) debe coincidir
        # exactamente con compute_utcp_point (volumen completo) por t, y las
        # curvas deben ser lineales en t (comps = tasas, sin doble escalado).
        from utcp_time_optimizer import sweep_utcp
        from utcp_engine import compute_utcp_point

        g = np.linspace(0.001, 0.008, 27).reshape(3, 3, 3)
        comps = {"boro": g * 0.5, "thn": g * 0.2, "fast": g * 0.3, "gamma": g * 0.1}
        p14 = [
            0.35,
            0.035,
            3.090888558919538,
            7.034e-17,
            2.468388470867692,
            1.0205e-16,
            2.468388470867692,
            1.0205e-16,
            38.961,
            216.4502,
            0.8,
            0.2,
            0.2,
            0.8,
        ]
        tumor = {
            "name": "GTVX",
            "mask": g > 0.005,
            "volume_cc": 1.0,
            "bratio": 3.5,
            "isoe_params14": p14,
            "tcp_params": dict(TCP_PROVENZANO_REF),
        }
        brain = {
            "name": "BRAIN",
            "mask": g <= 0.005,
            "volume_cc": 2.0,
            "bratio": 1.0,
            "isoe_params14": p14,
            "ntcp_params": dict(NTCP_BRAIN),
            "include": True,
        }
        tvals = [30.0, 45.0, 60.0]
        res = sweep_utcp(comps, 1.0, tumor, [brain], tvals, 20.0)
        for k, t in enumerate(tvals):
            ref = compute_utcp_point(comps, tumor, [brain], t, 20.0)
            self.assertAlmostEqual(res["tcp"][k], ref["tcp"], places=12)
            self.assertAlmostEqual(res["utcp"][k], ref["utcp"], places=12)
            self.assertAlmostEqual(
                res["ntcp"]["BRAIN"][k], ref["ntcp"]["BRAIN"]["ntcp"], places=12
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
