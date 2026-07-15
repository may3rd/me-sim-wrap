import json
import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.unitops.hydraulics import minor_loss_pressure_drop, orifice_pressure_drop, pipe_pressure_drop


class HydraulicsTest(unittest.TestCase):
    def test_flange_tap_orifice_matches_dwsim_iso_5167_equations(self):
        result = orifice_pressure_drop(0.2, 0.1, 10.0, 1000.0, 0.001, "flange")

        self.assertTrue(math.isclose(result.reynolds, 63661.82836771071, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.discharge_coefficient, 0.6073745138683921, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.orifice_drop_pa, 2557.7410450297325, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.overall_drop_pa, 1871.5167680391855, rel_tol=1e-12))

    def test_minor_loss_uses_dwsim_loss_coefficient_equation(self):
        self.assertEqual(minor_loss_pressure_drop(0.78, 1000.0, 2.5), 2437.5)
        with self.assertRaises(ValidationError):
            minor_loss_pressure_drop(-0.01, 1000.0, 2.5)

    def test_single_phase_pipe_matches_dwsim_darcy_weibach_equations(self):
        result = pipe_pressure_drop(0.1, 100.0, 10.0, 4.5e-5, 0.01, 1000.0, 0.001)

        self.assertTrue(math.isclose(result.velocity_m_s, 1.2732395447351625, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.reynolds, 127323.95447351626, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.friction_factor, 0.019597695586833035, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.friction_drop_pa, 15885.293708161134, rel_tol=1e-12))
        self.assertEqual(result.static_drop_pa, 98000.0)
        self.assertTrue(math.isclose(result.total_drop_pa, 113885.29370816113, rel_tol=1e-12))

    def test_single_phase_pipe_matches_captured_dwsim_segment_drops(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        pipe = next(item for item in golden["inputs"]["objects_before"] if item["tag"] == "PIPE-1")
        values = {item["property"]: item["value"]["value"] for item in pipe["properties"]}
        friction = static = 0.0
        for index in range(1, 6):
            prefix = f"HydraulicSegment,1,Results,{index},"
            result = pipe_pressure_drop(0.10226, 20.0, 2.0, 4.5e-5, values[prefix + "VolumetricFlowLiquid"], values[prefix + "DensityLiquid"], values[prefix + "ViscosityLiquid"])
            friction += result.friction_drop_pa
            static += result.static_drop_pa

        self.assertTrue(math.isclose(friction, values["PressureDropFriction"], rel_tol=1e-5))
        self.assertTrue(math.isclose(static, values["PressureDropStatic"], rel_tol=1e-5))

    def test_single_phase_pipe_rejects_invalid_geometry_and_properties(self):
        with self.assertRaises(ValidationError):
            pipe_pressure_drop(0.0, 100.0, 0.0, 0.0, 0.01, 1000.0, 0.001)
        with self.assertRaises(ValidationError):
            pipe_pressure_drop(0.1, 100.0, 101.0, 0.0, 0.01, 1000.0, 0.001)
        with self.assertRaises(ValidationError):
            pipe_pressure_drop(0.1, 100.0, 0.0, 0.0, 0.01, 1000.0, 0.0)
