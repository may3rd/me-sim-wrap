import json
import math
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.unitops.hydraulics import api_rp520_liquid_required_area, api_rp520_two_phase_required_area, api_rp520_vapor_required_area, beggs_brill_pressure_drop, lockhart_martinelli_pressure_drop, minor_loss_pressure_drop, orifice_pressure_drop, pipe_pressure_drop


class HydraulicsTest(unittest.TestCase):
    def test_api_rp520_two_phase_matches_dwsim_psv_capture(self):
        result = api_rp520_two_phase_required_area(1_100_000.0, 100_000.0, 0.0020051724826423986, 0.17275978005510137, 9.179643677024487, 49.87102150345734, 44.51374508131255, 10.0, 0.85, 1.0, 1.0)

        self.assertTrue(math.isclose(result.required_area_in2, 0.01905898823975082, rel_tol=1e-12))
        self.assertEqual(result.standard_orifice, "D")
        self.assertEqual(result.standard_area_in2, 0.11)

    def test_api_rp520_liquid_matches_dwsim_psv_capture(self):
        result = api_rp520_liquid_required_area(1_100_000.0, 100_000.0, 0.0001560675811547515, 640.7480609367764, 0.00020498569099240772, 10.0, 0.85, 1.0)

        self.assertTrue(math.isclose(result.required_area_in2, 0.004803176512159658, rel_tol=1e-12))
        self.assertEqual(result.standard_orifice, "D")
        self.assertEqual(result.standard_area_in2, 0.11)

    def test_api_rp520_vapor_matches_dwsim_psv_capture(self):
        result = api_rp520_vapor_required_area(300.0, 1_100_000.0, 100_000.0, 0.1, 0.9764896861575462, 16.04246, 1.3467422057951541, 10.0, 0.85, 1.0, 1.0)

        self.assertTrue(result.choked)
        self.assertTrue(math.isclose(result.required_area_in2, 0.08013876894324415, rel_tol=1e-12))
        self.assertEqual(result.standard_orifice, "D")
        self.assertEqual(result.standard_area_in2, 0.11)

    def test_api_rp520_vapor_rejects_nonphysical_heat_capacity_ratio(self):
        with self.assertRaises(ValidationError):
            api_rp520_vapor_required_area(300.0, 1_100_000.0, 100_000.0, 0.1, 0.98, 16.04, 1.0, 10.0, 0.85, 1.0, 1.0)

    def test_beggs_brill_matches_dwsim_inclined_two_phase_equations(self):
        result = beggs_brill_pressure_drop(0.05, 100.0, 10.0, 4.5e-5, 0.001, 0.002, 10.0, 800.0, 1e-5, 0.001, 0.03)

        self.assertEqual(result.flow_regime, "Intermittent")
        self.assertTrue(math.isclose(result.liquid_holdup, 0.6859528681345057, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.friction_drop_pa, 41930.32301561977, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.static_drop_pa, 54122.915748215484, rel_tol=1e-12))

    def test_beggs_brill_matches_captured_dwsim_two_phase_pipe(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-two-phase-beggs-brill-pr-eos.json").read_text(encoding="utf-8-sig"))
        pipe = next(item for item in golden["outputs"]["objects_after"] if item["tag"] == "PIPE-001")
        values = {item["property"]: item["value"]["value"] for item in pipe["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-two-phase-beggs-brill-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("SelectedFlowPackage") == "Beggs_Brill")
        diameter_m = float(source.findtext("./Profile/Sections/Section/DI")) * 0.0254
        prefix = "HydraulicSegment,1,Results,1,"
        result = beggs_brill_pressure_drop(
            diameter_m, 2.0, 0.2, 4.5e-5,
            values[prefix + "VolumetricFlowVapor"], values[prefix + "VolumetricFlowLiquid"],
            values[prefix + "DensityVapor"], values[prefix + "DensityLiquid"],
            values[prefix + "ViscosityVapor"], values[prefix + "ViscosityLiquid"],
            values[prefix + "SurfaceTension"],
        )

        self.assertEqual(result.flow_regime, values[prefix + "FlowRegime"])
        self.assertTrue(math.isclose(result.liquid_holdup, values[prefix + "LiquidHoldup"], rel_tol=1e-6))
        self.assertTrue(math.isclose(result.friction_drop_pa, values[prefix + "PressureDropFriction"], rel_tol=1e-6))
        self.assertTrue(math.isclose(result.static_drop_pa, values[prefix + "PressureDropHydrostatic"], rel_tol=1e-6))

    def test_lockhart_martinelli_matches_dwsim_two_phase_equations(self):
        result = lockhart_martinelli_pressure_drop(0.05, 100.0, 10.0, 4.5e-5, 0.001, 0.002, 10.0, 800.0, 1e-5, 0.001)

        self.assertTrue(math.isclose(result.vapor_reynolds, 25464.790894703252, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.liquid_reynolds, 40743.66543152521, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.martinelli_parameter, 17.1853059594934, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.liquid_holdup, 0.6792871801523932, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.friction_drop_pa, 44286.52902360688, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.static_drop_pa, 52593.33333333333, rel_tol=1e-12))

    def test_lockhart_martinelli_matches_captured_dwsim_two_phase_pipe(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-two-phase-lockhart-martinelli-pr-eos.json").read_text(encoding="utf-8-sig"))
        pipe = next(item for item in golden["outputs"]["objects_after"] if item["tag"] == "PIPE-001")
        values = {item["property"]: item["value"]["value"] for item in pipe["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-two-phase-lockhart-martinelli-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("SelectedFlowPackage") == "Lockhart_Martinelli")
        diameter_m = float(source.findtext("./Profile/Sections/Section/DI")) * 0.0254
        prefix = "HydraulicSegment,1,Results,1,"
        result = lockhart_martinelli_pressure_drop(
            diameter_m, 2.0, 0.2, 4.5e-5,
            values[prefix + "VolumetricFlowVapor"], values[prefix + "VolumetricFlowLiquid"],
            values[prefix + "DensityVapor"], values[prefix + "DensityLiquid"],
            values[prefix + "ViscosityVapor"], values[prefix + "ViscosityLiquid"],
        )

        self.assertTrue(math.isclose(result.liquid_holdup, values[prefix + "LiquidHoldup"], rel_tol=1e-12))
        self.assertTrue(math.isclose(result.friction_drop_pa, values[prefix + "PressureDropFriction"], rel_tol=1e-12))
        self.assertTrue(math.isclose(result.static_drop_pa, values[prefix + "PressureDropHydrostatic"], rel_tol=1e-12))

    def test_flange_tap_orifice_matches_dwsim_iso_5167_equations(self):
        result = orifice_pressure_drop(0.2, 0.1, 10.0, 1000.0, 0.001, "flange")

        self.assertTrue(math.isclose(result.reynolds, 63661.82836771071, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.discharge_coefficient, 0.6073745138683921, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.orifice_drop_pa, 2557.7410450297325, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.overall_drop_pa, 1871.5167680391855, rel_tol=1e-12))

    def test_orifice_matches_captured_dwsim_liquid_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-orifice-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["inputs"]["objects_before"]}
        properties = {item["property"]: item["value"]["value"] for item in objects["ORIF-FEED"]["properties"]}
        orifice = {item["property"]: item["value"]["value"] for item in objects["ORIF-001"]["properties"]}
        result = orifice_pressure_drop(
            orifice["PROP_OP_2"] / 1000.0, orifice["PROP_OP_1"] / 1000.0,
            properties["PROP_MS_2"], properties["PROP_MS_5"], properties["PROP_MS_38"], "flange",
        )

        self.assertTrue(math.isclose(result.orifice_drop_pa, orifice["PROP_OP_6"], rel_tol=1e-12))
        self.assertTrue(math.isclose(result.overall_drop_pa, orifice["PROP_OP_5"], rel_tol=1e-12))

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
