import json
import math
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.thermo.flash import flash_enthalpy, tp_flash
from mesim.thermo.ideal import load_correlations
from mesim.unitops.hydraulics import TwoPhasePipeSegment, api_rp520_liquid_required_area, api_rp520_two_phase_required_area, api_rp520_vapor_required_area, beggs_brill_pressure_drop, beggs_brill_pressure_drop_profile, dwsim_terrain_thermal_conductivity, liquid_pipe_supplied_state_profile, lockhart_martinelli_pressure_drop, lockhart_martinelli_pressure_drop_profile, minor_loss_pressure_drop, orifice_pressure_drop, pipe_absorbed_solar_radiation, pipe_defined_heat_pr_profile, pipe_defined_htc_gradient_profile, pipe_defined_htc_heat_transfer, pipe_defined_htc_profile, pipe_estimated_htc_air, pipe_estimated_htc_air_pr_gradient_profile, pipe_estimated_htc_air_profile, pipe_estimated_htc_soil, pipe_estimated_htc_soil_pr_profile, pipe_estimated_htc_water, pipe_estimated_htc_water_pr_profile, pipe_irradiated_heat_transfer, pipe_pressure_drop, pipe_pressure_drop_profile, pipe_solar_irradiation_source


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
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        values = {item["property"]: item["value"]["value"] for item in objects["PIPE-001"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["TPPIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["TPPIPE-PROD"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-two-phase-beggs-brill-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("SelectedFlowPackage") == "Beggs_Brill")
        section = source.find("./Profile/Sections/Section")
        diameter_m = float(section.findtext("DI")) * 0.0254
        increments = int(section.findtext("Incrementos"))
        profile = beggs_brill_pressure_drop_profile(
            feed["PROP_MS_1"], diameter_m, 4.5e-5,
            tuple(TwoPhasePipeSegment(
                float(section.findtext("Comprimento")) / increments,
                float(section.findtext("Elevacao")) / increments,
                values[f"HydraulicSegment,1,Results,{index},VolumetricFlowVapor"],
                values[f"HydraulicSegment,1,Results,{index},VolumetricFlowLiquid"],
                values[f"HydraulicSegment,1,Results,{index},DensityVapor"],
                values[f"HydraulicSegment,1,Results,{index},DensityLiquid"],
                values[f"HydraulicSegment,1,Results,{index},ViscosityVapor"],
                values[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"],
                values[f"HydraulicSegment,1,Results,{index},SurfaceTension"],
            ) for index in range(1, increments + 1)),
        )

        self.assertEqual(len(profile.segment_results), increments)
        for index, result in enumerate(profile.segment_results, 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            self.assertEqual(result.flow_regime, values[prefix + "FlowRegime"])
            self.assertTrue(math.isclose(result.liquid_holdup, values[prefix + "LiquidHoldup"], rel_tol=1e-6))
            self.assertTrue(math.isclose(result.friction_drop_pa, values[prefix + "PressureDropFriction"], rel_tol=1e-6))
            self.assertTrue(math.isclose(result.static_drop_pa, values[prefix + "PressureDropHydrostatic"], rel_tol=1e-6))
        self.assertTrue(math.isclose(profile.friction_drop_pa, values["PressureDropFriction"], rel_tol=1e-6))
        self.assertTrue(math.isclose(profile.static_drop_pa, values["PressureDropStatic"], rel_tol=1e-6))
        self.assertTrue(math.isclose(profile.outlet_pressure_pa, product["PROP_MS_1"], rel_tol=1e-6))

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
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        values = {item["property"]: item["value"]["value"] for item in objects["PIPE-001"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["TPPIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["TPPIPE-PROD"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-two-phase-lockhart-martinelli-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("SelectedFlowPackage") == "Lockhart_Martinelli")
        section = source.find("./Profile/Sections/Section")
        diameter_m = float(section.findtext("DI")) * 0.0254
        increments = int(section.findtext("Incrementos"))
        profile = lockhart_martinelli_pressure_drop_profile(
            feed["PROP_MS_1"], diameter_m, 4.5e-5,
            tuple(TwoPhasePipeSegment(
                float(section.findtext("Comprimento")) / increments,
                float(section.findtext("Elevacao")) / increments,
                values[f"HydraulicSegment,1,Results,{index},VolumetricFlowVapor"],
                values[f"HydraulicSegment,1,Results,{index},VolumetricFlowLiquid"],
                values[f"HydraulicSegment,1,Results,{index},DensityVapor"],
                values[f"HydraulicSegment,1,Results,{index},DensityLiquid"],
                values[f"HydraulicSegment,1,Results,{index},ViscosityVapor"],
                values[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"],
            ) for index in range(1, increments + 1)),
        )

        self.assertEqual(len(profile.segment_results), increments)
        for index, result in enumerate(profile.segment_results, 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            self.assertTrue(math.isclose(result.liquid_holdup, values[prefix + "LiquidHoldup"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.friction_drop_pa, values[prefix + "PressureDropFriction"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.static_drop_pa, values[prefix + "PressureDropHydrostatic"], rel_tol=1e-12))
        self.assertTrue(math.isclose(profile.friction_drop_pa, values["PressureDropFriction"], rel_tol=1e-12))
        self.assertTrue(math.isclose(profile.static_drop_pa, values["PressureDropStatic"], rel_tol=1e-12))
        self.assertTrue(math.isclose(profile.outlet_pressure_pa, product["PROP_MS_1"], rel_tol=1e-12))

    def test_two_phase_pipe_profiles_reject_empty_segments(self):
        with self.assertRaises(ValidationError):
            beggs_brill_pressure_drop_profile(1_100_000.0, 0.05, 4.5e-5, ())
        with self.assertRaises(ValidationError):
            lockhart_martinelli_pressure_drop_profile(1_100_000.0, 0.05, 4.5e-5, ())

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
        objects = {item["tag"]: item for item in golden["inputs"]["objects_before"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        profile = pipe_pressure_drop_profile(
            feed["PROP_MS_1"], 0.10226, 4.5e-5, (20.0,) * 5, (2.0,) * 5,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VolumetricFlowLiquid"] for index in range(1, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(1, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(1, 6)),
        )

        self.assertEqual(len(profile.segment_results), 5)
        for index, result in enumerate(profile.segment_results, 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            self.assertTrue(math.isclose(result.friction_drop_pa, pipe[prefix + "PressureDropFriction"], rel_tol=1e-5))
            self.assertTrue(math.isclose(result.static_drop_pa, pipe[prefix + "PressureDropHydrostatic"], rel_tol=1e-5))
        self.assertTrue(math.isclose(profile.friction_drop_pa, pipe["PressureDropFriction"], rel_tol=1e-5))
        self.assertTrue(math.isclose(profile.static_drop_pa, pipe["PressureDropStatic"], rel_tol=1e-5))
        self.assertTrue(math.isclose(profile.outlet_pressure_pa, product["PROP_MS_1"], rel_tol=1e-5))

    def test_defined_htc_pipe_matches_captured_dwsim_thermal_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        segment_elevation_m = float(section.findtext("Elevacao")) / increments
        overall_htc = pipe["ThermalProfile,OverallHTC"]
        external_temperature = pipe["ThermalProfile,ExternalTemperatureDefinedHTC"]

        for index in range(1, 6):
            prefix = f"HydraulicSegment,1,Results,{index},"
            result = pipe_defined_htc_heat_transfer(
                pipe[prefix + "InitialTemperature"], external_temperature, overall_htc,
                outer_diameter_m, segment_length_m, feed["PROP_MS_2"],
                pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
            )
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[prefix + "HeatTransfer"] * 1_000.0, rel_tol=3e-3))

        profile = liquid_pipe_supplied_state_profile(
            feed["PROP_MS_1"], pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            inner_diameter_m, outer_diameter_m, 4.5e-5,
            (segment_length_m,) * increments, (segment_elevation_m,) * increments,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VolumetricFlowLiquid"] for index in range(1, increments + 1)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(1, increments + 1)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(1, increments + 1)),
            external_temperature, overall_htc, (segment_length_m,) * 4, feed["PROP_MS_2"],
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
        )

        self.assertEqual(len(profile.pressure.segment_results), increments)
        self.assertTrue(math.isclose(profile.pressure.friction_drop_pa, pipe["PressureDropFriction"], rel_tol=1e-5))
        self.assertTrue(math.isclose(profile.pressure.static_drop_pa, pipe["PressureDropStatic"], rel_tol=1e-5))
        self.assertTrue(math.isclose(profile.pressure.outlet_pressure_pa, product["PROP_MS_1"], rel_tol=1e-5))
        self.assertEqual(len(profile.thermal.segment_results), 4)
        self.assertTrue(math.isclose(profile.thermal.total_area_m2, 4.0 * math.pi * outer_diameter_m * segment_length_m, rel_tol=1e-12))
        self.assertTrue(math.isclose(profile.thermal.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.thermal.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3))

    def test_defined_htc_gradient_uses_dwsim_segment_start_positions(self):
        profile = pipe_defined_htc_gradient_profile(
            300.0, 350.0, 0.1, 25.0, 0.1, (10.0, 20.0), 10.0, (2_000.0, 2_000.0),
        )

        self.assertEqual(profile.segment_start_distances_m, (0.0, 10.0))
        self.assertEqual(profile.external_temperatures_k, (350.0, 351.0))
        self.assertTrue(math.isclose(profile.segment_results[0].outlet_temperature_k, 300.19596451359183, rel_tol=1e-12))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, 300.5934156452265, rel_tol=1e-12))
        self.assertTrue(math.isclose(profile.heat_transfer_w, 11868.31290453023, rel_tol=1e-12))

        constant = pipe_defined_htc_profile(300.0, 350.0, 25.0, 0.1, (10.0, 20.0), 10.0, (2_000.0, 2_000.0))
        zero_gradient = pipe_defined_htc_gradient_profile(
            300.0, 350.0, 0.0, 25.0, 0.1, (10.0, 20.0), 10.0, (2_000.0, 2_000.0),
        )
        self.assertEqual(zero_gradient.segment_results, constant.segment_results)
        self.assertEqual(zero_gradient.outlet_temperature_k, constant.outlet_temperature_k)
        self.assertEqual(zero_gradient.heat_transfer_w, constant.heat_transfer_w)

    def test_defined_htc_gradient_matches_captured_dwsim_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-gradient-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-gradient-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        increments = int(section.findtext("Incrementos"))
        segment_length_m = float(section.findtext("Comprimento")) / increments
        segment_elevation_m = float(section.findtext("Elevacao")) / increments
        profile = liquid_pipe_supplied_state_profile(
            feed["PROP_MS_1"], pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            float(section.findtext("DI")) * 0.0254, float(section.findtext("DE")) * 0.0254, 4.5e-5,
            (segment_length_m,) * increments, (segment_elevation_m,) * increments,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VolumetricFlowLiquid"] for index in range(1, increments + 1)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(1, increments + 1)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(1, increments + 1)),
            pipe["ThermalProfile,ExternalTemperatureDefinedHTC"], pipe["ThermalProfile,OverallHTC"],
            (segment_length_m,) * 4, feed["PROP_MS_2"],
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            pipe["ThermalProfile,ExternalTemperatureGradientDefinedHTC"], segment_length_m,
        )

        self.assertEqual(profile.thermal.segment_start_distances_m, (20.0, 40.0, 60.0, 80.0))
        self.assertEqual(profile.thermal.external_temperatures_k, (352.0, 354.0, 356.0, 358.0))
        for index, result in enumerate(profile.thermal.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.pressure.outlet_pressure_pa, product["PROP_MS_1"], rel_tol=1e-5))
        self.assertTrue(math.isclose(profile.thermal.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.thermal.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3))

    def test_defined_heat_pipe_matches_captured_dwsim_case(self):
        root_path = Path(__file__).parents[1]
        golden = json.loads((root_path / "tests/golden/u3-pipe-thermal-defined-heat-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(root_path / "tests/u3-pipe-thermal-defined-heat-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        specified_heat_transfer_w = float(thermal.findtext("Calor_trocado")) * 1_000.0

        self.assertEqual(thermal.findtext("TipoPerfil"), "Definir_Q")
        self.assertEqual(specified_heat_transfer_w, 10_000.0)
        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "false")
        catalog = {compound.id: compound for compound in load_compounds(root_path / "data/compounds/v1.json")}
        compounds = (catalog["N-pentane"], catalog["Ethane"])
        composition = (0.952380952380952, 0.0476190476190476)
        interactions = load_pr_interactions(root_path / "data/interactions/pr-v1.json")
        correlations = load_correlations(root_path / "data/correlations/ideal-v1.json")
        molar_flow_kmol_s = feed["PROP_MS_3"] / 1_000.0
        profile = pipe_defined_heat_pr_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            pipe["HydraulicSegment,1,Results,2,InitialPressure"],
            compounds, composition, interactions, correlations, molar_flow_kmol_s,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},InitialPressure"] for index in range(3, 7)),
            outer_diameter_m, (segment_length_m,) * 4,
            specified_heat_transfer_w, increments, (250.0, 350.0),
        )

        for index, result in enumerate(profile.segment_results, 2):
            self.assertEqual(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0)
            self.assertTrue(math.isclose(
                result.outlet_temperature_k,
                pipe[f"HydraulicSegment,1,Results,{index + 1},InitialTemperature"],
                rel_tol=2e-8,
            ))
        self.assertEqual(profile.heat_transfer_w, 8_000.0)
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=2e-8))

        inlet_flash = tp_flash(compounds, composition, interactions, feed["PROP_MS_0"], feed["PROP_MS_1"])
        outlet_flash = tp_flash(compounds, composition, interactions, product["PROP_MS_0"], product["PROP_MS_1"])
        stream_energy_w = molar_flow_kmol_s * (
            flash_enthalpy(compounds, correlations, outlet_flash)
            - flash_enthalpy(compounds, correlations, inlet_flash)
        )
        self.assertTrue(math.isclose(stream_energy_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=1e-5))
        self.assertGreater(profile.heat_transfer_w - stream_energy_w, 400.0)

    def test_estimated_htc_air_matches_captured_dwsim_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-estimated-htc-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-estimated-htc-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        external_air_velocity = float(thermal.findtext("Velocidade"))

        self.assertEqual(thermal.findtext("TipoPerfil"), "Estimar_CGTC")
        self.assertEqual(thermal.findtext("Meio"), "0")
        self.assertEqual(thermal.findtext("Incluir_cti"), "true")
        self.assertEqual(thermal.findtext("Incluir_paredes"), "true")
        self.assertEqual(thermal.findtext("Incluir_cte"), "true")
        self.assertEqual(thermal.findtext("Incluir_isolamento"), "false")
        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            result = pipe_estimated_htc_air(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                external_temperature, inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], external_air_velocity,
            )
            self.assertTrue(math.isclose(result.internal_htc_w_m2_k, pipe[prefix + "HTCinternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.wall_htc_w_m2_k, pipe[prefix + "HTCpipewall"], rel_tol=5e-6))
            self.assertTrue(math.isclose(result.external_htc_w_m2_k, pipe[prefix + "HTCexternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.overall_htc_w_m2_k, pipe[prefix + "HTCoverall"], rel_tol=1e-7))

        profile = pipe_estimated_htc_air_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"], external_temperature,
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            feed["PROP_MS_2"],
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            external_air_velocity,
        )

        self.assertEqual(len(profile.segment_results), 4)
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3))

    def test_estimated_htc_air_gradient_matches_captured_dwsim_case(self):
        root_path = Path(__file__).parents[1]
        golden = json.loads((root_path / "tests/golden/u3-pipe-thermal-estimated-htc-gradient-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(root_path / "tests/u3-pipe-thermal-estimated-htc-gradient-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        base_external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        external_temperature_gradient = float(thermal.findtext("AmbientTemperatureGradient_EstimateHTC"))
        external_air_velocity = float(thermal.findtext("Velocidade"))

        self.assertEqual(thermal.findtext("TipoPerfil"), "Estimar_CGTC")
        self.assertEqual(thermal.findtext("Meio"), "0")
        self.assertEqual(base_external_temperature, 350.0)
        self.assertEqual(external_temperature_gradient, 0.1)
        self.assertEqual(external_air_velocity, 2.0)
        self.assertEqual(thermal.findtext("Incluir_isolamento"), "false")
        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "false")
        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            result = pipe_estimated_htc_air(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                pipe[prefix + "ExternalTemperature"],
                inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], external_air_velocity,
            )
            self.assertTrue(math.isclose(result.internal_htc_w_m2_k, pipe[prefix + "HTCinternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.wall_htc_w_m2_k, pipe[prefix + "HTCpipewall"], rel_tol=1e-5))
            self.assertTrue(math.isclose(result.external_htc_w_m2_k, pipe[prefix + "HTCexternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.overall_htc_w_m2_k, pipe[prefix + "HTCoverall"], rel_tol=1e-7))

        catalog = {compound.id: compound for compound in load_compounds(root_path / "data/compounds/v1.json")}
        profile = pipe_estimated_htc_air_pr_gradient_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            pipe["HydraulicSegment,1,Results,2,InitialPressure"],
            base_external_temperature, external_temperature_gradient, segment_length_m,
            (catalog["N-pentane"], catalog["Ethane"]),
            (0.952380952380952, 0.0476190476190476),
            load_pr_interactions(root_path / "data/interactions/pr-v1.json"),
            load_correlations(root_path / "data/correlations/ideal-v1.json"),
            feed["PROP_MS_3"] / 1_000.0,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},InitialPressure"] for index in range(3, 7)),
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            external_air_velocity,
        )

        self.assertEqual(profile.segment_start_distances_m, (20.0, 40.0, 60.0, 80.0))
        self.assertEqual(profile.external_temperatures_k, (352.0, 354.0, 356.0, 358.0))
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_estimated_htc_water_matches_captured_dwsim_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-estimated-htc-water-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-estimated-htc-water-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        external_water_velocity = float(thermal.findtext("Velocidade"))
        external_water_density = 993.38984745012783
        external_water_viscosity = 0.00069354225499842137
        external_water_heat_capacity = 4178.7265325108107
        external_water_conductivity = 0.62609581196362363

        self.assertEqual(thermal.findtext("TipoPerfil"), "Estimar_CGTC")
        self.assertEqual(thermal.findtext("Meio"), "1")
        self.assertEqual(external_temperature, 310.0)
        self.assertEqual(external_water_velocity, 2.0)
        self.assertEqual(thermal.findtext("Incluir_cti"), "true")
        self.assertEqual(thermal.findtext("Incluir_paredes"), "true")
        self.assertEqual(thermal.findtext("Incluir_cte"), "true")
        self.assertEqual(thermal.findtext("Incluir_isolamento"), "false")
        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "false")
        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            result = pipe_estimated_htc_water(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], external_water_velocity,
                external_water_heat_capacity, external_water_conductivity,
                external_water_viscosity, external_water_density,
            )
            self.assertTrue(math.isclose(result.internal_htc_w_m2_k, pipe[prefix + "HTCinternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.wall_htc_w_m2_k, pipe[prefix + "HTCpipewall"], rel_tol=1e-5))
            self.assertTrue(math.isclose(result.external_htc_w_m2_k, pipe[prefix + "HTCexternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.overall_htc_w_m2_k, pipe[prefix + "HTCoverall"], rel_tol=5e-6))
        self.assertTrue(math.isclose(result.external_reynolds, 327433.60641469847, rel_tol=1e-12))
        self.assertTrue(math.isclose(result.external_prandtl, 4.6288816615621506, rel_tol=1e-12))

        catalog = {compound.id: compound for compound in load_compounds(Path(__file__).parents[1] / "data/compounds/v1.json")}
        compounds = (catalog["N-pentane"], catalog["Ethane"])
        profile = pipe_estimated_htc_water_pr_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            pipe["HydraulicSegment,1,Results,2,InitialPressure"], external_temperature,
            compounds, (0.952380952380952, 0.0476190476190476),
            load_pr_interactions(Path(__file__).parents[1] / "data/interactions/pr-v1.json"),
            load_correlations(Path(__file__).parents[1] / "data/correlations/ideal-v1.json"),
            feed["PROP_MS_3"] / 1_000.0,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},InitialPressure"] for index in range(3, 7)),
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            external_water_velocity, external_water_heat_capacity,
            external_water_conductivity, external_water_viscosity,
            external_water_density,
        )

        self.assertEqual(len(profile.segment_results), 4)
        self.assertEqual(profile.segment_absorbed_radiation_w, (0.0,) * 4)
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=4e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_estimated_htc_dry_soil_matches_captured_dwsim_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-estimated-htc-dry-soil-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-estimated-htc-dry-soil-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        burial_depth_m = float(thermal.findtext("Velocidade"))
        soil_conductivity = 0.5

        self.assertEqual(thermal.findtext("TipoPerfil"), "Estimar_CGTC")
        self.assertEqual(thermal.findtext("Meio"), "4")
        self.assertEqual(external_temperature, 310.0)
        self.assertEqual(burial_depth_m, 1.0)
        self.assertEqual(thermal.findtext("Incluir_cti"), "true")
        self.assertEqual(thermal.findtext("Incluir_paredes"), "true")
        self.assertEqual(thermal.findtext("Incluir_cte"), "true")
        self.assertEqual(thermal.findtext("Incluir_isolamento"), "false")
        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "false")
        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            result = pipe_estimated_htc_soil(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], burial_depth_m, soil_conductivity,
            )
            self.assertTrue(math.isclose(result.internal_htc_w_m2_k, pipe[prefix + "HTCinternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.wall_htc_w_m2_k, pipe[prefix + "HTCpipewall"], rel_tol=5e-6))
            self.assertTrue(math.isclose(result.external_htc_w_m2_k, pipe[prefix + "HTCexternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.overall_htc_w_m2_k, pipe[prefix + "HTCoverall"], rel_tol=1e-7))
            self.assertEqual((result.external_reynolds, result.external_prandtl), (0.0, 0.0))

        catalog = {compound.id: compound for compound in load_compounds(Path(__file__).parents[1] / "data/compounds/v1.json")}
        compounds = (catalog["N-pentane"], catalog["Ethane"])
        profile = pipe_estimated_htc_soil_pr_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            pipe["HydraulicSegment,1,Results,2,InitialPressure"], external_temperature,
            compounds, (0.952380952380952, 0.0476190476190476),
            load_pr_interactions(Path(__file__).parents[1] / "data/interactions/pr-v1.json"),
            load_correlations(Path(__file__).parents[1] / "data/correlations/ideal-v1.json"),
            feed["PROP_MS_3"] / 1_000.0,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},InitialPressure"] for index in range(3, 7)),
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            burial_depth_m, soil_conductivity,
        )

        self.assertEqual(len(profile.segment_results), 4)
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_estimated_htc_moist_soil_matches_captured_dwsim_case(self):
        self.assertEqual(
            {name: dwsim_terrain_thermal_conductivity(name) for name in ("gravel", "stones", "dry_soil", "moist_soil")},
            {"gravel": 1.1, "stones": 1.95, "dry_soil": 0.5, "moist_soil": 2.2},
        )
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-estimated-htc-moist-soil-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-estimated-htc-moist-soil-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        burial_depth_m = float(thermal.findtext("Velocidade"))
        soil_conductivity = dwsim_terrain_thermal_conductivity("moist_soil")

        self.assertEqual(thermal.findtext("TipoPerfil"), "Estimar_CGTC")
        self.assertEqual(thermal.findtext("Meio"), "5")
        self.assertEqual(burial_depth_m, 1.0)
        self.assertEqual(thermal.findtext("Incluir_isolamento"), "false")
        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "false")
        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            result = pipe_estimated_htc_soil(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], burial_depth_m, soil_conductivity,
            )
            self.assertTrue(math.isclose(result.internal_htc_w_m2_k, pipe[prefix + "HTCinternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.wall_htc_w_m2_k, pipe[prefix + "HTCpipewall"], rel_tol=5e-6))
            self.assertTrue(math.isclose(result.external_htc_w_m2_k, pipe[prefix + "HTCexternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.overall_htc_w_m2_k, pipe[prefix + "HTCoverall"], rel_tol=1e-7))

        catalog = {compound.id: compound for compound in load_compounds(Path(__file__).parents[1] / "data/compounds/v1.json")}
        compounds = (catalog["N-pentane"], catalog["Ethane"])
        profile = pipe_estimated_htc_soil_pr_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"],
            pipe["HydraulicSegment,1,Results,2,InitialPressure"], external_temperature,
            compounds, (0.952380952380952, 0.0476190476190476),
            load_pr_interactions(Path(__file__).parents[1] / "data/interactions/pr-v1.json"),
            load_correlations(Path(__file__).parents[1] / "data/correlations/ideal-v1.json"),
            feed["PROP_MS_3"] / 1_000.0,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},InitialPressure"] for index in range(3, 7)),
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            burial_depth_m, soil_conductivity,
        )

        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_estimated_htc_insulation_matches_captured_dwsim_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-estimated-htc-insulated-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-estimated-htc-insulated-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        external_air_velocity = float(thermal.findtext("Velocidade"))
        insulation_thickness_m = float(thermal.findtext("Espessura"))
        insulation_conductivity = float(thermal.findtext("Condtermica"))

        self.assertEqual(thermal.findtext("TipoPerfil"), "Estimar_CGTC")
        self.assertEqual(thermal.findtext("Incluir_isolamento"), "true")
        self.assertEqual(thermal.findtext("Material"), "4")
        self.assertEqual(insulation_thickness_m, 0.025)
        self.assertEqual(insulation_conductivity, 0.035)
        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "false")
        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            result = pipe_estimated_htc_air(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                external_temperature, inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], external_air_velocity,
                insulation_thickness_m, insulation_conductivity,
            )
            self.assertTrue(math.isclose(result.internal_htc_w_m2_k, pipe[prefix + "HTCinternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.wall_htc_w_m2_k, pipe[prefix + "HTCpipewall"], rel_tol=5e-6))
            self.assertTrue(math.isclose(result.insulation_htc_w_m2_k, pipe[prefix + "HTCinsulation"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.external_htc_w_m2_k, pipe[prefix + "HTCexternal"], rel_tol=1e-12))
            self.assertTrue(math.isclose(result.overall_htc_w_m2_k, pipe[prefix + "HTCoverall"], rel_tol=1e-7))

        profile = pipe_estimated_htc_air_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"], external_temperature,
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            feed["PROP_MS_2"],
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            external_air_velocity, insulation_thickness_m, insulation_conductivity,
        )

        self.assertEqual(len(profile.segment_results), 4)
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_pipe_irradiation_matches_captured_dwsim_case(self):
        golden = json.loads((Path(__file__).parents[1] / "tests/golden/u3-pipe-thermal-estimated-htc-insulated-irradiated-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        with ZipFile(Path(__file__).parents[1] / "tests/u3-pipe-thermal-estimated-htc-insulated-irradiated-liquid-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        external_air_velocity = float(thermal.findtext("Velocidade"))
        insulation_thickness_m = float(thermal.findtext("Espessura"))
        insulation_conductivity = float(thermal.findtext("Condtermica"))
        solar_irradiation = float(thermal.findtext("SolarRadiationValue_kWh_m2"))
        absorption_efficiency = float(thermal.findtext("SolarRadiationAbsorptionEfficiency"))

        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "true")
        self.assertEqual(thermal.findtext("UseGlobalSolarRadiation"), "false")
        self.assertEqual(solar_irradiation, 0.01)
        self.assertEqual(absorption_efficiency, 0.1)
        absorbed_radiation = pipe_absorbed_solar_radiation(
            solar_irradiation, absorption_efficiency, outer_diameter_m,
            segment_length_m, feed["PROP_MS_4"],
        )
        self.assertTrue(math.isclose(absorbed_radiation, 2789.8341649201708, rel_tol=1e-12))

        for index in range(1, increments + 1):
            prefix = f"HydraulicSegment,1,Results,{index},"
            next_prefix = f"HydraulicSegment,1,Results,{index + 1},"
            htc = pipe_estimated_htc_air(
                (pipe[prefix + "InitialTemperature"] + pipe[next_prefix + "InitialTemperature"]) / 2.0,
                external_temperature, inner_diameter_m, outer_diameter_m, 4.5e-5,
                pipe[prefix + "VelocityLiquid"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                pipe[prefix + "ThermalConductivityLiquid"], pipe[prefix + "ViscosityLiquid"],
                pipe[prefix + "DensityLiquid"], external_air_velocity,
                insulation_thickness_m, insulation_conductivity,
            )
            result = pipe_irradiated_heat_transfer(
                pipe[prefix + "InitialTemperature"], external_temperature,
                htc.overall_htc_w_m2_k, outer_diameter_m, segment_length_m,
                feed["PROP_MS_2"], pipe[prefix + "HeatCapacityLiquid"] * 1_000.0,
                absorbed_radiation,
            )
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[prefix + "HeatTransfer"] * 1_000.0, rel_tol=3e-3))

        profile = pipe_estimated_htc_air_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"], external_temperature,
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            feed["PROP_MS_2"],
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            external_air_velocity, insulation_thickness_m, insulation_conductivity,
            solar_irradiation, absorption_efficiency, feed["PROP_MS_4"],
        )

        self.assertEqual(profile.segment_absorbed_radiation_w, (absorbed_radiation,) * 4)
        self.assertTrue(math.isclose(profile.absorbed_radiation_w, 4.0 * absorbed_radiation, rel_tol=1e-12))
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_pipe_global_irradiation_matches_captured_dwsim_case(self):
        root_path = Path(__file__).parents[1]
        golden = json.loads((root_path / "tests/golden/u3-pipe-thermal-estimated-htc-global-irradiated-liquid-pr-eos.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        pipe = {item["property"]: item["value"]["value"] for item in objects["PIPE-1"]["properties"]}
        feed = {item["property"]: item["value"]["value"] for item in objects["PIPE-FEED"]["properties"]}
        product = {item["property"]: item["value"]["value"] for item in objects["PIPE-PRODUCT"]["properties"]}
        energy = {item["property"]: item["value"]["value"] for item in objects["E1"]["properties"]}
        case_path = root_path / "tests/u3-pipe-thermal-estimated-htc-global-irradiated-liquid-pr-eos.dwxmz"
        with ZipFile(case_path) as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        source = next(item for item in root.findall("./SimulationObjects/SimulationObject") if item.findtext("CalculateHeatBalance") == "true")
        section = source.find("./Profile/Sections/Section")
        thermal = source.find("./ThermalProfile")
        weather = next(item for item in root.iter() if item.tag.split("}")[-1] == "CurrentWeather")
        increments = int(section.findtext("Incrementos"))
        inner_diameter_m = float(section.findtext("DI")) * 0.0254
        outer_diameter_m = float(section.findtext("DE")) * 0.0254
        segment_length_m = float(section.findtext("Comprimento")) / increments
        external_temperature = float(thermal.findtext("Temp_amb_estimar"))
        external_air_velocity = float(thermal.findtext("Velocidade"))
        insulation_thickness_m = float(thermal.findtext("Espessura"))
        insulation_conductivity = float(thermal.findtext("Condtermica"))
        local_irradiation = float(thermal.findtext("SolarRadiationValue_kWh_m2"))
        global_irradiation = float(weather.findtext("SolarIrradiation_kWh_m2"))
        absorption_efficiency = float(thermal.findtext("SolarRadiationAbsorptionEfficiency"))

        self.assertEqual(thermal.findtext("IncludeSolarRadiation"), "true")
        self.assertEqual(thermal.findtext("UseGlobalSolarRadiation"), "true")
        self.assertEqual(local_irradiation, 0.01)
        self.assertEqual(global_irradiation, 1.0)
        self.assertEqual(absorption_efficiency, 0.001)
        solar_irradiation = pipe_solar_irradiation_source(
            local_irradiation, True, global_irradiation,
        )
        self.assertEqual(solar_irradiation, global_irradiation)
        absorbed_radiation = pipe_absorbed_solar_radiation(
            solar_irradiation, absorption_efficiency, outer_diameter_m,
            segment_length_m, feed["PROP_MS_4"],
        )
        self.assertTrue(math.isclose(absorbed_radiation, 2789.8341649201708, rel_tol=1e-12))

        profile = pipe_estimated_htc_air_profile(
            pipe["HydraulicSegment,1,Results,2,InitialTemperature"], external_temperature,
            inner_diameter_m, outer_diameter_m, 4.5e-5, (segment_length_m,) * 4,
            feed["PROP_MS_2"],
            tuple(pipe[f"HydraulicSegment,1,Results,{index},VelocityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},HeatCapacityLiquid"] * 1_000.0 for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ThermalConductivityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},ViscosityLiquid"] for index in range(2, 6)),
            tuple(pipe[f"HydraulicSegment,1,Results,{index},DensityLiquid"] for index in range(2, 6)),
            external_air_velocity, insulation_thickness_m, insulation_conductivity,
            solar_irradiation, absorption_efficiency, feed["PROP_MS_4"],
        )

        self.assertEqual(profile.segment_absorbed_radiation_w, (absorbed_radiation,) * 4)
        for index, result in enumerate(profile.segment_results, 2):
            self.assertTrue(math.isclose(result.heat_transfer_w, pipe[f"HydraulicSegment,1,Results,{index},HeatTransfer"] * 1_000.0, rel_tol=3e-3))
        self.assertTrue(math.isclose(profile.outlet_temperature_k, product["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(profile.heat_transfer_w, -energy["PROP_ES_0"] * 1_000.0, rel_tol=2e-3, abs_tol=35.0))

    def test_estimated_htc_air_rejects_invalid_inputs(self):
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_air(300.0, 350.0, 0.1, 0.09, 4.5e-5, 2.0, 2_000.0, 0.1, 0.001, 700.0, 2.0)
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_air(300.0, 350.0, 0.1, 0.11, 4.5e-5, 2.0, 2_000.0, 0.1, 0.001, 700.0, 0.0)
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_air_profile(
                300.0, 350.0, 0.1, 0.11, 4.5e-5, (20.0,), 10.0,
                (2.0,), (2_000.0, 2_000.0), (0.1,), (0.001,), (700.0,), 2.0,
            )
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_air(300.0, 350.0, 0.1, 0.11, 4.5e-5, 2.0, 2_000.0, 0.1, 0.001, 700.0, 2.0, 0.025)
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_air_pr_gradient_profile(
                inlet_temperature_k=300.0, inlet_pressure_pa=500_000.0,
                base_external_temperature_k=350.0,
                external_temperature_gradient_k_m=math.nan, start_distance_m=20.0,
                compounds=(), composition=(), interactions={}, correlations=(),
                molar_flow_kmol_s=0.2, outlet_pressures_pa=(490_000.0,),
                inner_diameter_m=0.1, outer_diameter_m=0.11, roughness_m=4.5e-5,
                segment_lengths_m=(20.0,), liquid_velocities_m_s=(2.0,),
                heat_capacities_j_kg_k=(2_000.0,),
                thermal_conductivities_w_m_k=(0.1,), viscosities_pa_s=(0.001,),
                densities_kg_m3=(700.0,), external_air_velocity_m_s=2.0,
            )
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_water(
                300.0, 0.1, 0.11, 4.5e-5, 2.0, 2_000.0, 0.1, 0.001,
                700.0, 2.0, 4_000.0, 0.6, 0.0, 1_000.0,
            )
        with self.assertRaises(ValidationError):
            pipe_estimated_htc_soil(
                300.0, 0.1, 0.11, 4.5e-5, 2.0, 2_000.0, 0.1, 0.001,
                700.0, 0.05, 0.5,
            )
        with self.assertRaises(ValidationError):
            dwsim_terrain_thermal_conductivity("mud")
        with self.assertRaises(ValidationError):
            pipe_absorbed_solar_radiation(1.0, 1.1, 0.11, 20.0, 0.02)
        with self.assertRaises(ValidationError):
            pipe_solar_irradiation_source(0.01, "true", 1.0)
        with self.assertRaises(ValidationError):
            pipe_irradiated_heat_transfer(300.0, 350.0, 1.0, 0.11, 20.0, 10.0, 2_000.0, 1_000_000.0)

    def test_defined_htc_pipe_rejects_invalid_inputs(self):
        with self.assertRaises(ValidationError):
            pipe_defined_htc_heat_transfer(300.0, 350.0, 25.0, 0.0, 20.0, 10.0, 2_000.0)
        with self.assertRaises(ValidationError):
            pipe_defined_htc_heat_transfer(300.0, 350.0, -1.0, 0.1, 20.0, 10.0, 2_000.0)
        with self.assertRaises(ValidationError):
            pipe_defined_htc_profile(300.0, 350.0, 25.0, 0.1, (), 10.0, ())
        with self.assertRaises(ValidationError):
            pipe_defined_htc_profile(300.0, 350.0, 25.0, 0.1, (20.0,), 10.0, (2_000.0, 2_000.0))
        with self.assertRaises(ValidationError):
            pipe_defined_htc_gradient_profile(300.0, 350.0, math.nan, 25.0, 0.1, (20.0,), 10.0, (2_000.0,))
        with self.assertRaises(ValidationError):
            pipe_defined_htc_gradient_profile(300.0, 350.0, 0.1, 25.0, 0.1, (20.0,), 10.0, (2_000.0,), -1.0)
        with self.assertRaises(ValidationError):
            pipe_defined_heat_pr_profile(
                300.0, 500_000.0, (), (), {}, (), 0.2, (480_000.0,),
                0.11, (), 10_000.0, 5, (250.0, 350.0),
            )
        with self.assertRaises(ValidationError):
            liquid_pipe_supplied_state_profile(
                500_000.0, 300.0, 0.1, 0.11, 4.5e-5, (20.0,), (2.0,),
                (0.01,), (1_000.0,), (0.001,), 350.0, 25.0, (20.0,), 9.0, (2_000.0,),
            )

    def test_single_phase_pipe_rejects_invalid_geometry_and_properties(self):
        with self.assertRaises(ValidationError):
            pipe_pressure_drop(0.0, 100.0, 0.0, 0.0, 0.01, 1000.0, 0.001)
        with self.assertRaises(ValidationError):
            pipe_pressure_drop(0.1, 100.0, 101.0, 0.0, 0.01, 1000.0, 0.001)
        with self.assertRaises(ValidationError):
            pipe_pressure_drop(0.1, 100.0, 0.0, 0.0, 0.01, 1000.0, 0.0)
        with self.assertRaises(ValidationError):
            pipe_pressure_drop_profile(500_000.0, 0.1, 4.5e-5, (), (), (), (), ())
        with self.assertRaises(ValidationError):
            pipe_pressure_drop_profile(500_000.0, 0.1, 4.5e-5, (20.0,), (2.0,), (0.01,), (700.0,), (0.001, 0.001))
        with self.assertRaises(ValidationError):
            pipe_pressure_drop_profile(100.0, 0.1, 4.5e-5, (20.0,), (2.0,), (0.01,), (700.0,), (0.001,))
