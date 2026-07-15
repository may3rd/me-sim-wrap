import json
import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.streams import StreamState, flash_stream
from mesim.thermo.ideal import load_correlations
from mesim.unitops.basic import cooler, equilibrium_separator, heat_exchanger, heat_exchanger_efficiency, heat_exchanger_ua, heater, mix_streams, split_stream, valve
from mesim.unitops.pressure import compressor, expander, pump
from mesim.unitops.separation import component_separator


ROOT = Path(__file__).parents[1]


class BasicUnitOperationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.correlations = load_correlations(ROOT / "data/correlations/ideal-v1.json")
        cls.first = flash_stream(
            StreamState(180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (0.7, 0.3)),
            cls.compounds, cls.interactions, cls.correlations,
        )
        cls.second = flash_stream(
            StreamState(190.0, 600_000.0, 3.0, ("Methane", "Ethane"), (0.2, 0.8)),
            cls.compounds, cls.interactions, cls.correlations,
        )

    def test_mixer_closes_component_and_enthalpy_flows(self):
        outlet = mix_streams(
            (self.first, self.second), self.compounds, self.interactions, self.correlations,
            500_000.0, (140.0, 240.0),
        )

        self.assertTrue(outlet.flash.report.converged)
        self.assertEqual(outlet.stream.molar_flow_kmol_s, 5.0)
        self.assertTrue(math.isclose(outlet.stream.composition[0], 0.4, abs_tol=1e-12))
        self.assertTrue(math.isclose(
            outlet.enthalpy_j_per_kmol * outlet.stream.molar_flow_kmol_s,
            self.first.enthalpy_j_per_kmol * self.first.stream.molar_flow_kmol_s
            + self.second.enthalpy_j_per_kmol * self.second.stream.molar_flow_kmol_s,
            rel_tol=1e-6,
        ))

    def test_mixer_rejects_invalid_pressure_or_compound_order(self):
        with self.assertRaises(ValidationError):
            mix_streams((self.first, self.second), self.compounds, self.interactions, self.correlations, 600_001.0, (140.0, 240.0))
        zero = flash_stream(
            StreamState(180.0, 500_000.0, 0.0, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        with self.assertRaises(ValidationError):
            mix_streams((zero,), self.compounds, self.interactions, self.correlations, 500_000.0, (140.0, 240.0))
        reversed_stream = flash_stream(
            StreamState(180.0, 500_000.0, 1.0, ("Ethane", "Methane"), (0.3, 0.7)),
            tuple(reversed(self.compounds)), self.interactions, tuple(reversed(self.correlations)),
        )
        with self.assertRaises(ValidationError):
            mix_streams((self.first, reversed_stream), self.compounds, self.interactions, self.correlations, 500_000.0, (140.0, 240.0))

    def test_splitter_preserves_state_and_closes_flow(self):
        outlets = split_stream(self.first.stream, (0.25, 0.75))

        self.assertEqual(tuple(outlet.molar_flow_kmol_s for outlet in outlets), (0.5, 1.5))
        self.assertEqual(tuple(outlet.composition for outlet in outlets), ((0.7, 0.3), (0.7, 0.3)))
        self.assertEqual(tuple(outlet.temperature_k for outlet in outlets), (180.0, 180.0))
        self.assertEqual(math.fsum(outlet.molar_flow_kmol_s for outlet in outlets), self.first.stream.molar_flow_kmol_s)

    def test_splitter_rejects_invalid_fractions_and_preserves_zero_flow(self):
        with self.assertRaises(ValidationError):
            split_stream(self.first.stream, (0.2, 0.7))
        zero = StreamState(180.0, 500_000.0, 0.0, ("Methane", "Ethane"), (0.7, 0.3))
        self.assertEqual(tuple(outlet.molar_flow_kmol_s for outlet in split_stream(zero, (0.5, 0.5))), (0.0, 0.0))

    def test_heater_and_cooler_close_energy_at_fixed_pressure(self):
        heated = heater(self.first, self.compounds, self.interactions, self.correlations, 200.0)
        cooled = cooler(self.second, self.compounds, self.interactions, self.correlations, 170.0)

        self.assertEqual((heated.outlet.stream.temperature_k, cooled.outlet.stream.temperature_k), (200.0, 170.0))
        self.assertEqual((heated.outlet.stream.pressure_pa, cooled.outlet.stream.pressure_pa), (500_000.0, 600_000.0))
        self.assertGreater(heated.energy.duty_w, 0.0)
        self.assertLess(cooled.energy.duty_w, 0.0)
        self.assertTrue(math.isclose(
            heated.energy.duty_w,
            self.first.stream.molar_flow_kmol_s * (heated.outlet.enthalpy_j_per_kmol - self.first.enthalpy_j_per_kmol),
            rel_tol=1e-12,
        ))

    def test_thermal_operations_reject_wrong_sign_and_allow_zero_flow(self):
        with self.assertRaises(ValidationError):
            heater(self.first, self.compounds, self.interactions, self.correlations, 170.0)
        with self.assertRaises(ValidationError):
            cooler(self.first, self.compounds, self.interactions, self.correlations, 200.0)
        zero = flash_stream(
            StreamState(180.0, 500_000.0, 0.0, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        self.assertEqual(heater(zero, self.compounds, self.interactions, self.correlations, 200.0).energy.duty_w, 0.0)

    def test_valve_isenthalpically_flashes_at_lower_pressure(self):
        outlet = valve(
            self.second, self.compounds, self.interactions, self.correlations, 300_000.0, (140.0, 240.0),
        )

        self.assertTrue(outlet.flash.report.converged)
        self.assertEqual(outlet.stream.pressure_pa, 300_000.0)
        self.assertEqual(outlet.stream.molar_flow_kmol_s, self.second.stream.molar_flow_kmol_s)
        self.assertTrue(math.isclose(outlet.enthalpy_j_per_kmol, self.second.enthalpy_j_per_kmol, rel_tol=1e-6))
        with self.assertRaises(ValidationError):
            valve(self.second, self.compounds, self.interactions, self.correlations, 600_000.0, (140.0, 240.0))
        zero = flash_stream(
            StreamState(190.0, 600_000.0, 0.0, ("Methane", "Ethane"), (0.2, 0.8)),
            self.compounds, self.interactions, self.correlations,
        )
        self.assertEqual(
            valve(zero, self.compounds, self.interactions, self.correlations, 300_000.0, (140.0, 240.0)).stream.molar_flow_kmol_s,
            0.0,
        )

    def test_pump_matches_dwsim_eos_density_case(self):
        golden = json.loads((ROOT / "tests/golden/u1-pump-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        compound_ids = ("N-pentane", "Ethane")
        compounds = tuple(
            {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}[compound_id]
            for compound_id in compound_ids
        )
        correlations = tuple(
            {correlation.compound_id: correlation for correlation in load_correlations(ROOT / "data/correlations/ideal-v1.json")}[compound_id]
            for compound_id in compound_ids
        )
        inlet = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, compound_ids, (0.95, 0.05)),
            compounds, self.interactions, correlations,
        )
        result = pump(inlet, compounds, self.interactions, correlations, 1_000_000.0, 0.75, (290.0, 330.0))
        molar_volume = (
            inlet.flash.liquid_state.compressibility * 8_314.46261815324 * inlet.stream.temperature_k
            / inlet.stream.pressure_pa
        )

        self.assertEqual(result.outlet.stream.pressure_pa, 1_000_000.0)
        self.assertEqual(result.outlet.stream.molar_flow_kmol_s, inlet.stream.molar_flow_kmol_s)
        self.assertTrue(math.isclose(
            result.outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol,
            molar_volume * 500_000.0 / 0.75,
            abs_tol=1.0,
        ))
        self.assertTrue(math.isclose(
            result.energy.duty_w,
            inlet.stream.molar_flow_kmol_s * (result.outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol),
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(result.outlet.stream.temperature_k, records["product"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.outlet.enthalpy_j_per_kmol, records["product"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.energy.duty_w, records["PUMP"]["PROP_PU_3"] * 1_000.0, rel_tol=1e-4))

    def test_pump_rejects_nonliquid_inlet_efficiency_and_pressure_drop(self):
        with self.assertRaises(ValidationError):
            pump(self.first, self.compounds, self.interactions, self.correlations, 600_000.0, 0.8, (140.0, 240.0))
        liquid_ids = ("N-pentane", "Ethane")
        compounds = tuple(
            {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}[compound_id]
            for compound_id in liquid_ids
        )
        correlations = tuple(
            {correlation.compound_id: correlation for correlation in load_correlations(ROOT / "data/correlations/ideal-v1.json")}[compound_id]
            for compound_id in liquid_ids
        )
        liquid = flash_stream(StreamState(300.0, 500_000.0, 1.0, liquid_ids, (0.95, 0.05)), compounds, self.interactions, correlations)
        for outlet_pressure_pa, efficiency in ((500_000.0, 0.8), (1_000_000.0, 0.0), (1_000_000.0, 1.1), (1_000_000.0, True)):
            with self.assertRaises(ValidationError):
                pump(liquid, compounds, self.interactions, correlations, outlet_pressure_pa, efficiency, (290.0, 330.0))

    def test_compressor_matches_dwsim_eos_reference(self):
        golden = json.loads((ROOT / "tests/golden/u1-compressor-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        inlet = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        result = compressor(inlet, self.compounds, self.interactions, self.correlations, 1_000_000.0, 0.75, (300.0, 400.0))

        self.assertTrue(math.isclose(result.outlet.stream.temperature_k, records["3"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.outlet.enthalpy_j_per_kmol, records["3"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.energy.duty_w, records["compressor"]["PROP_CO_3"] * 1_000.0, rel_tol=1e-4))

    def test_compressor_rejects_nonvapor_inlet_efficiency_and_pressure_drop(self):
        with self.assertRaises(ValidationError):
            compressor(self.first, self.compounds, self.interactions, self.correlations, 600_000.0, 0.75, (300.0, 400.0))
        vapor = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        for outlet_pressure_pa, efficiency in ((500_000.0, 0.75), (1_000_000.0, 0.0), (1_000_000.0, 1.1), (1_000_000.0, True)):
            with self.assertRaises(ValidationError):
                compressor(vapor, self.compounds, self.interactions, self.correlations, outlet_pressure_pa, efficiency, (300.0, 400.0))

    def test_expander_matches_dwsim_eos_reference(self):
        golden = json.loads((ROOT / "tests/golden/u1-expander-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        inlet = flash_stream(
            StreamState(300.0, 1_000_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        result = expander(inlet, self.compounds, self.interactions, self.correlations, 500_000.0, 0.75, (240.0, 300.0))

        self.assertTrue(math.isclose(result.outlet.stream.temperature_k, records["3"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.outlet.enthalpy_j_per_kmol, records["3"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.energy.duty_w, -records["X-1"]["PROP_TU_3"] * 1_000.0, rel_tol=1e-4))

    def test_expander_rejects_nonvapor_inlet_efficiency_and_pressure_rise(self):
        with self.assertRaises(ValidationError):
            expander(self.first, self.compounds, self.interactions, self.correlations, 300_000.0, 0.75, (240.0, 300.0))
        vapor = flash_stream(
            StreamState(300.0, 1_000_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        for outlet_pressure_pa, efficiency in ((1_000_000.0, 0.75), (500_000.0, 0.0), (500_000.0, 1.1), (500_000.0, True)):
            with self.assertRaises(ValidationError):
                expander(vapor, self.compounds, self.interactions, self.correlations, outlet_pressure_pa, efficiency, (240.0, 300.0))

    def test_component_separator_matches_dwsim_mass_fraction_reference(self):
        golden = json.loads((ROOT / "tests/golden/u1-component-separator-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        compounds = ({compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}["N-pentane"], self.compounds[1])
        correlations = tuple(
            {correlation.compound_id: correlation for correlation in load_correlations(ROOT / "data/correlations/ideal-v1.json")}[compound.id]
            for compound in compounds
        )
        inlet = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("N-pentane", "Ethane"), (0.95, 0.05)),
            compounds, self.interactions, correlations,
        )
        result = component_separator(inlet, compounds, self.interactions, correlations, (0.10, 0.90))

        self.assertTrue(math.isclose(result.specified.stream.molar_flow_kmol_s * 1_000.0, records["3"]["PROP_MS_3"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.remainder.stream.molar_flow_kmol_s * 1_000.0, records["4"]["PROP_MS_3"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.specified.enthalpy_j_per_kmol, records["3"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.remainder.enthalpy_j_per_kmol, records["4"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.energy.duty_w, -records["CS-1"]["PROP_CP_0"] * 1_000.0, rel_tol=1e-4))

    def test_component_separator_rejects_invalid_fractions(self):
        inlet = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        for fractions in ((0.5,), (-0.1, 0.2), (1.1, 0.2), (math.nan, 0.2)):
            with self.assertRaises(ValidationError):
                component_separator(inlet, self.compounds, self.interactions, self.correlations, fractions)

    def test_heat_exchanger_matches_dwsim_fixed_duty_reference(self):
        golden = json.loads((ROOT / "tests/golden/u2-heat-duty-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        hot = flash_stream(
            StreamState(400.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        cold = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )

        result = heat_exchanger(hot, cold, self.compounds, self.interactions, self.correlations, 2_000.0, (280.0, 420.0))

        self.assertTrue(math.isclose(result.hot_outlet.stream.temperature_k, records["5"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.cold_outlet.stream.temperature_k, records["4"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.hot_outlet.enthalpy_j_per_kmol, records["5"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.cold_outlet.enthalpy_j_per_kmol, records["4"]["PROP_MS_9"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.hot_outlet.enthalpy_j_per_kmol * hot.stream.molar_flow_kmol_s + result.cold_outlet.enthalpy_j_per_kmol * cold.stream.molar_flow_kmol_s, hot.enthalpy_j_per_kmol * hot.stream.molar_flow_kmol_s + cold.enthalpy_j_per_kmol * cold.stream.molar_flow_kmol_s, abs_tol=1e-3))

    def test_heat_exchanger_ua_matches_dwsim_reference(self):
        golden = json.loads((ROOT / "tests/golden/u2-heat-ua-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        hot = flash_stream(
            StreamState(400.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        cold = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )

        result = heat_exchanger_ua(hot, cold, self.compounds, self.interactions, self.correlations, 25.0, 1.0, (280.0, 420.0))

        self.assertTrue(math.isclose(result.heat_duty_w, records["HX-1"]["PROP_HX_2"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.hot_outlet.stream.temperature_k, records["5"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.cold_outlet.stream.temperature_k, records["4"]["PROP_MS_0"], rel_tol=1e-4))

    def test_heat_exchanger_efficiency_matches_dwsim_reference(self):
        golden = json.loads((ROOT / "tests/golden/u2-heat-efficiency-pr-eos.json").read_text(encoding="utf-8-sig"))
        records = {
            object_["tag"]: {property_["property"]: property_["value"]["value"] for property_ in object_["properties"]}
            for object_ in golden["outputs"]["objects_after"]
        }
        hot = flash_stream(
            StreamState(400.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        cold = flash_stream(
            StreamState(300.0, 500_000.0, 0.002, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )

        result = heat_exchanger_efficiency(hot, cold, self.compounds, self.interactions, self.correlations, 50.0)

        self.assertTrue(math.isclose(result.heat_duty_w, records["HX-1"]["PROP_HX_2"] * 1_000.0, rel_tol=1e-4))
        self.assertTrue(math.isclose(result.hot_outlet.stream.temperature_k, records["5"]["PROP_MS_0"], rel_tol=1e-4))
        self.assertTrue(math.isclose(result.cold_outlet.stream.temperature_k, records["4"]["PROP_MS_0"], rel_tol=1e-4))

    def test_equilibrium_separator_routes_existing_flash_without_reflashing(self):
        result = equilibrium_separator(self.first, self.compounds, self.correlations)

        self.assertEqual(self.first.flash.phase, "two-phase")
        self.assertIsNotNone(result.liquid)
        self.assertIsNotNone(result.vapor)
        self.assertTrue(math.isclose(
            result.liquid.molar_flow_kmol_s + result.vapor.molar_flow_kmol_s,
            self.first.stream.molar_flow_kmol_s,
            abs_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            result.liquid.molar_flow_kmol_s * result.liquid.enthalpy_j_per_kmol
            + result.vapor.molar_flow_kmol_s * result.vapor.enthalpy_j_per_kmol,
            self.first.stream.molar_flow_kmol_s * self.first.enthalpy_j_per_kmol,
            rel_tol=1e-12,
        ))
        self.assertEqual(result.liquid.composition, self.first.flash.liquid_composition)
        self.assertEqual(result.vapor.composition, self.first.flash.vapor_composition)
        zero = flash_stream(
            StreamState(180.0, 500_000.0, 0.0, ("Methane", "Ethane"), (0.7, 0.3)),
            self.compounds, self.interactions, self.correlations,
        )
        zero_result = equilibrium_separator(zero, self.compounds, self.correlations)
        self.assertEqual((zero_result.liquid.molar_flow_kmol_s, zero_result.vapor.molar_flow_kmol_s), (0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
