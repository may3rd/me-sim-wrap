import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.streams import StreamState, flash_stream
from mesim.thermo.ideal import load_correlations
from mesim.unitops.basic import cooler, equilibrium_separator, heater, mix_streams, split_stream, valve
from mesim.unitops.pressure import pump


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

    def test_pump_raises_liquid_pressure_at_constant_efficiency(self):
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
            StreamState(300.0, 500_000.0, 2.0, compound_ids, (0.95, 0.05)),
            compounds, self.interactions, correlations,
        )
        result = pump(inlet, compounds, self.interactions, correlations, 1_000_000.0, 0.8, (290.0, 330.0))
        molar_volume = (
            inlet.flash.liquid_state.compressibility * 8_314.46261815324 * inlet.stream.temperature_k
            / inlet.stream.pressure_pa
        )

        self.assertEqual(result.outlet.stream.pressure_pa, 1_000_000.0)
        self.assertEqual(result.outlet.stream.molar_flow_kmol_s, inlet.stream.molar_flow_kmol_s)
        self.assertTrue(math.isclose(
            result.outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol,
            molar_volume * 500_000.0 / 0.8,
            abs_tol=1.0,
        ))
        self.assertTrue(math.isclose(
            result.energy.duty_w,
            inlet.stream.molar_flow_kmol_s * (result.outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol),
            rel_tol=1e-12,
        ))

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
