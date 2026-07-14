import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.streams import StreamState, flash_stream
from mesim.thermo.ideal import load_correlations
from mesim.unitops.basic import cooler, heater, mix_streams, split_stream


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


if __name__ == "__main__":
    unittest.main()
