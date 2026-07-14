import math
import sys
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.streams import EnergyStream, StreamState, flash_stream
from mesim.thermo.ideal import load_correlations


ROOT = Path(__file__).parents[1]


class StreamStateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.correlations = load_correlations(ROOT / "data/correlations/ideal-v1.json")

    def test_preserves_ordered_compound_ids_and_is_immutable(self):
        stream = StreamState(180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (0.7, 0.3))

        self.assertEqual(stream.compound_ids, ("Methane", "Ethane"))
        self.assertEqual(stream.composition, (0.7, 0.3))
        with self.assertRaises(FrozenInstanceError):
            stream.pressure_pa = 1.0

    def test_rejects_invalid_state_without_flashing(self):
        cases = (
            (0.0, 500_000.0, 2.0, ("Methane",), (1.0,)),
            (180.0, 0.0, 2.0, ("Methane",), (1.0,)),
            (180.0, 500_000.0, -1.0, ("Methane",), (1.0,)),
            (180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (1.0,)),
            (180.0, 500_000.0, 2.0, ("Methane", "Methane"), (0.7, 0.3)),
            (180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (0.7, 0.4)),
            (180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (math.nan, math.nan)),
        )
        for values in cases:
            with self.subTest(values=values):
                with self.assertRaises(ValidationError):
                    StreamState(*values)

    def test_flash_is_explicit_and_rejects_compound_order_mismatch(self):
        stream = StreamState(180.0, 500_000.0, 2.0, ("Methane", "Ethane"), (0.7, 0.3))
        result = flash_stream(stream, self.compounds, self.interactions, self.correlations)

        self.assertTrue(result.flash.report.converged)
        self.assertEqual(result.flash.phase, "two-phase")
        self.assertTrue(math.isfinite(result.enthalpy_j_per_kmol))
        with self.assertRaises(ValidationError):
            flash_stream(stream, tuple(reversed(self.compounds)), self.interactions, self.correlations)

    def test_energy_stream_is_immutable_signed_si_power(self):
        heating = EnergyStream(25_000.0)
        cooling = EnergyStream(-25_000.0)

        self.assertEqual((heating.duty_w, cooling.duty_w), (25_000.0, -25_000.0))
        with self.assertRaises(FrozenInstanceError):
            heating.duty_w = 0.0
        for invalid in (True, math.nan, math.inf, "100"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValidationError):
                    EnergyStream(invalid)


if __name__ == "__main__":
    unittest.main()
