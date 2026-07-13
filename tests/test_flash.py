import math
import sys
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.thermo.flash import rachford_rice


class RachfordRiceTest(unittest.TestCase):
    def test_classifies_single_phase_limits(self):
        liquid = rachford_rice((0.5, 0.5), (0.5, 0.8))
        self.assertEqual((liquid.phase, liquid.vapor_fraction), ("liquid", 0.0))
        self.assertEqual((liquid.liquid_composition, liquid.vapor_composition), ((0.5, 0.5), ()))
        self.assertTrue(liquid.report.converged)

        vapor = rachford_rice((0.5, 0.5), (2.0, 3.0))
        self.assertEqual((vapor.phase, vapor.vapor_fraction), ("vapor", 1.0))
        self.assertEqual((vapor.liquid_composition, vapor.vapor_composition), ((), (0.5, 0.5)))
        self.assertTrue(vapor.report.converged)

    def test_solves_known_two_phase_vector_and_preserves_balance(self):
        result = rachford_rice((0.5, 0.5), (2.0, 0.5))

        self.assertTrue(result.report.converged)
        self.assertEqual(result.phase, "two-phase")
        self.assertTrue(math.isclose(result.vapor_fraction, 0.5, abs_tol=1e-12))
        self.assertTrue(math.isclose(sum(result.liquid_composition), 1.0, abs_tol=1e-12))
        self.assertTrue(math.isclose(sum(result.vapor_composition), 1.0, abs_tol=1e-12))
        for feed, liquid, vapor in zip((0.5, 0.5), result.liquid_composition, result.vapor_composition):
            self.assertTrue(math.isclose(feed, 0.5 * liquid + 0.5 * vapor, abs_tol=1e-12))
        self.assertEqual(result.report.algorithm, "Rachford-Rice bisection")
        self.assertIsNone(result.report.failure_reason)
        with self.assertRaises(FrozenInstanceError):
            result.vapor_fraction = 0.0

    def test_reports_degeneracy_and_iteration_exhaustion(self):
        degenerate = rachford_rice((0.5, 0.5), (1.0, 1.0))
        self.assertFalse(degenerate.report.converged)
        self.assertEqual(degenerate.phase, "indeterminate")
        self.assertIn("indeterminate", degenerate.report.failure_reason)

        exhausted = rachford_rice((0.5, 0.5), (2.0, 0.5), max_iterations=1, tolerance=1e-30)
        self.assertFalse(exhausted.report.converged)
        self.assertEqual(exhausted.report.iterations, 1)
        self.assertIsNotNone(exhausted.report.failure_reason)

    def test_rejects_invalid_caller_input(self):
        cases = (
            ((1.0,), ()),
            ((0.5, 0.6), (2.0, 0.5)),
            ((-0.1, 1.1), (2.0, 0.5)),
            ((math.nan, math.nan), (2.0, 0.5)),
            ((0.5, 0.5), (0.0, 1.0)),
            ((0.5, 0.5), (math.inf, 1.0)),
        )
        for composition, k_values in cases:
            with self.subTest(composition=composition, k_values=k_values):
                with self.assertRaises(ValidationError):
                    rachford_rice(composition, k_values)
        with self.assertRaises(ValidationError):
            rachford_rice((0.5, 0.5), (2.0, 0.5), max_iterations=0)
        with self.assertRaises(ValidationError):
            rachford_rice((0.5, 0.5), (2.0, 0.5), tolerance=0.0)


if __name__ == "__main__":
    unittest.main()
