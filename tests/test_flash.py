import math
import sys
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.thermo.flash import pr_stability, rachford_rice, tp_flash


ROOT = Path(__file__).parents[1]


class RachfordRiceTest(unittest.TestCase):
    """Vectors use Rachford and Rice, DOI 10.2118/952327-G, residual in Phase 6 Task 1."""

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
        with self.assertRaises(ValidationError):
            rachford_rice((0.5, 0.5), (2.0, "invalid"))
        with self.assertRaises(ValidationError):
            rachford_rice((0.5, 0.5), (2.0, 0.5), tolerance="invalid")

    def test_positive_finite_extreme_k_values_do_not_leak_arithmetic_errors(self):
        result = rachford_rice((0.5, 0.5), (1e308, 5e-324))

        self.assertTrue(result.report.converged)
        self.assertEqual(result.phase, "two-phase")
        self.assertTrue(math.isclose(result.vapor_fraction, 0.5, abs_tol=1e-12))


class PRStabilityTest(unittest.TestCase):
    """Vectors use Michelsen Part I, DOI 10.1016/0378-3812(82)85001-2.

    Pure-state limits are also traceable to tests/golden/pr-t1.json, DWSIM 9.0.4.
    """

    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")

    def test_dwsim_backed_pure_limits_are_stable_with_zero_fraction(self):
        vapor = pr_stability(self.compounds, (1.0, 0.0), self.interactions, 228.672, 459_900)
        liquid = pr_stability(self.compounds, (1.0, 0.0), self.interactions, 133.392, 2_299_500)

        for result in (vapor, liquid):
            self.assertTrue(result.report.converged)
            self.assertTrue(result.stable)
            self.assertTrue(result.vapor_like.report.converged)
            self.assertTrue(result.liquid_like.report.converged)
            self.assertGreaterEqual(min(result.vapor_like.tangent_plane_distance, result.liquid_like.tangent_plane_distance), -1e-10)

    def test_distinguishes_stable_vapor_stable_liquid_and_unstable_feed(self):
        vapor = pr_stability(self.compounds, (0.7, 0.3), self.interactions, 200.0, 500_000)
        liquid = pr_stability(self.compounds, (0.7, 0.3), self.interactions, 150.0, 1_000_000)
        unstable = pr_stability(self.compounds, (0.7, 0.3), self.interactions, 180.0, 500_000)

        self.assertTrue(vapor.stable)
        self.assertEqual(vapor.feed_phase, "vapor")
        self.assertTrue(liquid.stable)
        self.assertEqual(liquid.feed_phase, "liquid")
        self.assertFalse(unstable.stable)
        self.assertLess(unstable.liquid_like.tangent_plane_distance, -0.5)

    def test_near_critical_result_is_deterministic(self):
        first = pr_stability(self.compounds, (0.7, 0.3), self.interactions, 250.0, 5_000_000)
        second = pr_stability(self.compounds, (0.7, 0.3), self.interactions, 250.0, 5_000_000)

        self.assertEqual(first, second)
        self.assertTrue(first.report.converged)
        self.assertTrue(first.stable)
        self.assertLess(first.report.residual, 1e-10)

    def test_reports_trial_iteration_exhaustion(self):
        result = pr_stability(
            self.compounds,
            (0.7, 0.3),
            self.interactions,
            180.0,
            500_000,
            max_iterations=1,
            tolerance=1e-30,
        )

        self.assertFalse(result.report.converged)
        self.assertIsNotNone(result.report.failure_reason)


class PRTPFlashTest(unittest.TestCase):
    """Vectors use Michelsen Part II, DOI 10.1016/0378-3812(82)85002-4."""

    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")

    def test_returns_stable_liquid_vapor_and_near_critical_states(self):
        liquid = tp_flash(self.compounds, (0.7, 0.3), self.interactions, 150.0, 1_000_000)
        vapor = tp_flash(self.compounds, (0.7, 0.3), self.interactions, 200.0, 500_000)
        near_critical = tp_flash(self.compounds, (0.7, 0.3), self.interactions, 250.0, 5_000_000)

        self.assertEqual((liquid.phase, liquid.vapor_fraction), ("liquid", 0.0))
        self.assertIsNotNone(liquid.liquid_state)
        self.assertEqual((vapor.phase, vapor.vapor_fraction), ("vapor", 1.0))
        self.assertIsNotNone(vapor.vapor_state)
        self.assertEqual(near_critical.phase, "single")
        self.assertTrue(all(result.report.converged for result in (liquid, vapor, near_critical)))

    def test_two_phase_flash_closes_fugacity_and_material_balance(self):
        feed = (0.7, 0.3)
        result = tp_flash(self.compounds, feed, self.interactions, 180.0, 500_000)

        self.assertTrue(result.report.converged)
        self.assertEqual(result.phase, "two-phase")
        self.assertGreater(result.vapor_fraction, 0.0)
        self.assertLess(result.vapor_fraction, 1.0)
        self.assertIsNotNone(result.liquid_state)
        self.assertIsNotNone(result.vapor_state)
        for index, overall in enumerate(feed):
            liquid = result.liquid_composition[index]
            vapor = result.vapor_composition[index]
            self.assertTrue(math.isclose(overall, (1 - result.vapor_fraction) * liquid + result.vapor_fraction * vapor, rel_tol=1e-10))
            liquid_fugacity = liquid * result.liquid_state.fugacity_coefficients[index]
            vapor_fugacity = vapor * result.vapor_state.fugacity_coefficients[index]
            self.assertTrue(math.isclose(liquid_fugacity, vapor_fugacity, rel_tol=1e-8))

    def test_zero_fraction_and_repeat_execution_are_deterministic(self):
        first = tp_flash(self.compounds, (1.0, 0.0), self.interactions, 228.672, 459_900)
        second = tp_flash(self.compounds, (1.0, 0.0), self.interactions, 228.672, 459_900)

        self.assertEqual(first, second)
        self.assertTrue(first.report.converged)

    def test_rejects_invalid_inputs_and_reports_iteration_exhaustion(self):
        with self.assertRaises(ValidationError):
            tp_flash(self.compounds, (0.7, 0.3), self.interactions, 0.0, 500_000)
        with self.assertRaises(ValidationError):
            tp_flash(self.compounds, (0.7, 0.4), self.interactions, 180.0, 500_000)

        result = tp_flash(
            self.compounds,
            (0.7, 0.3),
            self.interactions,
            180.0,
            500_000,
            max_iterations=1,
        )
        self.assertFalse(result.report.converged)
        self.assertIsNotNone(result.report.failure_reason)


if __name__ == "__main__":
    unittest.main()
