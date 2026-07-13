import math
import sys
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import ValidationError
from mesim.compounds import load_compounds, load_pr_interactions
from mesim.thermo.flash import (
    bubble_pressure,
    dew_pressure,
    flash_enthalpy,
    ideal_mixture_enthalpy,
    ph_flash,
    phase_enthalpy,
    pr_stability,
    rachford_rice,
    tp_flash,
)
from mesim.thermo.ideal import load_correlations


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


class PRBubbleDewPressureTest(unittest.TestCase):
    """Vectors use Michelsen Part II, DOI 10.1016/0378-3812(82)85002-4."""

    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")

    def test_pure_component_bubble_and_dew_pressures_agree(self):
        bubble = bubble_pressure(self.compounds, (1.0, 0.0), self.interactions, 150.0, (100_000.0, 5_000_000.0))
        dew = dew_pressure(self.compounds, (1.0, 0.0), self.interactions, 150.0, (100_000.0, 5_000_000.0))

        self.assertTrue(bubble.report.converged)
        self.assertTrue(dew.report.converged)
        self.assertTrue(math.isclose(bubble.pressure_pa, dew.pressure_pa, rel_tol=1e-8))

    def test_mixture_bubble_is_not_below_dew_pressure(self):
        bubble = bubble_pressure(self.compounds, (0.7, 0.3), self.interactions, 180.0, (100_000.0, 3_000_000.0))
        dew = dew_pressure(self.compounds, (0.7, 0.3), self.interactions, 180.0, (100_000.0, 500_000.0))

        self.assertTrue(bubble.report.converged)
        self.assertTrue(dew.report.converged)
        self.assertGreaterEqual(bubble.pressure_pa, dew.pressure_pa)
        self.assertTrue(math.isclose(math.fsum(bubble.liquid_composition), 1.0, abs_tol=1e-12))
        self.assertTrue(math.isclose(math.fsum(dew.vapor_composition), 1.0, abs_tol=1e-12))

    def test_requires_a_valid_bracket_and_reports_exhaustion(self):
        with self.assertRaises(ValidationError):
            bubble_pressure(self.compounds, (0.7, 0.3), self.interactions, 180.0, (0.0, 5_000_000.0))
        with self.assertRaises(ValidationError):
            dew_pressure(self.compounds, (0.7, 0.3), self.interactions, 180.0, (5_000_000.0, 100_000.0))
        with self.assertRaises(ValidationError):
            bubble_pressure(self.compounds, (0.7, 0.3), self.interactions, 180.0, (100_000.0, 5_000_000.0))

        result = bubble_pressure(
            self.compounds,
            (0.7, 0.3),
            self.interactions,
            180.0,
            (100_000.0, 5_000_000.0),
            max_iterations=1,
        )
        self.assertFalse(result.report.converged)
        self.assertIsNotNone(result.report.failure_reason)


class PRPHFlashTest(unittest.TestCase):
    """Caloric vectors use the Phase 6 reference state: 298.15 K and 101325 Pa."""

    @classmethod
    def setUpClass(cls):
        compounds = {compound.id: compound for compound in load_compounds(ROOT / "data/compounds/v1.json")}
        cls.compounds = (compounds["Methane"], compounds["Ethane"])
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.correlations = load_correlations(ROOT / "data/correlations/ideal-v1.json")

    def test_ideal_and_phase_enthalpy_use_the_explicit_reference(self):
        composition = (0.7, 0.3)
        ideal = ideal_mixture_enthalpy(self.compounds, composition, self.correlations, 298.15)
        flash = tp_flash(self.compounds, composition, self.interactions, 200.0, 500_000.0)
        phase = phase_enthalpy(self.compounds, composition, self.correlations, 200.0, flash.vapor_state)

        self.assertEqual(ideal, 0.0)
        self.assertTrue(math.isclose(phase, ideal_mixture_enthalpy(self.compounds, composition, self.correlations, 200.0) + flash.vapor_state.departure_enthalpy_j_per_kmol, rel_tol=1e-12))

    def test_total_two_phase_enthalpy_is_weighted_by_vapor_fraction(self):
        composition = (0.7, 0.3)
        flash = tp_flash(self.compounds, composition, self.interactions, 180.0, 500_000.0)
        total = flash_enthalpy(self.compounds, self.correlations, flash)
        liquid = phase_enthalpy(self.compounds, flash.liquid_composition, self.correlations, 180.0, flash.liquid_state)
        vapor = phase_enthalpy(self.compounds, flash.vapor_composition, self.correlations, 180.0, flash.vapor_state)

        self.assertEqual(flash.phase, "two-phase")
        self.assertTrue(math.isclose(total, (1.0 - flash.vapor_fraction) * liquid + flash.vapor_fraction * vapor, rel_tol=1e-12))

    def test_ph_flash_round_trips_single_and_phase_crossing_states(self):
        composition = (0.7, 0.3)
        for temperature_k in (150.0, 180.0, 200.0):
            with self.subTest(temperature_k=temperature_k):
                source = tp_flash(self.compounds, composition, self.interactions, temperature_k, 500_000.0)
                target = flash_enthalpy(self.compounds, self.correlations, source)
                result = ph_flash(self.compounds, composition, self.interactions, self.correlations, 500_000.0, target, (140.0, 210.0))

                self.assertTrue(result.report.converged)
                self.assertTrue(math.isclose(result.temperature_k, temperature_k, rel_tol=1e-7, abs_tol=1e-5))
                self.assertTrue(math.isclose(result.enthalpy_j_per_kmol, target, rel_tol=1e-6, abs_tol=1e-3))

    def test_ph_flash_rejects_invalid_brackets_and_reports_unreachable_target(self):
        with self.assertRaises(ValidationError):
            ph_flash(self.compounds, (0.7, 0.3), self.interactions, self.correlations, 500_000.0, 0.0, (0.0, 210.0))

        result = ph_flash(self.compounds, (0.7, 0.3), self.interactions, self.correlations, 500_000.0, 1e12, (140.0, 210.0))
        self.assertFalse(result.report.converged)
        self.assertIn("does not enclose", result.report.failure_reason)

    def test_ph_flash_is_deterministic_and_reports_iteration_exhaustion(self):
        composition = (0.7, 0.3)
        source = tp_flash(self.compounds, composition, self.interactions, 180.0, 500_000.0)
        target = flash_enthalpy(self.compounds, self.correlations, source)
        first = ph_flash(self.compounds, composition, self.interactions, self.correlations, 500_000.0, target, (140.0, 210.0))
        second = ph_flash(self.compounds, composition, self.interactions, self.correlations, 500_000.0, target, (140.0, 210.0))
        exhausted = ph_flash(self.compounds, composition, self.interactions, self.correlations, 500_000.0, target, (140.0, 210.0), max_iterations=1)

        self.assertEqual(first, second)
        self.assertFalse(exhausted.report.converged)
        self.assertIsNotNone(exhausted.report.failure_reason)


if __name__ == "__main__":
    unittest.main()
