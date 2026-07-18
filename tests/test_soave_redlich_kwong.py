import json
import math
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds, load_pr_interactions, load_srk_interactions
from mesim.errors import ValidationError
from mesim.thermo.peng_robinson import R
from mesim.thermo.soave_redlich_kwong import (
    SoaveRedlichKwong,
    SoaveRedlichKwongMixture,
)
from mesim.thermo.systems import SOAVE_REDLICH_KWONG, create_thermo_system


class SoaveRedlichKwongTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.compound_catalog = {
            compound.id: compound
            for compound in load_compounds(ROOT / "data/compounds/v1.json")
        }
        cls.compounds = (
            cls.compound_catalog["Methane"],
            cls.compound_catalog["Ethane"],
        )
        cls.interactions = load_srk_interactions(ROOT / "data/interactions/srk-v1.json")

    def test_parameters_follow_the_dwsim_srk_equation(self):
        compound = self.compounds[0]
        temperature_k = 180.0
        model = SoaveRedlichKwong(compound)
        parameters = model.parameters(temperature_k)
        omega = compound.acentric_factor.value
        kappa = 0.48 + 1.574 * omega - 0.176 * omega**2
        alpha = (
            1.0
            + kappa
            * (1.0 - math.sqrt(temperature_k / compound.critical_temperature.value))
        ) ** 2
        expected_a = (
            0.42748
            * alpha
            * R**2
            * compound.critical_temperature.value**2
            / compound.critical_pressure.value
        )
        expected_b = (
            0.08664
            * R
            * compound.critical_temperature.value
            / compound.critical_pressure.value
        )
        self.assertEqual(parameters.kappa, kappa)
        self.assertEqual(parameters.alpha, alpha)
        self.assertTrue(math.isclose(parameters.a_pa_m6_per_kmol2, expected_a, rel_tol=1.0e-15))
        self.assertEqual(parameters.b_m3_per_kmol, expected_b)

    def test_source_backed_interactions_are_strict_and_model_specific(self):
        self.assertEqual(self.interactions.model, "Soave-Redlich-Kwong")
        self.assertEqual(len(self.interactions.pairs), 52)
        self.assertEqual(self.interactions.get("Methane", "Ethane"), -0.0089)
        with self.assertRaises(ValidationError):
            self.interactions.get("Methane", "Acetone")
        with self.assertRaises(ValidationError):
            create_thermo_system(
                SOAVE_REDLICH_KWONG,
                compounds=self.compounds,
                interactions=load_pr_interactions(ROOT / "data/interactions/pr-v1.json"),
            )

    def test_phase_fugacities_match_repeatable_dwsim_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/srk-methane-ethane-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (ROOT / "tests/golden/srk-methane-ethane-state-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(golden["case_id"], "srk-methane-ethane-state")
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.SRKPropertyPackage",
        )
        inputs = golden["inputs"]
        outputs = golden["outputs"]
        system = create_thermo_system(
            SOAVE_REDLICH_KWONG,
            compounds=tuple(self.compound_catalog[name] for name in inputs["compounds"]),
            interactions=self.interactions,
        )
        for phase in ("liquid", "vapor"):
            state = system.state(
                tuple(inputs["composition"]),
                inputs["temperature_k"],
                inputs["pressure_pa"],
                phase,
            )
            for actual, expected in zip(
                state.fugacity_coefficients,
                outputs[f"{phase}_fugacity_coefficients"],
            ):
                self.assertTrue(math.isclose(actual, expected, rel_tol=5.0e-12))
            self.assertGreater(state.compressibility, 0.0)
            self.assertGreater(state.density_kg_per_m3, 0.0)

        flash = system.tp_flash(
            tuple(inputs["composition"]),
            inputs["temperature_k"],
            inputs["pressure_pa"],
        )
        self.assertTrue(flash.report.converged)
        self.assertEqual(flash.phase, "two-phase")
        # DWSIM stopped this reference after three external iterations at its
        # configured 1e-6 criterion; MeSim converges the fugacity residual tighter.
        self.assertTrue(math.isclose(
            flash.vapor_fraction,
            outputs["vapor_fraction"],
            rel_tol=0.0,
            abs_tol=1.0e-7,
        ))
        for actual, expected in zip(
            flash.liquid_composition, outputs["liquid_composition"]
        ):
            self.assertTrue(math.isclose(actual, expected, rel_tol=0.0, abs_tol=2.0e-6))
        for actual, expected in zip(
            flash.vapor_composition, outputs["vapor_composition"]
        ):
            self.assertTrue(math.isclose(actual, expected, rel_tol=0.0, abs_tol=1.0e-7))
        for actual, expected in zip(
            flash.equilibrium_ratios, outputs["equilibrium_ratios"]
        ):
            self.assertTrue(math.isclose(actual, expected, rel_tol=2.0e-5))

    def test_correlation_generator_has_no_srk_drift(self):
        completed = subprocess.run(
            [sys.executable, str(ROOT / "scripts/extract_chemsep_correlations.py")],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)

    def test_invalid_domains_fail_before_eos_arithmetic(self):
        with self.assertRaises(ValidationError):
            SoaveRedlichKwongMixture(
                self.compounds, (0.5, 0.4), self.interactions
            )
        with self.assertRaises(ValidationError):
            SoaveRedlichKwong(self.compounds[0]).state(180.0, 500_000.0, "solid")


if __name__ == "__main__":
    unittest.main()
