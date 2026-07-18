import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.errors import OutOfRangeError, ValidationError
from mesim.thermo.ideal import load_correlations
from mesim.thermo.raoult import (
    raoult_bubble_pressure,
    raoult_dew_pressure,
    raoult_fugacity_coefficients,
    raoult_tp_flash,
)
from mesim.thermo.systems import IDEAL_RAOULT, create_thermo_system


# DWSIM's captured Nested Loops PT flash uses a 1e-6 external-loop criterion.
# K-values are compared tightly; phase-split outputs use this reference tolerance.
DWSIM_REFERENCE_FLASH_PHASE_SPLIT_ABS_TOL = 3.0e-7


class IdealRaoultTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        catalog = {
            record.compound_id: record
            for record in load_correlations(ROOT / "data/correlations/ideal-v1.json")
        }
        cls.correlations = (catalog["Ethylene oxide"], catalog["Water"])

    def test_fugacity_coefficients_follow_dwsim_ideal_contract(self):
        temperature_k = 328.15
        pressure_pa = 101_325.0
        expected = tuple(
            record.vapor_pressure(temperature_k).value / pressure_pa
            for record in self.correlations
        )
        self.assertEqual(
            raoult_fugacity_coefficients(
                self.correlations, temperature_k, pressure_pa, "liquid"
            ),
            expected,
        )
        self.assertEqual(
            raoult_fugacity_coefficients(
                self.correlations, temperature_k, pressure_pa, "vapor"
            ),
            (1.0, 1.0),
        )

    def test_direct_bubble_and_dew_equations_close_material_balance(self):
        temperature_k = 328.15
        liquid = (0.35, 0.65)
        vapor_pressures = tuple(
            record.vapor_pressure(temperature_k).value for record in self.correlations
        )
        expected_bubble = math.fsum(
            fraction * pressure
            for fraction, pressure in zip(liquid, vapor_pressures)
        )
        bubble = raoult_bubble_pressure(
            self.correlations, liquid, temperature_k
        )
        self.assertEqual(bubble.pressure_pa, expected_bubble)
        self.assertTrue(bubble.report.converged)
        self.assertTrue(math.isclose(math.fsum(bubble.vapor_composition), 1.0))

        dew = raoult_dew_pressure(
            self.correlations, bubble.vapor_composition, temperature_k
        )
        self.assertTrue(math.isclose(dew.pressure_pa, expected_bubble, rel_tol=1.0e-14))
        for actual, expected in zip(dew.liquid_composition, liquid):
            self.assertTrue(math.isclose(actual, expected, rel_tol=1.0e-14))

    def test_tp_flash_handles_liquid_two_phase_and_vapor_regions(self):
        temperature_k = 328.15
        composition = (0.35, 0.65)
        bubble = raoult_bubble_pressure(
            self.correlations, composition, temperature_k
        ).pressure_pa
        dew = raoult_dew_pressure(
            self.correlations, composition, temperature_k
        ).pressure_pa
        self.assertGreater(bubble, dew)

        liquid = raoult_tp_flash(
            self.correlations, composition, temperature_k, bubble * 1.1
        )
        vapor = raoult_tp_flash(
            self.correlations, composition, temperature_k, dew * 0.9
        )
        two_phase = raoult_tp_flash(
            self.correlations, composition, temperature_k, (bubble + dew) / 2.0
        )
        self.assertEqual((liquid.phase, liquid.vapor_fraction), ("liquid", 0.0))
        self.assertEqual((vapor.phase, vapor.vapor_fraction), ("vapor", 1.0))
        self.assertEqual(two_phase.phase, "two-phase")
        self.assertTrue(two_phase.report.converged)
        self.assertTrue(0.0 < two_phase.vapor_fraction < 1.0)

    def test_domain_and_temperature_range_are_explicit(self):
        with self.assertRaises(ValidationError):
            raoult_bubble_pressure(self.correlations, (0.2, 0.2), 328.15)
        with self.assertRaises(ValidationError):
            raoult_fugacity_coefficients(
                self.correlations, 328.15, 101_325.0, "solid"
            )
        with self.assertRaises(OutOfRangeError):
            raoult_tp_flash(
                self.correlations, (0.5, 0.5), 150.0, 101_325.0
            )
        extrapolated = raoult_tp_flash(
            self.correlations,
            (0.5, 0.5),
            150.0,
            101_325.0,
            allow_extrapolation=True,
        )
        self.assertTrue(extrapolated.report.warnings)

    def test_two_phase_flash_matches_repeatable_dwsim_golden(self):
        path = ROOT / "tests/golden/ideal-raoult-ethylene-oxide-water-tp.json"
        repeat_path = ROOT / "tests/golden/ideal-raoult-ethylene-oxide-water-tp-repeat.json"
        golden = json.loads(path.read_text(encoding="utf-8-sig"))
        repeat = json.loads(repeat_path.read_text(encoding="utf-8-sig"))
        self.assertEqual(golden, repeat)
        self.assertEqual(golden["schema_version"], "dwsim-thermo-package-golden-1")
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.RaoultPropertyPackage",
        )

        catalog = {
            record.compound_id: record
            for record in load_correlations(ROOT / "data/correlations/ideal-v1.json")
        }
        inputs = golden["inputs"]
        outputs = golden["outputs"]
        correlations = tuple(catalog[name] for name in inputs["compounds"])
        system = create_thermo_system(IDEAL_RAOULT, correlations=correlations)
        result = system.tp_flash(
            tuple(inputs["composition"]),
            inputs["temperature_k"],
            inputs["pressure_pa"],
        )

        self.assertEqual(
            system.fugacity_coefficients(
                inputs["temperature_k"], inputs["pressure_pa"], "vapor"
            ),
            tuple(outputs["vapor_fugacity_coefficients"]),
        )
        for actual, expected in zip(
            result.equilibrium_ratios, outputs["liquid_fugacity_coefficients"]
        ):
            self.assertTrue(math.isclose(actual, expected, rel_tol=1.0e-14))
        self.assertEqual(result.phase, "two-phase")
        self.assertTrue(math.isclose(
            result.vapor_fraction,
            outputs["vapor_fraction"],
            rel_tol=0.0,
            abs_tol=DWSIM_REFERENCE_FLASH_PHASE_SPLIT_ABS_TOL,
        ))
        self.assertTrue(math.isclose(
            1.0 - result.vapor_fraction,
            outputs["liquid_fraction"],
            rel_tol=0.0,
            abs_tol=DWSIM_REFERENCE_FLASH_PHASE_SPLIT_ABS_TOL,
        ))
        for calculated, captured in (
            (result.liquid_composition, outputs["liquid_composition"]),
            (result.vapor_composition, outputs["vapor_composition"]),
        ):
            for actual, expected in zip(calculated, captured):
                self.assertTrue(math.isclose(
                    actual,
                    expected,
                    rel_tol=0.0,
                    abs_tol=DWSIM_REFERENCE_FLASH_PHASE_SPLIT_ABS_TOL,
                ))

    def test_dwsim_probe_is_deterministic_and_uses_explicit_optional_nulls(self):
        script = (ROOT / "scripts/capture_dwsim_thermo_package.ps1").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("captured_utc", script)
        self.assertIn(
            "new object[] { composition, pressure, temperature, propertyPackage, false, null }",
            script,
        )
        self.assertIn("[DwsimThermoProbe]::FullTypeName($propertyPackage)", script)


if __name__ == "__main__":
    unittest.main()
