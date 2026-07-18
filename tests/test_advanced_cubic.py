import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds, load_pr_interactions, load_srk_interactions
from mesim.thermo.advanced_cubic import load_advanced_cubic_data
from mesim.thermo.systems import (
    PENG_ROBINSON_1978_ADVANCED,
    PengRobinson1978AdvancedSystem,
    SOAVE_REDLICH_KWONG_ADVANCED,
    SoaveRedlichKwongAdvancedSystem,
    create_thermo_system,
)


class AdvancedCubicTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = {
            compound.id: compound
            for path in (
                ROOT / "data/compounds/v1.json",
                ROOT / "data/compounds/advanced-eos-v1.json",
            )
            for compound in load_compounds(path)
        }
        cls.data = load_advanced_cubic_data(
            ROOT / "data/interactions/prsrk-advanced-v1.json"
        )
        cls.interactions = load_pr_interactions(
            ROOT / "data/interactions/pr-v1.json"
        )
        cls.srk_interactions = load_srk_interactions(
            ROOT / "data/interactions/srk-v1.json"
        )

    def test_full_mercury_source_domain_and_polynomial_are_frozen(self):
        self.assertEqual(self.data.source_revision, "9.0.5.0")
        self.assertEqual(len(self.data.interactions), 13)
        self.assertIsNone(self.data.interaction("mercury", "Methane"))
        record = self.data.interaction("Mercury", "N-pentane")
        temperature = 400.0
        expected = 1.0036e-4 * temperature**2 - 5.5429e-2 * temperature + 7.5340
        self.assertTrue(
            math.isclose(record.value(temperature), expected, rel_tol=3.0e-15)
        )

    def test_scoped_mercury_eos_record_is_source_backed(self):
        mercury = self.catalog["Mercury"]
        self.assertEqual(mercury.cas, "7439-97-6")
        self.assertEqual(mercury.critical_temperature.value, 1735.0)
        self.assertEqual(mercury.critical_pressure.value, 160_803_000.0)
        self.assertEqual(mercury.acentric_factor.value, -0.16445)

    def test_advanced_pr78_matches_repeatable_dwsim_mercury_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/pr78-advanced-mercury-methane-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (
                ROOT
                / "tests/golden/pr78-advanced-mercury-methane-state-repeat.json"
            ).read_text(encoding="utf-8-sig")
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.AdvancedEOS.PengRobinson1978AdvancedPropertyPackage",
        )
        inputs = golden["inputs"]
        system = create_thermo_system(
            PENG_ROBINSON_1978_ADVANCED,
            compounds=tuple(self.catalog[name] for name in inputs["compounds"]),
            interactions=self.interactions,
            advanced_data=self.data,
        )
        self.assertIsInstance(system, PengRobinson1978AdvancedSystem)
        evaluated = self.data.evaluated(
            inputs["temperature_k"], inputs["pressure_pa"], self.interactions
        )
        self.assertEqual(evaluated.get("Mercury", "Methane"), 0.4)
        for phase in ("liquid", "vapor"):
            state = system.state(
                tuple(inputs["composition"]),
                inputs["temperature_k"],
                inputs["pressure_pa"],
                phase,
            )
            for actual, expected in zip(
                state.fugacity_coefficients,
                golden["outputs"][f"{phase}_fugacity_coefficients"],
            ):
                self.assertTrue(math.isclose(actual, expected, rel_tol=1.0e-14))

    def test_advanced_srk_matches_configured_repeatable_dwsim_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/srk-advanced-mercury-n-pentane-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (
                ROOT
                / "tests/golden/srk-advanced-mercury-n-pentane-state-repeat.json"
            ).read_text(encoding="utf-8-sig")
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.AdvancedEOS.SoaveRedlichKwongAdvancedPropertyPackage",
        )
        inputs = golden["inputs"]
        system = create_thermo_system(
            SOAVE_REDLICH_KWONG_ADVANCED,
            compounds=tuple(self.catalog[name] for name in inputs["compounds"]),
            interactions=self.srk_interactions,
            advanced_data=self.data,
        )
        self.assertIsInstance(system, SoaveRedlichKwongAdvancedSystem)
        evaluated = self.data.evaluated(
            inputs["temperature_k"], inputs["pressure_pa"], self.srk_interactions
        )
        self.assertTrue(
            math.isclose(
                evaluated.get("Mercury", "N-pentane"),
                -0.0623,
                rel_tol=1.0e-14,
            )
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
                golden["outputs"][f"{phase}_fugacity_coefficients"],
            ):
                self.assertTrue(math.isclose(actual, expected, rel_tol=2.0e-14))


if __name__ == "__main__":
    unittest.main()
