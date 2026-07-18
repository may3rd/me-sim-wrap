import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds, load_pr_interactions, load_srk_interactions
from mesim.errors import ValidationError
from mesim.thermo.systems import (
    PENG_ROBINSON_LEE_KESLER,
    PengRobinsonLeeKeslerSystem,
    create_thermo_system,
)


class PengRobinsonLeeKeslerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = {
            compound.id: compound
            for compound in load_compounds(ROOT / "data/compounds/v1.json")
        }
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        cls.golden = json.loads(
            (ROOT / "tests/golden/pr-lk-methane-ethane-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        cls.repeat = json.loads(
            (ROOT / "tests/golden/pr-lk-methane-ethane-state-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        inputs = cls.golden["inputs"]
        cls.system = create_thermo_system(
            PENG_ROBINSON_LEE_KESLER,
            compounds=tuple(cls.catalog[name] for name in inputs["compounds"]),
            interactions=cls.interactions,
        )

    def test_golden_is_repeatable_and_identifies_direct_desktop_package(self):
        self.assertEqual(self.golden, self.repeat)
        self.assertEqual(self.golden["case_id"], "pr-lk-methane-ethane-state")
        self.assertEqual(
            self.golden["source"]["property_package"],
            "Peng-Robinson / Lee-Kesler (PR/LK)",
        )
        self.assertEqual(
            self.golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.PengRobinsonLKPropertyPackage",
        )
        self.assertEqual(
            self.golden["source"]["property_package_construction"],
            "direct-class-over-case-compound-domain",
        )

    def test_feed_phase_fugacities_match_dwsim_pr_lk(self):
        inputs = self.golden["inputs"]
        for phase in ("liquid", "vapor"):
            state = self.system.state(
                tuple(inputs["composition"]),
                inputs["temperature_k"],
                inputs["pressure_pa"],
                phase,
            )
            for actual, expected in zip(
                state.fugacity_coefficients,
                self.golden["outputs"][f"{phase}_fugacity_coefficients"],
            ):
                self.assertTrue(math.isclose(actual, expected, rel_tol=2.0e-6))

    def test_tp_flash_matches_dwsim_phase_split(self):
        inputs = self.golden["inputs"]
        result = self.system.tp_flash(
            tuple(inputs["composition"]),
            inputs["temperature_k"],
            inputs["pressure_pa"],
        )
        self.assertTrue(result.report.converged)
        self.assertEqual(result.phase, "two-phase")
        self.assertTrue(
            math.isclose(
                result.vapor_fraction,
                self.golden["outputs"]["vapor_fraction"],
                rel_tol=5.0e-7,
            )
        )
        for actual, expected in zip(
            result.liquid_composition,
            self.golden["outputs"]["liquid_composition"],
        ):
            self.assertTrue(math.isclose(actual, expected, rel_tol=2.0e-5))
        for actual, expected in zip(
            result.vapor_composition,
            self.golden["outputs"]["vapor_composition"],
        ):
            self.assertTrue(math.isclose(actual, expected, rel_tol=2.0e-5))

    def test_system_rejects_srk_interaction_identity(self):
        with self.assertRaises(ValidationError):
            create_thermo_system(
                PENG_ROBINSON_LEE_KESLER,
                compounds=(self.catalog["Methane"], self.catalog["Ethane"]),
                interactions=load_srk_interactions(
                    ROOT / "data/interactions/srk-v1.json"
                ),
            )

    def test_runtime_type_is_distinct_from_classic_pr(self):
        self.assertIsInstance(self.system, PengRobinsonLeeKeslerSystem)
        self.assertFalse(hasattr(self.system, "enthalpy"))


if __name__ == "__main__":
    unittest.main()
