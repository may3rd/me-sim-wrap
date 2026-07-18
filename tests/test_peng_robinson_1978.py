import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds, load_pr_interactions, load_srk_interactions
from mesim.errors import ValidationError
from mesim.thermo.peng_robinson import PengRobinson
from mesim.thermo.peng_robinson_1978 import PengRobinson1978
from mesim.thermo.systems import PENG_ROBINSON_1978, create_thermo_system


class PengRobinson1978Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = {
            compound.id: compound
            for compound in load_compounds(ROOT / "data/compounds/v1.json")
        }
        cls.interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")

    def test_high_acentric_factor_branch_matches_dwsim_equation(self):
        methanol = self.catalog["Methanol"]
        omega = methanol.acentric_factor.value
        self.assertGreater(omega, 0.491)
        expected = (
            0.379642
            + 1.48503 * omega
            - 0.164423 * omega**2
            + 0.016666 * omega**3
        )
        pr78 = PengRobinson1978(methanol).parameters(300.0)
        classic = PengRobinson(methanol).parameters(300.0)
        self.assertEqual(pr78.kappa, expected)
        self.assertNotEqual(pr78.kappa, classic.kappa)
        self.assertNotEqual(pr78.alpha, classic.alpha)

    def test_feed_phase_states_match_repeatable_dwsim_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/pr78-ethane-methanol-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (ROOT / "tests/golden/pr78-ethane-methanol-state-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(golden["case_id"], "pr78-ethane-methanol-state")
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.PengRobinson1978PropertyPackage",
        )
        inputs = golden["inputs"]
        system = create_thermo_system(
            PENG_ROBINSON_1978,
            compounds=tuple(self.catalog[name] for name in inputs["compounds"]),
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
                golden["outputs"][f"{phase}_fugacity_coefficients"],
            ):
                self.assertTrue(math.isclose(actual, expected, rel_tol=1.0e-13))

    def test_system_rejects_srk_interaction_identity(self):
        with self.assertRaises(ValidationError):
            create_thermo_system(
                PENG_ROBINSON_1978,
                compounds=(self.catalog["Methane"], self.catalog["Ethane"]),
                interactions=load_srk_interactions(
                    ROOT / "data/interactions/srk-v1.json"
                ),
            )


if __name__ == "__main__":
    unittest.main()
