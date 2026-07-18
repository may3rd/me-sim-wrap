import json
import math
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.errors import ValidationError
from mesim.thermo.systems import (
    UNIFAC_LL_1_PROPANOL_WATER,
    UnifacLLSystem,
    create_thermo_system,
)
from mesim.thermo.unifac import load_unifac_data, unifac_activity_coefficients


class UnifacLLTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / "data/interactions/unifac-ll-v1.json"
        cls.data = load_unifac_data(cls.path)

    def test_complete_ll_interaction_domain_and_executable_basis_are_frozen(self):
        self.assertEqual(self.data.model, "UNIFAC-LL")
        self.assertEqual(len(self.data.groups), 119)
        self.assertEqual(len(self.data.interactions), 1467)
        self.assertEqual(
            self.data.interactions_sha256,
            "8f720e3b2e372399a0c7f7c709f650562788fb52b8d9fe95166ece8f4af0bda7",
        )
        propanol = self.data.compound("1-propanol")
        water = self.data.compound("Water")
        self.assertTrue(
            math.isclose(
                math.fsum(item.value for item in propanol.group_surface_fractions),
                1.374340949033392,
                rel_tol=1.0e-15,
            )
        )
        self.assertTrue(
            math.isclose(
                math.fsum(item.value for item in water.group_surface_fractions),
                0.9776536312849162,
                rel_tol=1.0e-15,
            )
        )
        self.assertEqual(self.data.interaction(1, 7), 310.7)
        self.assertEqual(self.data.interaction(7, 1), -131.9)

    def test_activity_coefficients_match_repeatable_dwsim_golden_exactly(self):
        golden = json.loads(
            (ROOT / "tests/golden/unifac-ll-1-propanol-water-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (ROOT / "tests/golden/unifac-ll-1-propanol-water-state-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.UNIFACLLPropertyPackage",
        )
        inputs = golden["inputs"]
        actual = unifac_activity_coefficients(
            self.data, tuple(inputs["compounds"]), tuple(inputs["composition"]),
            inputs["temperature_k"],
        )
        self.assertEqual(actual, tuple(golden["outputs"]["activity_coefficients"]))

    def test_registered_ll_system_rejects_original_unifac_data(self):
        system = create_thermo_system(
            UNIFAC_LL_1_PROPANOL_WATER,
            data=self.data,
            compound_ids=("1-propanol", "Water"),
        )
        self.assertIsInstance(system, UnifacLLSystem)
        self.assertEqual(
            system.activity_coefficients((0.5, 0.5), 350.0),
            (0.38278298255146526, 0.6363828049999769),
        )
        original = load_unifac_data(ROOT / "data/interactions/unifac-v1.json")
        with self.assertRaises(ValidationError):
            create_thermo_system(
                UNIFAC_LL_1_PROPANOL_WATER,
                data=original,
                compound_ids=("1-propanol", "Water"),
            )

    def test_invalid_ll_interaction_count_is_rejected(self):
        document = json.loads(self.path.read_text(encoding="utf-8-sig"))
        document["interactions"].pop()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaises(ValidationError):
                load_unifac_data(path)


if __name__ == "__main__":
    unittest.main()
