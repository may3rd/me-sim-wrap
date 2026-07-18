import json
import math
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.errors import MissingCompoundData, ValidationError
from mesim.thermo.systems import (
    UNIFAC_1_PROPANOL_WATER,
    UnifacSystem,
    create_thermo_system,
)
from mesim.thermo.unifac import load_unifac_data, unifac_activity_coefficients


class UnifacTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / "data/interactions/unifac-v1.json"
        cls.data = load_unifac_data(cls.path)

    def test_complete_runtime_group_and_interaction_domains_are_frozen(self):
        self.assertEqual(self.data.source_revision, "9.0.5.0")
        self.assertEqual(len(self.data.groups), 119)
        self.assertEqual(len(self.data.interactions), 1403)
        self.assertEqual(
            self.data.groups_sha256,
            "d7e2c18c1938f690b13e5566fcdd0fbb57f6837c8cdbece95317e9274158173d",
        )
        self.assertEqual(
            self.data.interactions_sha256,
            "e61e2324d3a2c10d62e70c85e9e0f6b5a7c5c2c187b3db5b47927f9f3c10d819",
        )
        propanol = self.data.compound("1-propanol")
        water = self.data.compound("Water")
        self.assertTrue(math.isclose(propanol.q, 3.128, rel_tol=1.0e-15))
        self.assertTrue(math.isclose(propanol.r, 3.2499, rel_tol=1.0e-15))
        self.assertEqual(
            tuple(item.secondary_id for item in propanol.group_surface_fractions),
            (1, 2, 15),
        )
        self.assertEqual(tuple(item.secondary_id for item in water.group_surface_fractions), (17,))
        self.assertEqual(self.data.interaction(1, 7), 1318.0)
        self.assertEqual(self.data.interaction(7, 5), -229.1)
        self.assertEqual(self.data.group(4).q, 0.0)

    def test_activity_coefficients_match_repeatable_dwsim_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/unifac-1-propanol-water-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (ROOT / "tests/golden/unifac-1-propanol-water-state-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.UNIFACPropertyPackage",
        )
        inputs = golden["inputs"]
        actual = unifac_activity_coefficients(
            self.data, tuple(inputs["compounds"]), tuple(inputs["composition"]),
            inputs["temperature_k"],
        )
        for value, expected in zip(actual, golden["outputs"]["activity_coefficients"]):
            self.assertTrue(math.isclose(value, expected, rel_tol=2.0e-13))

    def test_registered_system_preserves_unifac_result(self):
        system = create_thermo_system(
            UNIFAC_1_PROPANOL_WATER,
            data=self.data,
            compound_ids=("1-propanol", "Water"),
        )
        self.assertIsInstance(system, UnifacSystem)
        expected = unifac_activity_coefficients(
            self.data, system.compound_ids, (0.5, 0.5), 350.0
        )
        self.assertEqual(system.activity_coefficients((0.5, 0.5), 350.0), expected)

    def test_invalid_case_key_and_group_count_are_rejected(self):
        with self.assertRaises(MissingCompoundData):
            self.data.compound("1-Propanol")
        with self.assertRaises(ValidationError):
            unifac_activity_coefficients(
                self.data, ("1-propanol", "Water"), (0.6, 0.5), 350.0
            )
        document = json.loads(self.path.read_text(encoding="utf-8-sig"))
        document["groups"].pop()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaises(ValidationError):
                load_unifac_data(path)


if __name__ == "__main__":
    unittest.main()
