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
    UNIQUAC_1_PROPANOL_WATER,
    UniquacSystem,
    create_thermo_system,
)
from mesim.thermo.uniquac import load_uniquac_data, uniquac_activity_coefficients


class UniquacTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / "data/interactions/uniquac-v1.json"
        cls.data = load_uniquac_data(cls.path)

    def test_complete_source_table_and_resolved_binary_basis_are_frozen(self):
        self.assertEqual(self.data.source_revision, "9.0.5.0")
        self.assertEqual(len(self.data.source_interactions), 376)
        source_pairs = [
            (record.first_chemsep_id, record.second_chemsep_id)
            for record in self.data.source_interactions
        ]
        self.assertEqual(len(set(source_pairs)), 359)
        self.assertEqual(
            self.data.resource_sha256,
            "c5e36cf466e3f683a4db23555d33e460180f6df4d84968fac7e223f88c1df81f",
        )
        propanol = self.data.compound("1-propanol")
        water = self.data.compound("Water")
        self.assertEqual((propanol.chemsep_id, propanol.q, propanol.r), ("1103", 2.51, 2.78))
        self.assertEqual((water.chemsep_id, water.q, water.r), ("1921", 1.4, 0.92))
        self.assertEqual(self.data.interaction.a12, 190.5947)
        self.assertEqual(self.data.interaction.a21, 290.554)

    def test_activity_coefficients_match_repeatable_dwsim_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/uniquac-1-propanol-water-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (
                ROOT / "tests/golden/uniquac-1-propanol-water-state-repeat.json"
            ).read_text(encoding="utf-8-sig")
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.PropertyPackages.UNIQUACPropertyPackage",
        )
        inputs = golden["inputs"]
        actual = uniquac_activity_coefficients(
            self.data,
            tuple(inputs["compounds"]),
            tuple(inputs["composition"]),
            inputs["temperature_k"],
        )
        for value, expected in zip(actual, golden["outputs"]["activity_coefficients"]):
            self.assertTrue(math.isclose(value, expected, rel_tol=2.0e-13))

    def test_registered_system_preserves_uniquac_result(self):
        system = create_thermo_system(
            UNIQUAC_1_PROPANOL_WATER,
            data=self.data,
            compound_ids=("1-propanol", "Water"),
        )
        self.assertIsInstance(system, UniquacSystem)
        expected = uniquac_activity_coefficients(
            self.data, system.compound_ids, (0.5, 0.5), 350.0
        )
        self.assertEqual(system.activity_coefficients((0.5, 0.5), 350.0), expected)

    def test_invalid_domain_and_source_count_are_rejected(self):
        with self.assertRaises(MissingCompoundData):
            self.data.compound("1-Propanol")
        with self.assertRaises(ValidationError):
            uniquac_activity_coefficients(
                self.data, ("1-propanol", "Water"), (0.5, 0.6), 350.0
            )
        document = json.loads(self.path.read_text(encoding="utf-8-sig"))
        document["source_interactions"].pop()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaises(ValidationError):
                load_uniquac_data(path)


if __name__ == "__main__":
    unittest.main()
