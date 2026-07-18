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
    WILSON_ACETONE_METHANOL,
    WilsonSystem,
    create_thermo_system,
)
from mesim.thermo.wilson import load_wilson_data, wilson_activity_coefficients


class WilsonTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / "data/interactions/wilson-v1.json"
        cls.data = load_wilson_data(cls.path)

    def test_complete_dwsim_table_and_scoped_volume_basis_are_frozen(self):
        self.assertEqual(self.data.source_revision, "9.0.5.0")
        self.assertEqual(len(self.data.interactions), 364)
        self.assertEqual(
            self.data.resource_sha256,
            "51523fafba912f74f7e484b74ae48bf8b94fd27e672ce518dd0ddfb2d9be710e",
        )
        acetone = self.data.compound("Acetone")
        methanol = self.data.compound("Methanol")
        self.assertEqual(acetone.cas, "67-64-1")
        self.assertEqual(methanol.cas, "67-56-1")
        self.assertTrue(
            math.isclose(
                acetone.molar_volume_m3_per_kmol,
                0.07380259101031655,
                rel_tol=1.0e-15,
            )
        )
        self.assertEqual(self.data.energy(acetone.cas, methanol.cas), -161.8813)
        self.assertEqual(self.data.energy(methanol.cas, acetone.cas), 583.1054)

    def test_activity_coefficients_match_repeatable_dwsim_golden(self):
        golden = json.loads(
            (ROOT / "tests/golden/wilson-acetone-methanol-state.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (
                ROOT
                / "tests/golden/wilson-acetone-methanol-state-repeat.json"
            ).read_text(encoding="utf-8-sig")
        )
        self.assertEqual(golden, repeat)
        self.assertEqual(
            golden["source"]["property_package_class"],
            "DWSIM.Thermodynamics.WilsonPropertyPackage",
        )
        inputs = golden["inputs"]
        actual = wilson_activity_coefficients(
            self.data,
            tuple(inputs["compounds"]),
            tuple(inputs["composition"]),
            inputs["temperature_k"],
        )
        for value, expected in zip(
            actual, golden["outputs"]["activity_coefficients"]
        ):
            self.assertTrue(math.isclose(value, expected, rel_tol=2.0e-15))

    def test_registered_system_preserves_wilson_result(self):
        system = create_thermo_system(
            WILSON_ACETONE_METHANOL,
            data=self.data,
            compound_ids=("Acetone", "Methanol"),
        )
        self.assertIsInstance(system, WilsonSystem)
        expected = wilson_activity_coefficients(
            self.data, system.compound_ids, (0.4, 0.6), 330.0
        )
        self.assertEqual(system.activity_coefficients((0.4, 0.6), 330.0), expected)

    def test_invalid_domain_and_units_are_rejected(self):
        with self.assertRaises(MissingCompoundData):
            self.data.compound("acetone")
        with self.assertRaises(ValidationError):
            wilson_activity_coefficients(
                self.data, ("Acetone", "Methanol"), (0.4, 0.7), 330.0
            )
        document = json.loads(self.path.read_text(encoding="utf-8-sig"))
        document["interactions"][0]["A12"]["unit"] = "J/mol"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaises(ValidationError):
                load_wilson_data(path)


if __name__ == "__main__":
    unittest.main()
