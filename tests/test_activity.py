import copy
import json
import math
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim import MissingCompoundData, ValidationError
from mesim.thermo.activity import (
    load_nrtl_vle_data,
    nrtl_activity_coefficients,
    nrtl_bubble_pressure,
    nrtl_dew_pressure,
)


ROOT = Path(__file__).parents[1]


class NRTLActivityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_nrtl_vle_data(ROOT / "data/correlations/nrtl-acetone-methanol-v1.json")
        cls.primary = json.loads(
            (ROOT / "tests/golden/t4-nrtl-acetone-methanol-vle.json").read_text(encoding="utf-8-sig")
        )
        cls.repeated = json.loads(
            (ROOT / "tests/golden/t4-nrtl-acetone-methanol-vle-repeat.json").read_text(encoding="utf-8-sig")
        )
        cls.record = cls.primary["outputs"]["objects_after"][0]

    @classmethod
    def value(cls, property_id: str) -> float:
        for prop in cls.record["properties"]:
            if prop["property"] == property_id:
                return prop["value"]["value"]
        raise AssertionError(f"golden lacks {property_id}")

    def test_reference_capture_is_repeatable_clean_and_physically_ordered(self):
        primary = copy.deepcopy(self.primary)
        repeated = copy.deepcopy(self.repeated)
        for document in (primary, repeated):
            document["source"].pop("captured_utc", None)
        self.assertEqual(primary, repeated)
        self.assertEqual(self.primary["source"]["property_package"], "NRTL")
        self.assertTrue(self.primary["source"]["bubble_dew_calculation"])
        self.assertTrue(self.primary["outputs"]["solve"]["success"])
        self.assertEqual(self.primary["outputs"]["solve"]["errors"], [])
        self.assertEqual(self.record["tag"], "HP Azeotrope")
        self.assertFalse(self.record["error"])
        self.assertTrue(all(prop["read_error"] is None for prop in self.record["properties"]))
        self.assertGreater(self.value("PROP_MS_126"), self.value("PROP_MS_127"))

    def test_saved_nrtl_equation_vector_and_dwsim_activity_coefficients(self):
        ids = ("Acetone", "Methanol")
        composition = tuple(self.value(f"PROP_MS_102/{compound}") for compound in ids)
        temperature_k = self.value("PROP_MS_0")
        gamma = nrtl_activity_coefficients(self.data, ids, composition, temperature_k)

        # Independently evaluated from the saved directed A12/A21/alpha records
        # and DWSIM's documented 1.98721 cal/mol/K NRTL equation constant.
        expected = (1.0653504217516205, 1.3002235244580072)
        for calculated, equation_value in zip(gamma, expected):
            self.assertTrue(math.isclose(calculated, equation_value, rel_tol=1e-12))
        for calculated, compound in zip(gamma, ids):
            captured = self.value(f"Activity Coefficient, Liquid Phase 1 / {compound}")
            self.assertTrue(math.isclose(calculated, captured, rel_tol=2.5e-3))

    def test_modified_raoult_bubble_and_dew_pressures_match_dwsim(self):
        ids = ("Acetone", "Methanol")
        composition = tuple(self.value(f"PROP_MS_102/{compound}") for compound in ids)
        temperature_k = self.value("PROP_MS_0")
        bubble = nrtl_bubble_pressure(self.data, ids, composition, temperature_k)
        dew = nrtl_dew_pressure(self.data, ids, composition, temperature_k)

        self.assertTrue(bubble.converged)
        self.assertTrue(dew.converged)
        self.assertGreater(bubble.pressure_pa, dew.pressure_pa)
        self.assertTrue(math.isclose(bubble.pressure_pa, self.value("PROP_MS_126"), rel_tol=2e-3))
        self.assertTrue(math.isclose(dew.pressure_pa, self.value("PROP_MS_127"), rel_tol=2e-3))
        self.assertTrue(math.isclose(math.fsum(bubble.vapor_composition), 1.0, abs_tol=1e-12))
        self.assertTrue(math.isclose(math.fsum(dew.liquid_composition), 1.0, abs_tol=1e-12))
        self.assertLessEqual(bubble.residual, 1e-12)
        self.assertLessEqual(dew.residual, 1e-12)

    def test_pure_limit_and_failure_contracts(self):
        bubble = nrtl_bubble_pressure(self.data, ("Acetone", "Methanol"), (1.0, 0.0), 388.0)
        dew = nrtl_dew_pressure(self.data, ("Acetone", "Methanol"), (1.0, 0.0), 388.0)
        self.assertTrue(math.isclose(bubble.pressure_pa, dew.pressure_pa, rel_tol=1e-12))
        self.assertTrue(math.isclose(bubble.activity_coefficients[0], 1.0, abs_tol=1e-12))

        exhausted = nrtl_dew_pressure(
            self.data,
            ("Acetone", "Methanol"),
            (0.7, 0.3),
            388.0,
            max_iterations=1,
            tolerance=1e-30,
        )
        self.assertFalse(exhausted.converged)
        self.assertEqual(exhausted.iterations, 1)
        self.assertIsNotNone(exhausted.failure_reason)

        with self.assertRaises(ValidationError):
            nrtl_activity_coefficients(self.data, ("Acetone", "Methanol"), (0.8, 0.3), 388.0)
        with self.assertRaises(ValidationError):
            nrtl_activity_coefficients(self.data, ("Acetone", "Methanol"), (0.7, 0.3), 0.0)
        with self.assertRaises(ValidationError):
            nrtl_activity_coefficients(self.data, ("Acetone", "Methanol"), (0.7, 0.3), 1e308)
        with self.assertRaises(ValidationError):
            nrtl_dew_pressure(self.data, ("Acetone", "Methanol"), (0.7, 0.3), 388.0, max_iterations=0)

        missing = replace(self.data, interactions=())
        with self.assertRaises(MissingCompoundData):
            nrtl_activity_coefficients(missing, ("Acetone", "Methanol"), (0.7, 0.3), 388.0)

    def test_loader_rejects_invalid_units_duplicates_ranges_and_provenance(self):
        source = json.loads(
            (ROOT / "data/correlations/nrtl-acetone-methanol-v1.json").read_text(encoding="utf-8")
        )
        mutations = []

        bad_unit = copy.deepcopy(source)
        bad_unit["interactions"][0]["A12"]["unit"] = "kJ/mol"
        mutations.append(bad_unit)

        duplicate = copy.deepcopy(source)
        duplicate["interactions"].append(copy.deepcopy(duplicate["interactions"][0]))
        mutations.append(duplicate)

        bad_range = copy.deepcopy(source)
        bad_range["vapor_pressure_correlations"][0]["minimum_k"] = -1.0
        mutations.append(bad_range)

        local_time = copy.deepcopy(source)
        local_time["provenance"]["imported_utc"] = "2026-07-18T10:17:02"
        mutations.append(local_time)

        for index, document in enumerate(mutations):
            with self.subTest(index=index), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "nrtl.json"
                path.write_text(json.dumps(document), encoding="utf-8")
                with self.assertRaises(ValidationError):
                    load_nrtl_vle_data(path)


if __name__ == "__main__":
    unittest.main()
