import copy
import json
import math
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.errors import ValidationError
from mesim.flowsheet import RecycleConvergenceError, solve_recycle


ROOT = Path(__file__).parents[1]


class DirectSubstitutionRecycleTest(unittest.TestCase):
    def test_damped_direct_substitution_converges_with_complete_history(self):
        result = solve_recycle(
            lambda values: (1.0 + 0.5 * values[0], 10.0 + 0.25 * values[1]),
            initial_guess=(0.0, 0.0),
            scales=(2.0, 20.0),
            tolerances=(1.0e-10, 1.0e-10),
            damping=0.8,
        )

        self.assertTrue(math.isclose(result.values[0], 2.0, rel_tol=0.0, abs_tol=1.0e-9))
        self.assertTrue(math.isclose(result.values[1], 40.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-9))
        self.assertEqual(result.algorithm, "direct_substitution")
        self.assertEqual(tuple(item.iteration for item in result.history), tuple(range(1, len(result.history) + 1)))
        self.assertTrue(all(item.damping == 0.8 for item in result.history))
        self.assertTrue(all(
            later.scaled_norm <= earlier.scaled_norm
            for earlier, later in zip(result.history, result.history[1:])
        ))
        final = result.history[-1]
        self.assertEqual(result.values, final.calculated)
        self.assertTrue(all(abs(value) <= 1.0e-10 for value in final.residual))

    def test_failure_exposes_history_without_returning_partial_values(self):
        with self.assertRaises(RecycleConvergenceError) as caught:
            solve_recycle(
                lambda values: (values[0] + 1.0,),
                initial_guess=(0.0,), scales=(1.0,), tolerances=(1.0e-12,),
                max_iterations=3,
            )
        history = caught.exception.history
        self.assertEqual(len(history), 3)
        self.assertEqual(tuple(item.guess[0] for item in history), (0.0, 1.0, 2.0))
        self.assertTrue(all(item.residual == (1.0,) for item in history))

    def test_invalid_vectors_and_damping_are_rejected(self):
        with self.assertRaises(ValidationError):
            solve_recycle(lambda values: values, (), (), ())
        with self.assertRaises(ValidationError):
            solve_recycle(lambda values: values, (1.0,), (0.0,), (1.0e-3,))
        with self.assertRaises(ValidationError):
            solve_recycle(lambda values: values, (1.0,), (1.0,), (1.0e-3,), damping=1.1)
        with self.assertRaises(ValidationError):
            solve_recycle(lambda values: (math.nan,), (1.0,), (1.0,), (1.0e-3,))


class DWSIMMaterialRecycleTest(unittest.TestCase):
    def test_cavett_direct_substitution_recycle_is_repeatable_and_closed(self):
        golden = json.loads(
            (ROOT / "tests/golden/u8-material-recycle-cavett-pr-eos.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u8-material-recycle-cavett-pr-eos-repeat.json").read_text(encoding="utf-8-sig")
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertFalse(any(
            property_record["read_error"]
            for object_record in golden["outputs"]["objects_after"]
            for property_record in object_record["properties"]
        ))

        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("10", "1")
        }
        tear_properties = ("PROP_MS_0", "PROP_MS_1", "PROP_MS_2", "PROP_MS_7")
        calculated = tuple(properties["10"][name] for name in tear_properties)
        initial = tuple(properties["1"][name] for name in tear_properties)
        self.assertEqual(calculated, initial)
        for compound_property in (
            name for name in properties["10"] if name.startswith("PROP_MS_104/")
        ):
            self.assertEqual(properties["10"][compound_property], properties["1"][compound_property])

        result = solve_recycle(
            lambda _: calculated,
            initial_guess=initial,
            scales=(300.0, 2.0e6, 20.0, 300.0),
            tolerances=(1.0, 0.1, 0.01, 1.0),
            damping=1.0,
            max_iterations=100,
        )
        self.assertEqual(result.values, calculated)
        self.assertEqual(len(result.history), 1)
        self.assertEqual(result.history[0].residual, (0.0, 0.0, 0.0, 0.0))

        with ZipFile(ROOT / "tests/u8-material-recycle-cavett-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        saved = next(
            record for record in root.findall(".//SimulationObject")
            if record.findtext("Name") == "REC-f5ba6034-d804-4914-8b4b-c070058ba4c3"
        )
        self.assertEqual(saved.findtext("AccelerationMethod"), "None")
        self.assertEqual(saved.findtext("Converged"), "true")
        self.assertEqual(saved.findtext("IterationsTaken"), "1")
        self.assertEqual(saved.findtext("MaximumIterations"), "100")
        self.assertEqual(saved.findtext("ConvergenceParameters/Temperatura"), "1")
        self.assertEqual(saved.findtext("ConvergenceParameters/Pressao"), "0.1")
        self.assertEqual(saved.findtext("ConvergenceParameters/VazaoMassica"), "0.01")


if __name__ == "__main__":
    unittest.main()
