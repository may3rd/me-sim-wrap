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
from mesim.flowsheet import (
    AdjustConvergenceError,
    RecycleConvergenceError,
    solve_adjust,
    solve_energy_recycle,
    solve_recycle,
)


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


class DWSIMEnergyRecycleTest(unittest.TestCase):
    def test_scalar_energy_recycle_records_complete_watt_basis_history(self):
        result = solve_energy_recycle(
            lambda duty_w: 200_000.0 + 0.25 * duty_w,
            initial_duty_w=0.0,
            scale_w=300_000.0,
            tolerance_w=1.0e-6,
            damping=0.8,
        )

        self.assertTrue(math.isclose(result.values[0], 800_000.0 / 3.0, abs_tol=1.0e-6))
        self.assertEqual(result.algorithm, "direct_substitution")
        self.assertGreater(len(result.history), 1)
        self.assertEqual(
            tuple(item.iteration for item in result.history),
            tuple(range(1, len(result.history) + 1)),
        )
        self.assertTrue(all(item.damping == 0.8 for item in result.history))
        self.assertLessEqual(abs(result.history[-1].residual[0]), 1.0e-6)

    def test_turboexpander_energy_recycle_is_repeatable_and_closed(self):
        golden = json.loads(
            (ROOT / "tests/golden/u8-energy-recycle-turboexpander.json").read_text(
                encoding="utf-8-sig"
            )
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u8-energy-recycle-turboexpander-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertTrue(golden["outputs"]["solve"]["success"])
        self.assertFalse(golden["outputs"]["solve"]["errors"])
        self.assertFalse(any(
            property_record["read_error"]
            for object_record in golden["outputs"]["objects_after"]
            for property_record in object_record["properties"]
        ))

        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in record["properties"]}
            for tag, record in objects.items()
        }
        inlet_kw = properties["ESTR-004"]["PROP_ES_0"]
        outlet_kw = properties["ESTR-003"]["PROP_ES_0"]
        expander_kw = properties["X-1 (TE)"]["PROP_TU_3"]
        compressor_kw = properties["C-1 (TE)"]["PROP_CO_3"]
        tolerance_kw = properties["EREC-116"]["PROP_ER_1"]
        residual_kw = properties["EREC-116"]["PROP_ER_2"]

        self.assertEqual(inlet_kw, outlet_kw)
        self.assertTrue(math.isclose(
            residual_kw, expander_kw - compressor_kw, rel_tol=0.0, abs_tol=1.0e-10,
        ))
        self.assertLessEqual(abs(residual_kw), tolerance_kw)

        result = solve_energy_recycle(
            lambda _: inlet_kw * 1000.0,
            initial_duty_w=outlet_kw * 1000.0,
            scale_w=400_000.0,
            tolerance_w=tolerance_kw * 1000.0,
            max_iterations=100,
        )
        self.assertEqual(result.values, (inlet_kw * 1000.0,))
        self.assertEqual(len(result.history), 1)
        self.assertEqual(result.history[0].residual, (0.0,))

        with ZipFile(ROOT / "tests/u8-energy-recycle-turboexpander.dwxmz") as archive:
            root = ElementTree.fromstring(
                archive.read(next(name for name in archive.namelist() if name.endswith(".xml")))
            )
        saved = next(
            record for record in root.findall(".//SimulationObject")
            if record.findtext("Type") == "DWSIM.UnitOperations.SpecialOps.EnergyRecycle"
        )
        self.assertEqual(saved.findtext("AccelerationMethod"), "None")
        self.assertEqual(saved.findtext("IterationsTaken"), "1")
        self.assertEqual(saved.findtext("MaximumIterations"), "100")
        self.assertLessEqual(abs(float(saved.find("ConvHist").attrib["EnergyE"])), tolerance_kw)


class DWSIMAdjustTest(unittest.TestCase):
    def test_bounded_newton_adjust_records_complete_history(self):
        result = solve_adjust(
            lambda manipulated: 1.0 + 2.0 * manipulated,
            target=5.0,
            initial_guess=0.5,
            lower_bound=0.0,
            upper_bound=4.0,
            controlled_scale=5.0,
            tolerance=1.0e-12,
            step_size=0.1,
        )

        self.assertEqual(result.algorithm, "bounded_newton")
        self.assertTrue(math.isclose(result.manipulated, 2.0, abs_tol=1.0e-12))
        self.assertTrue(math.isclose(result.controlled, 5.0, abs_tol=1.0e-12))
        self.assertEqual(tuple(item.iteration for item in result.history), (1, 2))
        self.assertTrue(math.isclose(result.history[0].derivative, 2.0, abs_tol=1.0e-12))
        self.assertEqual(result.history[-1].step, 0.0)
        self.assertLessEqual(abs(result.history[-1].residual), 1.0e-12)

    def test_adjust_failure_returns_no_partial_result_and_preserves_history(self):
        with self.assertRaises(AdjustConvergenceError) as caught:
            solve_adjust(
                lambda _: 1.0,
                target=2.0,
                initial_guess=0.5,
                lower_bound=0.0,
                upper_bound=1.0,
                controlled_scale=2.0,
                tolerance=1.0e-12,
                step_size=0.1,
            )
        self.assertEqual(len(caught.exception.history), 1)
        self.assertEqual(caught.exception.history[0].residual, 1.0)
        self.assertEqual(caught.exception.history[0].derivative, 0.0)

    def test_adjust_rejects_invalid_bounds_and_controlled_values(self):
        arguments = dict(
            target=1.0,
            initial_guess=0.5,
            lower_bound=0.0,
            upper_bound=1.0,
            controlled_scale=1.0,
            tolerance=1.0e-6,
            step_size=0.1,
        )
        with self.assertRaises(ValidationError):
            solve_adjust(None, **arguments)
        with self.assertRaises(ValidationError):
            solve_adjust(lambda value: value, **(arguments | {"lower_bound": 1.0}))
        with self.assertRaises(ValidationError):
            solve_adjust(lambda _: math.nan, **arguments)

    def test_biodiesel_adjust_is_repeatable_and_meets_target(self):
        golden = json.loads(
            (ROOT / "tests/golden/u8-adjust-biodiesel-nrtl.json").read_text(encoding="utf-8-sig")
        )
        repeat = json.loads(
            (ROOT / "tests/golden/u8-adjust-biodiesel-nrtl-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        normalized_golden, normalized_repeat = copy.deepcopy(golden), copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertTrue(golden["outputs"]["solve"]["success"])
        self.assertFalse(golden["outputs"]["solve"]["errors"])
        self.assertFalse(any(
            property_record["read_error"]
            for object_record in golden["outputs"]["objects_after"]
            for property_record in object_record["properties"]
        ))

        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in record["properties"]}
            for tag, record in objects.items()
        }
        manipulated_kg_s = properties["EtOH"]["PROP_MS_2"]
        controlled_mol_s = properties["Etanol"]["PROP_MS_104/Ethanol_BD"]
        target_mol_s = properties["ADJ-000"]["AdjustValue"]
        tolerance_mol_s = properties["ADJ-000"]["Tolerance"]
        self.assertLessEqual(abs(controlled_mol_s - target_mol_s), tolerance_mol_s)

        result = solve_adjust(
            lambda _: controlled_mol_s,
            target=target_mol_s,
            initial_guess=manipulated_kg_s,
            lower_bound=0.0,
            upper_bound=0.1,
            controlled_scale=target_mol_s,
            tolerance=tolerance_mol_s,
            step_size=0.001,
        )
        self.assertEqual(result.manipulated, manipulated_kg_s)
        self.assertEqual(result.controlled, controlled_mol_s)
        self.assertEqual(len(result.history), 1)

        with ZipFile(ROOT / "tests/u8-adjust-biodiesel-nrtl.dwxmz") as archive:
            root = ElementTree.fromstring(
                archive.read(next(name for name in archive.namelist() if name.endswith(".xml")))
            )
        saved = next(
            record for record in root.findall(".//SimulationObject")
            if record.findtext("Type") == "DWSIM.UnitOperations.SpecialOps.Adjust"
        )
        self.assertEqual(saved.findtext("SimultaneousAdjust"), "true")
        self.assertEqual(saved.findtext("AdjustValue"), "1.66666666666667")
        self.assertEqual(saved.findtext("Tolerance"), "0.0001")
        self.assertEqual(saved.find("ManipulatedObjectData").attrib["Name"], "EtOH")
        self.assertEqual(saved.find("ManipulatedObjectData").attrib["Property"], "PROP_MS_2")
        self.assertEqual(saved.find("ControlledObjectData").attrib["Name"], "Etanol")
        self.assertEqual(
            saved.find("ControlledObjectData").attrib["Property"],
            "PROP_MS_104/Ethanol_BD",
        )


if __name__ == "__main__":
    unittest.main()
