import copy
import json
import math
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.compounds import load_compounds, load_pr_interactions
from mesim.errors import ValidationError
from mesim.unitops.columns import (
    EquilibriumStageState,
    StageStream,
    ColumnConvergenceError,
    column_balance_residuals,
    equilibrium_stage_residuals,
    fixed_k_sum_rates_absorber,
    shortcut_column,
)


class EquilibriumStageBalanceTest(unittest.TestCase):
    def test_closed_binary_stage_reports_zero_mesh_residuals(self):
        feed = StageStream(10.0, (0.5, 0.5), 180.0)
        state = EquilibriumStageState(
            temperature_k=350.0,
            pressure_pa=101325.0,
            liquid=StageStream(6.0, (0.7, 0.3), 100.0),
            vapor=StageStream(4.0, (0.2, 0.8), 300.0),
        )

        residuals = equilibrium_stage_residuals(
            state,
            incoming_streams=(feed,),
            equilibrium_ratios=(2.0 / 7.0, 8.0 / 3.0),
        )

        self.assertTrue(all(abs(value) <= 1.0e-15 for value in residuals.component_material_kmol_s))
        self.assertTrue(all(abs(value) <= 1.0e-15 for value in residuals.phase_equilibrium))
        self.assertEqual(residuals.liquid_summation, 0.0)
        self.assertEqual(residuals.vapor_summation, 0.0)
        self.assertEqual(residuals.energy_w, 0.0)
        self.assertTrue(residuals.is_closed(1.0e-12, 1.0e-12, 1.0e-12, 1.0e-9))

    def test_trial_state_preserves_each_independent_residual(self):
        state = EquilibriumStageState(
            350.0,
            101325.0,
            StageStream(6.0, (0.8, 0.3), 100.0),
            StageStream(4.0, (0.1, 0.8), 250.0),
        )
        residuals = equilibrium_stage_residuals(
            state,
            incoming_streams=(StageStream(10.0, (0.5, 0.5), 180.0),),
            equilibrium_ratios=(0.25, 2.0),
            heat_duty_w=25.0,
        )

        self.assertTrue(math.isclose(
            residuals.component_material_kmol_s[0], -0.2, abs_tol=1.0e-15,
        ))
        self.assertEqual(residuals.component_material_kmol_s[1], 0.0)
        self.assertTrue(math.isclose(residuals.phase_equilibrium[0], -0.1, abs_tol=1.0e-15))
        self.assertTrue(math.isclose(residuals.phase_equilibrium[1], 0.2, abs_tol=1.0e-15))
        self.assertTrue(math.isclose(residuals.liquid_summation, 0.1, abs_tol=1.0e-15))
        self.assertTrue(math.isclose(residuals.vapor_summation, -0.1, abs_tol=1.0e-15))
        self.assertEqual(residuals.energy_w, 225.0)
        self.assertFalse(residuals.is_closed(1.0e-6, 1.0e-6, 1.0e-6, 1.0e-3))

    def test_invalid_stage_vectors_and_tolerances_are_rejected(self):
        with self.assertRaises(ValidationError):
            StageStream(-1.0, (1.0,), 0.0)
        with self.assertRaises(ValidationError):
            EquilibriumStageState(
                300.0, 101325.0,
                StageStream(1.0, (1.0,), 0.0),
                StageStream(1.0, (0.5, 0.5), 0.0),
            )
        state = EquilibriumStageState(
            300.0, 101325.0,
            StageStream(1.0, (1.0,), 0.0),
            StageStream(1.0, (1.0,), 0.0),
        )
        with self.assertRaises(ValidationError):
            equilibrium_stage_residuals(state, (), (1.0,))
        with self.assertRaises(ValidationError):
            equilibrium_stage_residuals(state, (StageStream(2.0, (1.0,), 0.0),), (0.0,))


class ShortcutColumnParityTest(unittest.TestCase):
    NAMES = (
        "Methanol", "Water", "Hydrogen", "Carbon dioxide",
        "Carbon monoxide", "Nitrogen", "Oxygen", "Argon",
    )

    @classmethod
    def setUpClass(cls):
        cls.golden = json.loads(
            (ROOT / "tests/golden/u6-shortcut-column-methanol-pr-eos.json").read_text(
                encoding="utf-8-sig"
            )
        )
        cls.objects = {
            record["tag"]: record
            for record in cls.golden["outputs"]["objects_after"]
        }
        cls.properties = {
            tag: {
                item["property"]: item["value"]["value"]
                for item in record["properties"]
            }
            for tag, record in cls.objects.items()
        }

    def test_dwsim_reference_is_repeatable_solved_and_error_free(self):
        repeat = json.loads(
            (ROOT / "tests/golden/u6-shortcut-column-methanol-pr-eos-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        normalized_golden = copy.deepcopy(self.golden)
        normalized_repeat = copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertTrue(self.golden["outputs"]["solve"]["success"])
        self.assertFalse(self.golden["outputs"]["solve"]["errors"])
        self.assertFalse(any(
            record["error"] or prop["read_error"]
            for record in self.golden["outputs"]["objects_after"]
            for prop in record["properties"]
        ))

    def test_frozen_compounds_and_pr_pairs_cover_reference_case(self):
        compounds = {record.id: record for record in load_compounds(ROOT / "data/compounds/v1.json")}
        captured = {record["id"]: record for record in self.golden["inputs"]["compounds"]}
        self.assertEqual(set(captured), set(self.NAMES))
        for name in self.NAMES:
            expected = captured[name]
            actual = compounds[name]
            self.assertEqual(actual.cas, expected["cas"])
            self.assertEqual(actual.formula, expected["formula"])
            self.assertEqual(actual.molecular_weight.value, expected["molecular_weight"]["value"])
            self.assertEqual(actual.critical_temperature.value, expected["critical_temperature"]["value"])
            self.assertEqual(actual.critical_pressure.value, expected["critical_pressure"]["value"])
            self.assertEqual(actual.acentric_factor.value, expected["acentric_factor"]["value"])
            self.assertEqual(actual.normal_boiling_point.value, expected["normal_boiling_point"]["value"])

        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        for first_index, first in enumerate(self.NAMES):
            for second in self.NAMES[first_index + 1:]:
                self.assertTrue(math.isfinite(interactions.get(first, second)))

    def test_dwsim_fug_shortcut_results_match(self):
        feed = self.properties["MSTR-084"]
        column = self.properties["SC-078"]
        top = self.properties["Methanol Product"]
        bottom = self.properties["MSTR-080"]
        feed_flow_mol_s = feed["PROP_MS_3"]
        feed_component_flows = tuple(
            feed[f"PROP_MS_104/{name}"] for name in self.NAMES
        )
        composition = tuple(value / feed_flow_mol_s for value in feed_component_flows)

        # Calculate the reference relative-volatility vector from DWSIM's
        # stated Fenske distribution equation. This preserves the property
        # package's shortcut K-vector independently of our PR flash solver.
        minimum_stages = column["PROP_SC_6"]
        top_flow = top["PROP_MS_3"]
        bottom_flow = bottom["PROP_MS_3"]
        top_fractions = tuple(
            top[f"PROP_MS_104/{name}"] / top_flow for name in self.NAMES
        )
        bottom_fractions = tuple(
            bottom[f"PROP_MS_104/{name}"] / bottom_flow for name in self.NAMES
        )
        distribution_constant = math.log10(top_fractions[1] / bottom_fractions[1])
        relative_volatilities = tuple(
            1.0 if index == 1 else 10.0 ** (
                (math.log10(top_fractions[index] / bottom_fractions[index])
                 - distribution_constant)
                / minimum_stages
            )
            for index in range(len(self.NAMES))
        )
        feed_liquid_fraction = (
            column["PROP_SC_8"] - column["PROP_SC_9"]
        ) / feed_flow_mol_s

        result = shortcut_column(
            feed_flow_kmol_s=feed_flow_mol_s / 1000.0,
            feed_mole_fractions=composition,
            relative_volatilities=relative_volatilities,
            light_key_index=0,
            heavy_key_index=1,
            heavy_key_distillate_fraction=column["PROP_SC_1"],
            light_key_bottoms_fraction=column["PROP_SC_2"],
            reflux_ratio=column["PROP_SC_0"],
            feed_liquid_fraction=feed_liquid_fraction,
        )

        self.assertEqual(result.iterations, 1)
        self.assertTrue(math.isclose(
            result.distillate_flow_kmol_s * 1000.0, top_flow, abs_tol=2.0e-15,
        ))
        self.assertTrue(math.isclose(
            result.bottoms_flow_kmol_s * 1000.0, bottom_flow, abs_tol=3.0e-14,
        ))
        for actual, expected in zip(result.distillate_mole_fractions, top_fractions):
            self.assertTrue(math.isclose(actual, expected, rel_tol=2.0e-13, abs_tol=1.0e-15))
        for actual, expected in zip(result.bottoms_mole_fractions, bottom_fractions):
            self.assertTrue(math.isclose(actual, expected, rel_tol=3.0e-13, abs_tol=1.0e-60))
        self.assertTrue(math.isclose(result.minimum_stages, column["PROP_SC_6"], abs_tol=2.0e-14))
        self.assertTrue(math.isclose(result.minimum_reflux_ratio, column["PROP_SC_5"], abs_tol=2.0e-6))
        self.assertTrue(math.isclose(result.actual_stages, column["PROP_SC_14"], abs_tol=3.0e-6))
        self.assertTrue(math.isclose(result.feed_stage, column["PROP_SC_7"], abs_tol=3.0e-6))
        flow_fields = (
            (result.stripping_liquid_flow_kmol_s, "PROP_SC_8"),
            (result.rectifying_liquid_flow_kmol_s, "PROP_SC_9"),
            (result.stripping_vapor_flow_kmol_s, "PROP_SC_10"),
            (result.rectifying_vapor_flow_kmol_s, "PROP_SC_11"),
        )
        for actual, property_id in flow_fields:
            self.assertTrue(math.isclose(
                actual * 1000.0, column[property_id], abs_tol=3.0e-15,
            ))
        self.assertEqual(column["PROP_SC_12"], self.properties["ESTR-081"]["PROP_ES_0"])
        self.assertEqual(column["PROP_SC_13"], self.properties["ESTR-082"]["PROP_ES_0"])

    def test_saved_reference_mode_fields_are_authoritative(self):
        with ZipFile(ROOT / "tests/u6-shortcut-column-methanol-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(
                archive.read(next(name for name in archive.namelist() if name.endswith(".xml")))
            )
        graphic_ids = {
            graphic.findtext("Tag"): graphic.findtext("Name")
            for graphic in root.findall(".//GraphicObject")
        }
        saved_objects = {
            record.findtext("Name"): record for record in root.findall(".//SimulationObject")
        }
        saved = saved_objects[graphic_ids["SC-078"]]
        self.assertEqual(saved.findtext("m_lightkey"), "Methanol")
        self.assertEqual(saved.findtext("m_heavykey"), "Water")
        self.assertEqual(saved.findtext("m_lightkeymolarfrac"), "0.1")
        self.assertEqual(saved.findtext("m_heavykeymolarfrac"), "0.001")
        self.assertEqual(saved.findtext("m_refluxratio"), "1.5")
        self.assertEqual(saved.findtext("m_condenserpressure"), "101325")
        self.assertEqual(saved.findtext("m_boilerpressure"), "101325")
        self.assertEqual(saved.findtext("condtype"), "TotalCond")

    def test_invalid_or_unsupported_shortcut_cases_are_rejected(self):
        with self.assertRaises(ValidationError):
            shortcut_column(1.0, (0.5, 0.5), (1.0, 2.0), 0, 0, 0.01, 0.1, 1.5, 1.0)
        with self.assertRaises(ValidationError):
            shortcut_column(1.0, (0.5, 0.5), (1.0, 2.0), 0, 1, 0.01, 0.1, 1.5, 1.0)


class AbsorberGoldenGateTest(unittest.TestCase):
    NAMES = ("Metano", "Etano", "Propano", "nOctano", "nNonano", "nDecano")

    def test_repeatable_absorber_closes_material_and_preserves_saved_solver(self):
        path = ROOT / "tests/golden/u6-absorber-simple-pr-eos.json"
        golden = json.loads(path.read_text(encoding="utf-8-sig"))
        repeat = json.loads(
            (ROOT / "tests/golden/u6-absorber-simple-pr-eos-repeat.json").read_text(
                encoding="utf-8-sig"
            )
        )
        normalized_golden = copy.deepcopy(golden)
        normalized_repeat = copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertTrue(golden["outputs"]["solve"]["success"])
        self.assertFalse(golden["outputs"]["solve"]["errors"])
        self.assertFalse(any(
            record["error"] or prop["read_error"]
            for record in golden["outputs"]["objects_after"]
            for prop in record["properties"]
        ))

        objects = {record["tag"]: record for record in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in record["properties"]}
            for tag, record in objects.items()
        }
        stream_vector = lambda tag: tuple(
            properties[tag][f"PROP_MS_104/{name}"] / 1000.0 for name in self.NAMES
        )
        residuals = column_balance_residuals(
            (stream_vector("Gas"), stream_vector("Liquid Solvent")),
            (stream_vector("Residue Gas"), stream_vector("Solvent + Solute")),
        )
        self.assertTrue(residuals.is_closed(2.0e-12))
        self.assertGreater(
            1.0 - properties["Residue Gas"]["PROP_MS_104/Propano"]
            / properties["Gas"]["PROP_MS_104/Propano"],
            0.999999999,
        )

        with ZipFile(ROOT / "tests/u6-absorber-simple-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(
                archive.read(next(name for name in archive.namelist() if name.endswith(".xml")))
            )
        graphic_ids = {
            graphic.findtext("Tag"): graphic.findtext("Name")
            for graphic in root.findall(".//GraphicObject")
        }
        saved = next(
            record for record in root.findall(".//SimulationObject")
            if record.findtext("Name") == graphic_ids["ABS-000"]
        )
        self.assertEqual(saved.findtext("OperationMode"), "Absorber")
        self.assertEqual(saved.findtext("NumberOfStages"), "12")
        self.assertEqual(saved.findtext("SolvingMethodName"), "Burningham-Otto (Sum Rates)")
        self.assertEqual(saved.findtext("MaxIterations"), "100")
        self.assertEqual(saved.findtext("ExternalLoopTolerance"), "1E-10")
        self.assertEqual(saved.findtext("InternalLoopTolerance"), "1E-10")

    def test_column_balance_rejects_inconsistent_vectors(self):
        with self.assertRaises(ValidationError):
            column_balance_residuals(((1.0, 2.0),), ((1.0,),))

    def test_fixed_k_sum_rates_reproduces_dwsim_stage_profiles(self):
        golden = json.loads(
            (ROOT / "tests/golden/u6-absorber-simple-pr-eos.json").read_text(
                encoding="utf-8-sig"
            )
        )
        properties = {
            record["tag"]: {
                item["property"]: item["value"]["value"]
                for item in record["properties"]
            }
            for record in golden["outputs"]["objects_after"]
        }
        with ZipFile(ROOT / "tests/u6-absorber-simple-pr-eos.dwxmz") as archive:
            root = ElementTree.fromstring(
                archive.read(next(name for name in archive.namelist() if name.endswith(".xml")))
            )
        graphic_ids = {
            graphic.findtext("Tag"): graphic.findtext("Name")
            for graphic in root.findall(".//GraphicObject")
        }
        saved = next(
            record for record in root.findall(".//SimulationObject")
            if record.findtext("Name") == graphic_ids["ABS-000"]
        )
        results = saved.find("Results")

        def vector(tag):
            return tuple(float(value) / 1000.0 for value in results.findtext(tag).strip("{}").split(";"))

        def matrix(tag):
            return tuple(
                tuple(float(value) for value in row.text.strip("{}").split(";"))
                for row in results.find(tag)
            )

        feeds = [[0.0] * len(self.NAMES) for _ in range(12)]
        for stage, tag in ((0, "Liquid Solvent"), (11, "Gas")):
            feeds[stage] = [
                properties[tag][f"PROP_MS_104/{name}"] / 1000.0
                for name in self.NAMES
            ]
        profile = fixed_k_sum_rates_absorber(
            tuple(tuple(row) for row in feeds),
            matrix("Kf"),
            vector("L0"),
            vector("V0"),
            flow_tolerance_kmol_s=1.0e-15,
        )
        self.assertLess(len(profile.history), 100)
        for actual, expected in zip(profile.liquid_flows_kmol_s, vector("Lf")):
            self.assertTrue(math.isclose(actual, expected, abs_tol=5.2e-7))
        for actual, expected in zip(profile.vapor_flows_kmol_s, vector("Vf")):
            self.assertTrue(math.isclose(actual, expected, abs_tol=1.4e-7))
        for actual_row, expected_row in zip(profile.liquid_mole_fractions, matrix("xf")):
            for actual, expected in zip(actual_row, expected_row):
                self.assertTrue(math.isclose(actual, expected, abs_tol=5.0e-7))
        for actual_row, expected_row in zip(profile.vapor_mole_fractions, matrix("yf")):
            for actual, expected in zip(actual_row, expected_row):
                self.assertTrue(math.isclose(actual, expected, abs_tol=5.0e-6))

    def test_fixed_k_sum_rates_validates_and_preserves_failure_history(self):
        with self.assertRaises(ValidationError):
            fixed_k_sum_rates_absorber(((1.0,),), ((1.0,),), (1.0,), (1.0,))
        feeds = ((0.5, 0.0), (0.0, 0.5))
        with self.assertRaises(ColumnConvergenceError) as caught:
            fixed_k_sum_rates_absorber(
                feeds, ((2.0, 0.5), (2.0, 0.5)), (1.0, 1.0), (1.0, 1.0),
                flow_tolerance_kmol_s=1.0e-30, maximum_iterations=1,
            )
        self.assertEqual(len(caught.exception.history), 1)


class ReboiledAbsorberGoldenGateTest(unittest.TestCase):
    NAMES = ("Methanol", "Acetone")

    @classmethod
    def setUpClass(cls):
        cls.golden = json.loads(
            (ROOT / "tests/golden/u6-reboiled-absorber-acetone-nrtl.json").read_text(
                encoding="utf-8-sig"
            )
        )
        cls.properties = {
            record["tag"]: {
                item["property"]: item["value"]["value"]
                for item in record["properties"]
            }
            for record in cls.golden["outputs"]["objects_after"]
        }

    def test_reference_is_repeatable_solved_and_materially_closed(self):
        repeat = json.loads(
            (
                ROOT
                / "tests/golden/u6-reboiled-absorber-acetone-nrtl-repeat.json"
            ).read_text(encoding="utf-8-sig")
        )
        normalized_golden = copy.deepcopy(self.golden)
        normalized_repeat = copy.deepcopy(repeat)
        normalized_golden["source"].pop("captured_utc")
        normalized_repeat["source"].pop("captured_utc")
        self.assertEqual(normalized_golden, normalized_repeat)
        self.assertTrue(self.golden["outputs"]["solve"]["success"])
        self.assertFalse(self.golden["outputs"]["solve"]["errors"])
        self.assertFalse(any(
            record["error"] or prop["read_error"]
            for record in self.golden["outputs"]["objects_after"]
            for prop in record["properties"]
        ))

        stream_vector = lambda tag: tuple(
            self.properties[tag][f"PROP_MS_104/{name}"] / 1000.0
            for name in self.NAMES
        )
        residuals = column_balance_residuals(
            (stream_vector("HP Feed"),),
            (stream_vector("HP Azeotrope"), stream_vector("Acetone Rich Product")),
        )
        self.assertTrue(residuals.is_closed(1.5e-9))
        self.assertEqual(
            self.properties["HP Feed"]["PROP_MS_3"],
            self.properties["HP Azeotrope"]["PROP_MS_3"]
            + self.properties["Acetone Rich Product"]["PROP_MS_3"],
        )

    def test_saved_mode_and_vapor_product_mapping_are_authoritative(self):
        with ZipFile(
            ROOT / "tests/u6-reboiled-absorber-acetone-nrtl.dwxmz"
        ) as archive:
            root = ElementTree.fromstring(
                archive.read(
                    next(name for name in archive.namelist() if name.endswith(".xml"))
                )
            )
        graphic_ids = {
            graphic.findtext("Tag"): graphic.findtext("Name")
            for graphic in root.findall(".//GraphicObject")
        }
        saved = next(
            record
            for record in root.findall(".//SimulationObject")
            if record.findtext("Name") == graphic_ids["Acetone Column (6 atm)"]
        )
        saved_tags = {value: key for key, value in graphic_ids.items()}
        streams = {
            stream.findtext("StreamBehavior"): saved_tags[stream.findtext("StreamID")]
            for stream in saved.findall("./MaterialStreams/MaterialStream")
        }
        self.assertEqual(saved.findtext("ReboiledAbsorber"), "true")
        self.assertEqual(saved.findtext("NumberOfStages"), "20")
        self.assertEqual(saved.findtext("SolvingMethodName"), "Wang-Henke (Bubble Point)")
        self.assertEqual(saved.findtext("InitialEstimatesProvider"), "Internal 2 (Experimental)")
        self.assertEqual(saved.findtext("MaxIterations"), "1000")
        self.assertEqual(saved.findtext("ExternalLoopTolerance"), "0.001")
        self.assertEqual(saved.findtext("InternalLoopTolerance"), "0.001")
        self.assertNotIn("Distillate", streams)
        self.assertEqual(streams["OverheadVapor"], "HP Azeotrope")
        self.assertEqual(streams["BottomsLiquid"], "Acetone Rich Product")


if __name__ == "__main__":
    unittest.main()
