import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.errors import ValidationError
from mesim.reactions import load_reaction_data
from mesim.unitops.reactors import conversion_reactor


ROOT = Path(__file__).parents[1]
DWSIM_REACTOR_PHASE_AGGREGATE_REL_TOL = 2e-7


class ReactionDataTest(unittest.TestCase):
    def test_isomerisation_data_is_explicit_balanced_and_source_backed(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = data.reactions[0]
        thermo = {record.compound_id: record for record in data.thermochemistry}
        golden = json.loads((ROOT / "tests/golden/u4-conversion-reactor-isomerization.json").read_text(encoding="utf-8-sig"))
        captured = {record["id"]: record for record in golden["inputs"]["compounds"]}

        self.assertEqual(dict(thermo["N-butane"].elements), {"C": 4.0, "H": 10.0})
        self.assertEqual(dict(thermo["Isobutane"].elements), {"C": 4.0, "H": 10.0})
        self.assertEqual(dict(reaction.stoichiometry), {"N-butane": -1.0, "Isobutane": 1.0})
        self.assertEqual(reaction.reaction_heat_j_per_kmol, -9_200_000.0)
        self.assertEqual(reaction.conversion_fraction, 0.33)
        self.assertEqual(data.provenance.source_revision, "9.0.5.0")
        for compound_id in ("N-butane", "Isobutane"):
            record = thermo[compound_id]
            reference = captured[compound_id]
            molecular_weight = reference["molecular_weight"]["value"]
            self.assertEqual(dict(record.elements), reference["elements"])
            self.assertTrue(math.isclose(
                record.ideal_gas_formation_enthalpy_j_per_kmol,
                reference["ideal_gas_formation"]["enthalpy"]["value"] * molecular_weight * 1_000.0,
                rel_tol=1e-12,
            ))
            self.assertTrue(math.isclose(
                record.ideal_gas_formation_gibbs_energy_j_per_kmol,
                reference["ideal_gas_formation"]["gibbs_energy"]["value"] * molecular_weight * 1_000.0,
                rel_tol=1e-12,
            ))

    def test_reaction_loader_rejects_unbalanced_stoichiometry_and_inconsistent_heat(self):
        source = json.loads((ROOT / "data/reactions/v1.json").read_text())
        for mutation in ("stoichiometry", "heat"):
            broken = json.loads(json.dumps(source))
            if mutation == "stoichiometry":
                broken["reactions"][0]["stoichiometry"]["Isobutane"] = 2.0
            else:
                broken["reactions"][0]["reaction_heat"]["value"] = 0.0
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "reactions.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_reaction_data(path)


class ConversionReactorTest(unittest.TestCase):
    def test_conversion_reactor_matches_captured_dwsim_isomerisation(self):
        data = load_reaction_data(ROOT / "data/reactions/v1.json")
        reaction = data.reactions[0]
        golden = json.loads((ROOT / "tests/golden/u4-conversion-reactor-isomerization.json").read_text(encoding="utf-8-sig"))
        objects = {item["tag"]: item for item in golden["outputs"]["objects_after"]}
        properties = {
            tag: {item["property"]: item["value"]["value"] for item in objects[tag]["properties"]}
            for tag in ("R4", "Reactor Output1", "Reactor Output 2", "Conversion Reactor")
        }
        compounds = ("Propane", "Isobutane", "N-butane", "Isopentane")
        inlet = tuple(
            (compound, properties["R4"][f"PROP_MS_104/{compound}"] / 1_000.0)
            for compound in compounds
        )
        expected_outlet = {
            compound: math.fsum(
                properties[tag][f"PROP_MS_104/{compound}"] / 1_000.0
                for tag in ("Reactor Output1", "Reactor Output 2")
            )
            for compound in compounds
        }
        result = conversion_reactor(inlet, reaction)
        actual_outlet = dict(result.outlet_component_flows_kmol_s)

        self.assertEqual(properties["Conversion Reactor"]["isomerisation: Extent"], 33.0)
        self.assertEqual(properties["Conversion Reactor"]["N-butane: Conversion"], 33.0)
        self.assertEqual(result.extent_kmol_s, dict(inlet)["N-butane"] * 0.33)
        self.assertEqual(result.reaction_heat_w, result.extent_kmol_s * -9_200_000.0)
        for compound in compounds:
            self.assertTrue(math.isclose(
                actual_outlet[compound], expected_outlet[compound],
                rel_tol=DWSIM_REACTOR_PHASE_AGGREGATE_REL_TOL,
            ))
        inlet_total = math.fsum(flow for _, flow in inlet)
        expected_total = math.fsum(expected_outlet.values())
        self.assertTrue(math.isclose(result.total_molar_flow_kmol_s, inlet_total, rel_tol=0.0, abs_tol=1e-15))
        self.assertTrue(math.isclose(expected_total, inlet_total, rel_tol=1e-9))

        with ZipFile(ROOT / "tests/u4-conversion-reactor-isomerization.dwxmz") as archive:
            root = ElementTree.fromstring(archive.read(next(name for name in archive.namelist() if name.endswith(".xml"))))
        saved = root.find("./Reactions/Reaction")
        self.assertEqual(saved.findtext("ReactionType"), "Conversion")
        self.assertEqual(saved.findtext("BaseReactant"), reaction.base_reactant)
        self.assertEqual(float(saved.findtext("Expression")) / 100.0, reaction.conversion_fraction)
        self.assertTrue(math.isclose(float(saved.findtext("ReactionHeat")) * 1_000.0, reaction.reaction_heat_j_per_kmol, rel_tol=1e-12))

    def test_conversion_reactor_rejects_invalid_conversion_and_missing_base_flow(self):
        reaction = load_reaction_data(ROOT / "data/reactions/v1.json").reactions[0]
        with self.assertRaises(ValidationError):
            conversion_reactor((("N-butane", 1.0),), reaction, 1.01)
        with self.assertRaises(ValidationError):
            conversion_reactor((("Isobutane", 1.0),), reaction)


if __name__ == "__main__":
    unittest.main()
