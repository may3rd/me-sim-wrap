import json
import math
import sys
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mesim.compounds import load_compounds, load_pr_interactions
from mesim.errors import MissingCompoundData, ValidationError


ROOT = Path(__file__).parents[1]


class CompoundDataTest(unittest.TestCase):
    def test_catalog_matches_captured_dwsim_records_and_is_immutable(self):
        compounds = load_compounds(ROOT / "data/compounds/v1.json")
        captured = json.loads(
            (ROOT / "tests/golden/compound-catalog.json").read_text(encoding="utf-8-sig")
        )["inputs"]["compounds"]

        self.assertEqual(len(compounds), 408)
        self.assertTrue({c["id"] for c in captured}.issubset({c.id for c in compounds}))
        for compound in (item for item in compounds if item.id in {c["id"] for c in captured}):
            reference = next(c for c in captured if c["id"] == compound.id)
            for field in ("name", "cas", "formula"):
                self.assertEqual(getattr(compound, field), reference[field])
            for field in ("molecular_weight", "critical_temperature", "critical_pressure", "acentric_factor", "normal_boiling_point"):
                self.assertEqual(getattr(compound, field).value, reference[field]["value"])
                self.assertEqual(getattr(compound, field).unit, reference[field]["unit"])
            self.assertEqual(compound.provenance.database, "ChemSep")
            self.assertEqual(
                compound.provenance.source,
                "dwsim-windows/DWSIM.Thermodynamics/Assets/Databases/chemsep1.xml",
            )
            self.assertEqual(compound.provenance.source_revision, "9.0.4")
            self.assertTrue(compound.provenance.imported_utc.endswith("Z"))

        equilibrium = json.loads(
            (ROOT / "tests/golden/u4-equilibrium-reactor-steam-reforming-pr-eos.json").read_text(encoding="utf-8-sig")
        )["inputs"]["compounds"]
        equilibrium_by_id = {record["id"]: record for record in equilibrium}
        for compound in (item for item in compounds if item.id in {"Water", "Carbon monoxide", "Carbon dioxide", "Hydrogen"}):
            reference = equilibrium_by_id[compound.id]
            self.assertEqual(compound.formula, reference["formula"])
            for field in ("molecular_weight", "critical_temperature", "critical_pressure", "acentric_factor", "normal_boiling_point"):
                self.assertEqual(getattr(compound, field).value, reference[field]["value"])

        nrtl = json.loads(
            (ROOT / "tests/golden/t4-nrtl-acetone-methanol-vle.json").read_text(encoding="utf-8-sig")
        )["inputs"]["compounds"]
        acetone_reference = next(record for record in nrtl if record["id"] == "Acetone")
        acetone = next(compound for compound in compounds if compound.id == "Acetone")
        self.assertEqual((acetone.name, acetone.cas, acetone.formula), ("Acetone", "67-64-1", "CH3COCH3"))
        for field in ("molecular_weight", "critical_temperature", "critical_pressure", "acentric_factor", "normal_boiling_point"):
            self.assertEqual(getattr(acetone, field).value, acetone_reference[field]["value"])
        with self.assertRaises(FrozenInstanceError):
            compounds[0].id = "changed"

    def test_catalog_rejects_duplicate_cas_and_missing_pr_property(self):
        valid = json.loads((ROOT / "data/compounds/v1.json").read_text())
        for mutation in ("duplicate", "missing"):
            broken = json.loads(json.dumps(valid))
            if mutation == "duplicate":
                broken["compounds"][1]["cas"] = broken["compounds"][0]["cas"]
            else:
                del broken["compounds"][0]["critical_pressure"]
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "compounds.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_compounds(path)

    def test_catalog_rejects_invalid_pr_numbers(self):
        valid = json.loads((ROOT / "data/compounds/v1.json").read_text())
        for field, value in (("critical_temperature", 0), ("critical_pressure", -1), ("molecular_weight", math.inf), ("acentric_factor", True)):
            broken = json.loads(json.dumps(valid))
            broken["compounds"][0][field]["value"] = value
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "compounds.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_compounds(path)

    def test_catalog_rejects_invalid_provenance(self):
        valid = json.loads((ROOT / "data/compounds/v1.json").read_text())
        for field, value in (("source", ""), ("imported_utc", "not-a-time")):
            broken = json.loads(json.dumps(valid))
            broken["compounds"][0]["provenance"][field] = value
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "compounds.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_compounds(path)

    def test_pr_interactions_are_symmetric_and_missing_is_explicit(self):
        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        self.assertEqual(interactions.provenance.source_revision, "9.0.5.0")
        self.assertTrue(interactions.provenance.imported_utc.endswith("Z"))
        self.assertEqual(interactions.get("Methane", "N-pentane"), 0.023)
        self.assertEqual(interactions.get("N-pentane", "Methane"), 0.023)
        self.assertEqual(interactions.get("Methane", "Methane"), 0.0)
        self.assertEqual(interactions.get("Methane", "Water"), 0.5)
        self.assertEqual(len(interactions.pairs), 197)
        with self.assertRaises(MissingCompoundData):
            interactions.get("Methane", "Bromine")

    def test_pr_interactions_match_first_entries_loaded_by_dwsim(self):
        interactions = load_pr_interactions(ROOT / "data/interactions/pr-v1.json")
        root = ElementTree.parse(
            ROOT
            / "dwsim-windows/DWSIM.Thermodynamics/Assets/Databases/chemsep1.xml"
        ).getroot()
        names = {
            compound.find("LibraryIndex").attrib["value"]: compound.find(
                "CompoundID"
            ).attrib["value"]
            for compound in root.iter("compound")
        }
        expected = {}
        source = ROOT / interactions.provenance.source
        for line in source.read_text(encoding="utf-8-sig").splitlines()[1:]:
            fields = line.split(";")
            if len(fields) >= 3 and fields[0] in names and fields[1] in names:
                expected.setdefault(frozenset((names[fields[0]], names[fields[1]])), float(fields[2]))
        actual = {frozenset((first, second)): value for first, second, value in interactions.pairs}
        for pair, value in actual.items():
            if value == 0.0 and pair not in expected:
                continue
            self.assertEqual(value, expected[pair])

    def test_unknown_schema_and_missing_pair_policy_are_rejected(self):
        cases = (
            (ROOT / "data/compounds/v1.json", "schema_version", "unknown", load_compounds),
            (ROOT / "data/interactions/pr-v1.json", "missing_pair_policy", "unknown", load_pr_interactions),
        )
        for source, key, value, loader in cases:
            data = json.loads(source.read_text())
            data[key] = value
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "data.json"
                path.write_text(json.dumps(data))
                with self.assertRaises(ValidationError):
                    loader(path)

    def test_invalid_pr_model_pairs_and_values_are_rejected(self):
        valid = json.loads((ROOT / "data/interactions/pr-v1.json").read_text())
        mutations = (
            ("model", lambda data: data.update(model="SRK")),
            (
                "self pair",
                lambda data: data["pairs"][0].update(
                    compound_2=data["pairs"][0]["compound_1"]
                ),
            ),
            ("blank ID", lambda data: data["pairs"][0].update(compound_1="")),
            ("non-finite", lambda data: data["pairs"][0].update(kij=math.inf)),
            ("wrong unit", lambda data: data["pairs"][0].update(unit="Pa")),
            ("missing unit", lambda data: data["pairs"][0].pop("unit")),
            ("blank provenance", lambda data: data["provenance"].update(source="")),
            ("bad import time", lambda data: data["provenance"].update(imported_utc="not-a-time")),
        )
        for _, mutate in mutations:
            broken = json.loads(json.dumps(valid))
            mutate(broken)
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "interactions.json"
                path.write_text(json.dumps(broken))
                with self.assertRaises(ValidationError):
                    load_pr_interactions(path)


if __name__ == "__main__":
    unittest.main()
