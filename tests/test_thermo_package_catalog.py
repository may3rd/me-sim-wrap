import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mesim.errors import ValidationError
from mesim.thermo.package_catalog import load_thermodynamic_package_catalog
from mesim.thermo.systems import THERMO_SYSTEM_CONSTRUCTORS


class ThermodynamicPackageCatalogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = ROOT / "data/thermo-packages/dwsim-v1.json"
        cls.catalog = load_thermodynamic_package_catalog(cls.path)

    def test_catalog_freezes_all_dwsim_builtins(self):
        self.assertEqual(self.catalog.catalog_id, "dwsim-builtins-v1")
        self.assertEqual(self.catalog.source.version, "9.0.5.0")
        self.assertEqual(len(self.catalog.packages), 29)
        self.assertEqual(
            {record.id for record in self.catalog.packages},
            {
                "coolprop", "coolprop-incompressible-pure",
                "coolprop-incompressible-mixture", "steam-tables", "seawater",
                "peng-robinson", "prsv2-m", "prsv2-vl",
                "soave-redlich-kwong", "peng-robinson-lee-kesler", "unifac",
                "unifac-ll", "modified-unifac-dortmund", "modified-unifac-nist",
                "nrtl", "wilson", "uniquac", "chao-seader", "grayson-streed",
                "ideal-raoult", "lee-kesler-plocker", "reaktoro",
                "ideal-aqueous-electrolyte", "black-oil", "gerg-2008", "pc-saft",
                "peng-robinson-1978", "peng-robinson-1978-advanced",
                "soave-redlich-kwong-advanced",
            },
        )

    def test_extracted_model_ids_match_the_runtime_registry(self):
        extracted = {
            model_id
            for record in self.catalog.packages
            for model_id in record.mesim_model_ids
        }
        self.assertEqual(extracted, set(THERMO_SYSTEM_CONSTRUCTORS))
        self.assertEqual(self.catalog.package("peng-robinson").extraction_status, "partial")
        with self.assertRaises(ValidationError):
            self.catalog.package("not-a-package")

    def test_loader_rejects_an_untracked_or_duplicate_model_id(self):
        data = json.loads(self.path.read_text(encoding="utf-8"))
        data["packages"][0]["mesim_model_ids"] = ["peng-robinson-classic"]
        data["packages"][0]["extraction_status"] = "partial"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            with self.assertRaises(ValidationError):
                load_thermodynamic_package_catalog(path)


if __name__ == "__main__":
    unittest.main()
