import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class Phase0ArtifactsTest(unittest.TestCase):
    def test_golden_schema_declares_auditable_case_sections(self):
        schema = json.loads((ROOT / "tests/golden/schema.json").read_text())

        self.assertEqual(schema["$schema"], "https://json-schema.org/draft/2020-12/schema")
        required = set(schema["required"])
        self.assertTrue({"schema_version", "case_id", "case_kind", "source", "inputs", "outputs"} <= required)

    def test_capture_script_uses_automation_and_records_property_units(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        for token in (
            "DWSIM.Automation.Automation3",
            "CalculateFlowsheet4",
            "GetProperties",
            "GetPropertyUnit",
            "ConvertTo-Json",
        ):
            self.assertIn(token, script)

    def test_capture_script_accepts_generic_dictionary_entries(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("[object]$Entry", script)
        self.assertNotIn("[System.Collections.DictionaryEntry]$Entry", script)

    def test_capture_script_preloads_portable_thermoc_assembly(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn('Join-Path $engineDirectory "ThermoCS\\ThermoCS.dll"', script)
        self.assertIn('[Reflection.Assembly]::LoadFrom($thermoCAssemblyPath)', script)

    def test_capture_script_reads_dwsim_objects_through_clr_reflection(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn('::AsPSObject($Object).BaseObject', script)
        self.assertIn('[object].GetMethod("GetType").Invoke', script)
        self.assertIn('InvokeMember', script)
        self.assertNotIn('$Object.GetType()', script)
        self.assertNotIn('$Object.PSObject.Properties[$Name]', script)

    def test_compatibility_record_requires_source_revision(self):
        text = (ROOT / "docs/compatibility.md").read_text()

        self.assertIn("GPLv3", text)
        self.assertIn("source revision", text.lower())
        self.assertIn("golden-case-1", text)


if __name__ == "__main__":
    unittest.main()
