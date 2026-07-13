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

    def test_capture_script_preloads_portable_thermoc_assembly(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn('Join-Path $engineDirectory "ThermoCS\\ThermoCS.dll"', script)
        self.assertIn('[Reflection.Assembly]::LoadFrom($thermoCAssemblyPath)', script)

    def test_capture_script_reads_dwsim_objects_through_clr_reflection(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("Add-Type -TypeDefinition", script)
        self.assertIn("public static class DwsimCaptureReflection", script)
        self.assertIn('[DwsimCaptureReflection]::Get($object, "Name")', script)
        self.assertIn('$Object.PSObject.Properties[$Name]', script)
        self.assertNotIn("function Get-ClrBaseObject", script)

    def test_compatibility_record_requires_source_revision(self):
        text = (ROOT / "docs/compatibility.md").read_text()

        self.assertIn("GPLv3", text)
        self.assertIn("source revision", text.lower())
        self.assertIn("golden-case-1", text)

    def test_capture_script_reads_compounds_without_dictionary_entry_conversion(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("AvailableCompounds.ContainsKey", script)
        self.assertIn("$automation.AvailableCompounds[", script)
        self.assertIn("[object]$Constant", script)
        self.assertNotIn("[object]$Entry", script)
        self.assertNotIn("[System.Collections.DictionaryEntry]$Entry", script)

    def test_capture_script_preloads_portable_thermoc_assembly(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn('-Filter "ThermoCS.dll"', script)
        self.assertIn("-Recurse", script)
        self.assertIn("$thermoCAssemblyPath.FullName", script)

    def test_capture_script_reads_dwsim_objects_through_clr_reflection(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("Add-Type -TypeDefinition", script)
        self.assertIn("public static class DwsimCaptureReflection", script)
        self.assertIn("[DwsimCaptureReflection]::Get(", script)
        self.assertIn("$object,", script)
        self.assertIn('"Name"', script)
        self.assertIn("$Object.PSObject.Properties[$Name]", script)
        self.assertNotIn("function Get-ClrBaseObject", script)

if __name__ == "__main__":
    unittest.main()
