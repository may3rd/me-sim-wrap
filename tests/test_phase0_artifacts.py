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

    def test_capture_script_can_enable_dwsim_bubble_and_dew_calculation(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("[switch]$CalculateBubbleAndDewPoints", script)
        self.assertIn("SetFlashSetting", script)
        self.assertIn('"CalculateBubbleAndDewPoints"', script)
        self.assertIn("bubble_dew_calculation", script)
        self.assertIn("property_packages_updated", script)

        load = script.index("$flowsheet = $automation.LoadFlowsheet2(")
        setting = script.index("$propertyPackagesUpdated = [DwsimCaptureReflection]::SetFlashSetting(")
        before = script.index("$before = @(", load)
        solve = script.index("$automation.CalculateFlowsheet4(", load)
        self.assertLess(load, setting)
        self.assertLess(setting, before)
        self.assertLess(setting, solve)

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

    def test_capture_script_accepts_empty_unit_before_normalizing_it(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("[AllowEmptyString()]\n        [string]$Unit", script)

    def test_capture_script_serializes_non_finite_numbers_as_text(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("[double]::IsNaN", script)
        self.assertIn("[double]::IsInfinity", script)
        self.assertIn('value_type = "non_finite"', script)

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
        self.assertEqual(script.count('[object[]]@($PropertyName, $null)'), 2)

    def test_capture_script_records_saved_attached_utility_inputs(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn('"AttachedUtilities"', script)
        self.assertIn('"SaveData"', script)
        self.assertIn("Get-UtilityStates", script)

    def test_capture_script_reads_attached_utility_inputs_from_dwxmz(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("System.IO.Compression.ZipFile", script)
        self.assertIn("Get-SavedUtilityStates", script)
        self.assertIn("AttachedUtility", script)
        self.assertIn("ConvertFrom-Json", script)

    def test_capture_script_records_ideal_compound_reference_values(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("ideal_reference", script)
        self.assertIn('"GetIdealGasHeatCapacity"', script)
        self.assertIn('"GetVaporPressure"', script)

    def test_capture_script_uses_dwsim_canonical_compound_id(self):
        script = (ROOT / "scripts/capture_dwsim_reference.ps1").read_text()

        self.assertIn("id = $canonicalName", script)
        self.assertNotIn("id = $CompoundId", script)
        self.assertIn('"N-butane"', script)
        self.assertIn('"N-pentane"', script)
        self.assertNotIn('"n-Butane"', script)
        self.assertNotIn('"n-Pentane"', script)

if __name__ == "__main__":
    unittest.main()
