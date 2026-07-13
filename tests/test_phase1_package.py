import subprocess
import sys
import tempfile
import unittest
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


class Phase1PackageTest(unittest.TestCase):
    def test_package_exports_domain_errors(self):
        from mesim import ConvergenceError, MissingCompoundData, OutOfRangeError, ValidationError

        self.assertTrue(issubclass(MissingCompoundData, ValidationError))
        self.assertTrue(issubclass(OutOfRangeError, ValidationError))
        self.assertTrue(issubclass(ConvergenceError, RuntimeError))

    def test_quiet_validation_succeeds_without_enabled_cases(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "validate.py"), "--quiet"],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")

    def test_validation_compares_cases_without_capture_timestamp(self):
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.json"
            second = Path(directory) / "second.json"
            first.write_text(
                json.dumps({"source": {"captured_utc": "2026-07-13T00:00:00Z"}, "value": 1}),
                encoding="utf-8-sig",
            )
            second.write_text(
                json.dumps({"source": {"captured_utc": "2026-07-13T00:00:01Z"}, "value": 1}),
                encoding="utf-8-sig",
            )

            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "validate.py"), "--compare", str(first), str(second)],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.returncode, 0, result.stderr)

    def test_validation_rejects_flowsheet_property_read_errors(self):
        case = {
            "case_kind": "flowsheet",
            "inputs": {"objects_before": []},
            "outputs": {
                "solve": {"executed": True, "success": True, "errors": []},
                "objects_after": [{"tag": "8", "properties": [{"read_error": "method not found"}]}],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "broken.json"
            path.write_text(json.dumps(case), encoding="utf-8")
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "validate.py"), "--compare", str(path), str(path)],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )

        self.assertEqual(result.returncode, 1)
        self.assertIn("property read error", result.stderr)


if __name__ == "__main__":
    unittest.main()
