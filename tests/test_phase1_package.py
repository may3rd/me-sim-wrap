import subprocess
import sys
import unittest
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


if __name__ == "__main__":
    unittest.main()
