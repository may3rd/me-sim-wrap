import sys
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src"))

import mesim
from mesim.api import DATA, app


class ReleaseArtifactTest(unittest.TestCase):
    def test_package_and_api_versions_match_release(self):
        metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual(metadata["project"]["version"], "0.1.0")
        self.assertEqual(mesim.__version__, "0.1.0")
        self.assertEqual(app.version, "0.1.0")
        self.assertTrue((DATA / "compounds/v1.json").is_file())

    def test_release_documents_and_changelog_exist(self):
        for relative in (
            "docs/compatibility.md",
            "docs/model-limitations.md",
            "docs/data-provenance.md",
            "CHANGELOG.md",
        ):
            text = (ROOT / relative).read_text(encoding="utf-8")
            self.assertGreater(len(text), 500, relative)
        self.assertIn("## [0.1.0] - 2026-07-18", (ROOT / "CHANGELOG.md").read_text())

    def test_support_matrix_uses_only_normative_statuses(self):
        compatibility = (ROOT / "docs/compatibility.md").read_text(encoding="utf-8")
        section = compatibility.split("## Release support matrix", 1)[1].split("## Source revision", 1)[0]
        rows = [line for line in section.splitlines() if line.startswith("| ")][2:]
        self.assertGreaterEqual(len(rows), 10)
        statuses = {row.split("|")[2].strip() for row in rows}
        self.assertLessEqual(statuses, {"unsupported", "partial", "verified", "verified-with-difference"})

    def test_release_workflow_covers_required_architectures(self):
        workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
        for runner in ("ubuntu-24.04", "ubuntu-24.04-arm", "macos-14"):
            self.assertIn(runner, workflow)
        self.assertIn("python -m unittest discover -s tests", workflow)
        self.assertIn("python scripts/validate.py --quiet", workflow)
        self.assertIn("docker build", workflow)
        self.assertIn("mesim.__version__ == '0.1.0'", workflow)


if __name__ == "__main__":
    unittest.main()
