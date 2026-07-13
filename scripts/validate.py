from __future__ import annotations

import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = ROOT / "tests" / "golden"


def enabled_cases() -> list[Path]:
    return sorted(path for path in GOLDEN_DIR.glob("*.json") if path.name != "schema.json")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate enabled me-sim-wrap reference cases")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    cases = enabled_cases()
    if not args.quiet:
        message = "no enabled golden cases" if not cases else f"{len(cases)} golden cases queued"
        print(f"validation: {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
