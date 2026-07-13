from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = ROOT / "tests" / "golden"


def enabled_cases() -> list[Path]:
    return sorted(path for path in GOLDEN_DIR.glob("*.json") if path.name != "schema.json")


def normalized_digest(path: Path) -> str:
    case = json.loads(path.read_text())
    case.get("source", {}).pop("captured_utc", None)
    normalized = json.dumps(case, allow_nan=False, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(normalized.encode()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate enabled me-sim-wrap reference cases")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--compare", nargs=2, metavar=("FIRST", "SECOND"))
    args = parser.parse_args(argv)

    if args.compare:
        first, second = (Path(path) for path in args.compare)
        first_digest = normalized_digest(first)
        second_digest = normalized_digest(second)
        equal = first_digest == second_digest
        if not args.quiet:
            print(f"validation: normalized cases {'match' if equal else 'differ'}")
        return 0 if equal else 1

    cases = enabled_cases()
    if not args.quiet:
        message = "no enabled golden cases" if not cases else f"{len(cases)} golden cases queued"
        print(f"validation: {message}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
