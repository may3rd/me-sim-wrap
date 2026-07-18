"""Generate the frozen pure-component correlation datasets from ChemSep XML."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from xml.etree import ElementTree


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "dwsim-windows/DWSIM.Thermodynamics/Assets/Databases/chemsep1.xml"
CORRELATIONS = ROOT / "data/correlations"
COMPOUNDS = ROOT / "data/compounds/v1.json"
PR_SOURCE = ROOT / "dwsim-windows/DWSIM.Thermodynamics/Assets/pr_ip.dat"
PR_INTERACTIONS = ROOT / "data/interactions/pr-v1.json"
SOURCE_PATH = SOURCE.relative_to(ROOT).as_posix()
REVISION = "9.0.4"
PRIMARY_COMPOUNDS = (
    "Methane",
    "Ethane",
    "Propane",
    "N-butane",
    "N-pentane",
    "Nitrogen",
    "Argon",
    "Oxygen",
    "Methanol",
    "Water",
    "Carbon monoxide",
    "Carbon dioxide",
    "Hydrogen",
    "Acetone",
)
SUPPORTED_EQUATIONS = {
    "IdealGasHeatCapacityCp": {1, 16, 100},
    "VaporPressure": {10, 101},
    "LiquidDensity": {105, 106},
    "LiquidHeatCapacityCp": {3, 4, 16, 100},
    "HeatOfVaporization": {106},
    "LiquidViscosity": {10, 16, 101},
    "VaporViscosity": {2, 3, 16, 102},
    "LiquidThermalConductivity": {3, 16, 100},
    "VaporThermalConductivity": {3, 16, 102},
    "SurfaceTension": {2, 16, 106, 116},
}
APPROVED_ZERO_PAIRS = (
    ("Hydrogen", "Water"),
    ("Carbon monoxide", "Carbon dioxide"),
    ("Carbon monoxide", "Water"),
    ("Hydrogen", "Argon"),
    ("Hydrogen", "Oxygen"),
    ("Hydrogen", "Methanol"),
    ("Nitrogen", "Water"),
    ("Carbon monoxide", "Argon"),
    ("Carbon monoxide", "Oxygen"),
    ("Carbon monoxide", "Methanol"),
    ("Argon", "Water"),
    ("Argon", "Methanol"),
    ("Argon", "Carbon dioxide"),
    ("Oxygen", "Water"),
    ("Oxygen", "Methanol"),
    ("Oxygen", "Carbon dioxide"),
)


def _number(text: str) -> int | float:
    value = float(text)
    return int(value) if value.is_integer() else value


def _correlation(node, unit: str, *, include_e: bool = True) -> dict:
    values = {child.tag: child.attrib["value"] for child in node}
    keys = "ABCDE" if include_e else "ABCD"
    return {
        "equation": int(values["eqno"]),
        **{key: _number(values.get(key, "0")) for key in keys},
        "minimum_k": _number(values["Tmin"]),
        "maximum_k": _number(values["Tmax"]),
        "unit": unit,
    }


def _equation(compound, tag: str) -> int | None:
    node = compound.find(tag)
    equation = None if node is None else node.find("eqno")
    if equation is None or equation.attrib.get("value") in (None, ""):
        return None
    return int(equation.attrib["value"])


def _supported(compound) -> bool:
    for tag, equations in SUPPORTED_EQUATIONS.items():
        node = compound.find(tag)
        if node is None or _equation(compound, tag) not in equations:
            return False
        values = {child.tag: child.attrib.get("value") for child in node}
        if any(values.get(key) in (None, "") for key in ("A", "B", "C", "D", "Tmin", "Tmax")):
            return False
    return True


def _source_records() -> tuple[list[str], dict[str, object]]:
    root = ElementTree.parse(SOURCE).getroot()
    by_id = {}
    for compound in root.iter("compound"):
        identifier = compound.find("CompoundID")
        if identifier is not None and _supported(compound):
            by_id[identifier.attrib["value"]] = compound
    missing = [compound_id for compound_id in PRIMARY_COMPOUNDS if compound_id not in by_id]
    if missing:
        raise RuntimeError(f"ChemSep is missing catalog compounds: {', '.join(missing)}")
    compound_ids = list(PRIMARY_COMPOUNDS)
    compound_ids.extend(name for name in by_id if name not in PRIMARY_COMPOUNDS)
    return compound_ids, by_id


def _datasets() -> dict[Path, dict]:
    compound_ids, compounds = _source_records()
    source = {
        "source": SOURCE_PATH,
        "source_revision": REVISION,
    }
    catalog = {
        "schema_version": "compound-data-1",
        "compounds": [],
    }
    ideal = {
        "schema_version": "ideal-correlations-1",
        "provenance": {**source, "imported_utc": "2026-07-13T07:24:14Z"},
        "correlations": [],
    }
    transport = {
        "schema_version": "transport-correlations-3",
        "provenance": {**source, "imported_utc": "2026-07-18T00:00:00Z"},
        "correlations": [],
    }
    saturated = {
        "schema_version": "saturated-liquid-correlations-1",
        "provenance": {**source, "imported_utc": "2026-07-18T00:00:00Z"},
        "correlations": [],
    }
    for compound_id in compound_ids:
        compound = compounds[compound_id]
        catalog["compounds"].append(
            {
                "id": compound_id,
                "name": compound_id,
                "cas": compound.find("CAS").attrib["value"],
                "formula": compound.find("StructureFormula").attrib["value"],
                "molecular_weight": {
                    "value": _number(compound.find("MolecularWeight").attrib["value"]),
                    "unit": "kg/kmol",
                },
                "critical_temperature": {
                    "value": _number(compound.find("CriticalTemperature").attrib["value"]),
                    "unit": "K",
                },
                "critical_pressure": {
                    "value": _number(compound.find("CriticalPressure").attrib["value"]),
                    "unit": "Pa",
                },
                "acentric_factor": {
                    "value": _number(compound.find("AcentricityFactor").attrib["value"]),
                    "unit": "dimensionless",
                },
                "normal_boiling_point": {
                    "value": _number(
                        compound.find("NormalBoilingPointTemperature").attrib["value"]
                    ),
                    "unit": "K",
                },
                "provenance": {
                    "database": "ChemSep",
                    **source,
                    "imported_utc": "2026-07-18T00:00:00Z",
                },
            }
        )
        ideal["correlations"].append(
            {
                "compound_id": compound_id,
                "heat_capacity": _correlation(
                    compound.find("IdealGasHeatCapacityCp"), "J/kmol/K"
                ),
                "vapor_pressure": _correlation(compound.find("VaporPressure"), "Pa"),
            }
        )
        transport["correlations"].append(
            {
                "compound_id": compound_id,
                "critical_volume": {
                    "value": _number(compound.find("CriticalVolume").attrib["value"]),
                    "unit": "m3/kmol",
                },
                "vapor_viscosity": _correlation(
                    compound.find("VaporViscosity"), "Pa.s"
                ),
                "vapor_thermal_conductivity": _correlation(
                    compound.find("VaporThermalConductivity"),
                    "W/m/K",
                ),
                "liquid_viscosity": _correlation(
                    compound.find("LiquidViscosity"), "Pa.s"
                ),
                "liquid_thermal_conductivity": _correlation(
                    compound.find("LiquidThermalConductivity"), "W/m/K"
                ),
            }
        )
        saturated["correlations"].append(
            {
                "compound_id": compound_id,
                "critical_temperature": {
                    "value": _number(
                        compound.find("CriticalTemperature").attrib["value"]
                    ),
                    "unit": "K",
                },
                "liquid_density": _correlation(
                    compound.find("LiquidDensity"), "kmol/m3"
                ),
                "liquid_heat_capacity": _correlation(
                    compound.find("LiquidHeatCapacityCp"), "J/kmol/K"
                ),
                "heat_of_vaporization": _correlation(
                    compound.find("HeatOfVaporization"), "J/kmol"
                ),
                "surface_tension": _correlation(
                    compound.find("SurfaceTension"), "N/m"
                ),
            }
        )
    interactions = _pr_interactions(compound_ids, compounds)
    exclusions = _exclusions()
    return {
        COMPOUNDS: catalog,
        CORRELATIONS / "ideal-v1.json": ideal,
        CORRELATIONS / "transport-v1.json": transport,
        CORRELATIONS / "saturated-liquid-v1.json": saturated,
        CORRELATIONS / "chemsep-exclusions-v1.json": exclusions,
        PR_INTERACTIONS: interactions,
    }


def _exclusions() -> dict:
    root = ElementTree.parse(SOURCE).getroot()
    records = []
    source_count = 0
    for compound in root.iter("compound"):
        source_count += 1
        if _supported(compound):
            continue
        missing_fields = []
        unsupported_equations = []
        for tag, equations in SUPPORTED_EQUATIONS.items():
            node = compound.find(tag)
            if node is None:
                missing_fields.append(tag)
                continue
            values = {child.tag: child.attrib.get("value") for child in node}
            for key in ("eqno", "A", "B", "C", "D", "Tmin", "Tmax"):
                if values.get(key) in (None, ""):
                    missing_fields.append(f"{tag}.{key}")
            equation = _equation(compound, tag)
            if equation is not None and equation not in equations:
                unsupported_equations.append(
                    {"property": tag, "equation": equation}
                )
        records.append(
            {
                "compound_id": compound.find("CompoundID").attrib["value"],
                "missing_fields": sorted(set(missing_fields)),
                "unsupported_equations": unsupported_equations,
            }
        )
    return {
        "schema_version": "chemsep-exclusions-1",
        "provenance": {
            "source": SOURCE_PATH,
            "source_revision": REVISION,
            "imported_utc": "2026-07-18T00:00:00Z",
        },
        "source_compound_count": source_count,
        "supported_compound_count": source_count - len(records),
        "excluded_compounds": records,
    }


def _pr_interactions(compound_ids: list[str], compounds: dict[str, object]) -> dict:
    supported = set(compound_ids)
    by_index = {
        compound.find("LibraryIndex").attrib["value"]: compound_id
        for compound_id, compound in compounds.items()
    }
    pairs = []
    seen = set()
    for line in PR_SOURCE.read_text(encoding="utf-8-sig").splitlines():
        fields = line.split(";")
        if len(fields) < 3 or fields[0] not in by_index or fields[1] not in by_index:
            continue
        first, second = by_index[fields[0]], by_index[fields[1]]
        key = frozenset((first, second))
        if first not in supported or second not in supported or key in seen:
            continue
        seen.add(key)
        pairs.append(
            {
                "compound_1": first,
                "compound_2": second,
                "kij": _number(fields[2]),
                "unit": "dimensionless",
            }
        )
    for first, second in APPROVED_ZERO_PAIRS:
        key = frozenset((first, second))
        if key not in seen:
            seen.add(key)
            pairs.append(
                {
                    "compound_1": first,
                    "compound_2": second,
                    "kij": 0,
                    "unit": "dimensionless",
                }
            )
    return {
        "schema_version": "pr-interactions-1",
        "model": "Peng-Robinson",
        "missing_pair_policy": "error",
        "provenance": {
            "source": PR_SOURCE.relative_to(ROOT).as_posix(),
            "source_revision": "9.0.5.0",
            "selection": (
                "first entry loaded by DWSIM for each supported pair; explicit zero "
                "for absent accepted steam-reforming and methanol-column pairs"
            ),
            "imported_utc": "2026-07-18T00:00:00Z",
        },
        "pairs": pairs,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="replace generated datasets")
    args = parser.parse_args()
    drift = []
    for path, data in _datasets().items():
        rendered = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
        if args.write:
            path.write_text(rendered, encoding="utf-8")
        elif not path.exists() or path.read_text(encoding="utf-8-sig") != rendered:
            drift.append(path.relative_to(ROOT).as_posix())
    if drift:
        print("Generated correlation data is stale: " + ", ".join(drift))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
