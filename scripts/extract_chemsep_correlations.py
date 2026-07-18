"""Generate the frozen pure-component correlation datasets from ChemSep XML."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from xml.etree import ElementTree


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "dwsim-windows/DWSIM.Thermodynamics/Assets/Databases/chemsep1.xml"
CATALOG = ROOT / "data/compounds/v1.json"
CORRELATIONS = ROOT / "data/correlations"
SOURCE_PATH = SOURCE.relative_to(ROOT).as_posix()
REVISION = "9.0.4"


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


def _source_records() -> tuple[list[str], dict[str, object]]:
    catalog = json.loads(CATALOG.read_text(encoding="utf-8-sig"))
    compound_ids = [record["id"] for record in catalog["compounds"]]
    root = ElementTree.parse(SOURCE).getroot()
    by_id = {}
    for compound in root.iter("compound"):
        identifier = compound.find("CompoundID")
        if identifier is not None and identifier.attrib.get("value") in compound_ids:
            by_id[identifier.attrib["value"]] = compound
    missing = [compound_id for compound_id in compound_ids if compound_id not in by_id]
    if missing:
        raise RuntimeError(f"ChemSep is missing catalog compounds: {', '.join(missing)}")
    return compound_ids, by_id


def _datasets() -> dict[Path, dict]:
    compound_ids, compounds = _source_records()
    source = {
        "source": SOURCE_PATH,
        "source_revision": REVISION,
    }
    ideal = {
        "schema_version": "ideal-correlations-1",
        "provenance": {**source, "imported_utc": "2026-07-13T07:24:14Z"},
        "correlations": [],
    }
    transport = {
        "schema_version": "transport-correlations-2",
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
                    compound.find("VaporViscosity"), "Pa.s", include_e=False
                ),
                "vapor_thermal_conductivity": _correlation(
                    compound.find("VaporThermalConductivity"),
                    "W/m/K",
                    include_e=False,
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
    return {
        CORRELATIONS / "ideal-v1.json": ideal,
        CORRELATIONS / "transport-v1.json": transport,
        CORRELATIONS / "saturated-liquid-v1.json": saturated,
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
