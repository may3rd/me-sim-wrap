import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from .errors import MissingCompoundData, ValidationError


@dataclass(frozen=True)
class Property:
    value: float
    unit: str


@dataclass(frozen=True)
class Provenance:
    database: str
    source: str
    source_revision: str
    imported_utc: str


@dataclass(frozen=True)
class InteractionProvenance:
    source: str
    source_revision: str
    selection: str
    imported_utc: str


@dataclass(frozen=True)
class Compound:
    id: str
    name: str
    cas: str
    formula: str
    molecular_weight: Property
    critical_temperature: Property
    critical_pressure: Property
    acentric_factor: Property
    normal_boiling_point: Property
    provenance: Provenance


@dataclass(frozen=True)
class PRInteractions:
    pairs: tuple[tuple[str, str, float], ...]
    default_zero: bool
    provenance: InteractionProvenance

    def get(self, first: str, second: str) -> float:
        if first == second:
            return 0.0
        for left, right, value in self.pairs:
            if {first, second} == {left, right}:
                return value
        if self.default_zero:
            return 0.0
        raise MissingCompoundData(f"missing PR interaction: {first}/{second}")


_UNITS = {
    "molecular_weight": "kg/kmol",
    "critical_temperature": "K",
    "critical_pressure": "Pa",
    "acentric_factor": "dimensionless",
    "normal_boiling_point": "K",
}


def _validate_provenance(provenance: Provenance | InteractionProvenance) -> None:
    if not all(isinstance(value, str) and value for value in vars(provenance).values()):
        raise ValidationError("provenance fields must be non-empty strings")
    try:
        imported = datetime.fromisoformat(provenance.imported_utc)
    except ValueError as error:
        raise ValidationError("provenance imported_utc must be an ISO 8601 timestamp") from error
    if not provenance.imported_utc.endswith("Z") or imported.utcoffset() != timedelta(0):
        raise ValidationError("provenance imported_utc must be UTC")


def load_compounds(path: str | Path) -> tuple[Compound, ...]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if data["schema_version"] != "compound-data-1":
            raise ValidationError("unsupported compound data schema")
        records = data["compounds"]
        compounds = tuple(
            Compound(
                id=record["id"],
                name=record["name"],
                cas=record["cas"],
                formula=record["formula"],
                **{key: Property(**record[key]) for key in _UNITS},
                provenance=Provenance(**record["provenance"]),
            )
            for record in records
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid compound data: {error}") from error
    if len({c.id for c in compounds}) != len(compounds) or len({c.cas for c in compounds}) != len(compounds):
        raise ValidationError("compound IDs and CAS numbers must be unique")
    for compound in compounds:
        _validate_provenance(compound.provenance)
        if not all(isinstance(getattr(compound, field), str) and getattr(compound, field) for field in ("id", "name", "cas", "formula")):
            raise ValidationError("compound identifiers must be non-empty strings")
        for key, unit in _UNITS.items():
            prop = getattr(compound, key)
            if prop.unit != unit:
                raise ValidationError(f"{compound.id}.{key} must use {unit}")
            if isinstance(prop.value, bool) or not isinstance(prop.value, (int, float)) or not math.isfinite(prop.value):
                raise ValidationError(f"{compound.id}.{key} must be finite numeric data")
            if key != "acentric_factor" and prop.value <= 0:
                raise ValidationError(f"{compound.id}.{key} must be positive")
    return compounds


def load_pr_interactions(path: str | Path) -> PRInteractions:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if data["schema_version"] != "pr-interactions-1":
            raise ValidationError("unsupported PR interaction schema")
        if data["model"] != "Peng-Robinson":
            raise ValidationError("interaction model must be Peng-Robinson")
        if data["missing_pair_policy"] not in {"error", "zero"}:
            raise ValidationError("missing_pair_policy must be error or zero")
        if any(pair["unit"] != "dimensionless" for pair in data["pairs"]):
            raise ValidationError("PR interaction kij must use dimensionless units")
        pairs = tuple((p["compound_1"], p["compound_2"], p["kij"]) for p in data["pairs"])
        default_zero = data["missing_pair_policy"] == "zero"
        provenance = InteractionProvenance(**data["provenance"])
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid PR interaction data: {error}") from error
    if len({frozenset(pair[:2]) for pair in pairs}) != len(pairs):
        raise ValidationError("PR interaction pairs must be unique")
    _validate_provenance(provenance)
    for first, second, value in pairs:
        if not isinstance(first, str) or not first or not isinstance(second, str) or not second or first == second:
            raise ValidationError("PR interaction compound IDs must be non-empty and distinct")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValidationError("PR interaction kij must be finite numeric data")
    return PRInteractions(pairs, default_zero, provenance)
