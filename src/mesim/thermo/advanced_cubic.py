"""Temperature-dependent interaction layer for DWSIM advanced PR/SRK packages."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..compounds import PRInteractions
from ..errors import ValidationError


@dataclass(frozen=True, slots=True)
class AdvancedCubicInteraction:
    first: str
    second: str
    expression: str
    temperature_coefficients: tuple[float, ...]

    def value(self, temperature_k: float) -> float:
        if (
            isinstance(temperature_k, bool)
            or not isinstance(temperature_k, (int, float))
            or not math.isfinite(temperature_k)
            or temperature_k <= 0.0
        ):
            raise ValidationError("advanced cubic temperature must be positive")
        result = 0.0
        for coefficient in reversed(self.temperature_coefficients):
            result = result * temperature_k + coefficient
        if not math.isfinite(result):
            raise ValidationError("advanced cubic interaction is unrepresentable")
        return result


@dataclass(frozen=True, slots=True)
class AdvancedCubicData:
    source_revision: str
    interactions: tuple[AdvancedCubicInteraction, ...]

    def interaction(
        self, first: str, second: str
    ) -> AdvancedCubicInteraction | None:
        return next(
            (
                record
                for record in self.interactions
                if (record.first == first and record.second == second)
                or (record.first == second and record.second == first)
            ),
            None,
        )

    def evaluated(
        self,
        temperature_k: float,
        pressure_pa: float,
        base: PRInteractions,
    ) -> "EvaluatedAdvancedCubicInteractions":
        if (
            isinstance(pressure_pa, bool)
            or not isinstance(pressure_pa, (int, float))
            or not math.isfinite(pressure_pa)
            or pressure_pa <= 0.0
            or not isinstance(base, PRInteractions)
        ):
            raise ValidationError("advanced cubic state or base interactions are invalid")
        # The current DWSIM source table depends only on T. Pressure is retained
        # in this boundary because user-entered package expressions may use P.
        return EvaluatedAdvancedCubicInteractions(
            self, float(temperature_k), float(pressure_pa), base
        )


@dataclass(frozen=True, slots=True)
class EvaluatedAdvancedCubicInteractions:
    data: AdvancedCubicData
    temperature_k: float
    pressure_pa: float
    base: PRInteractions

    @property
    def model(self) -> str:
        return self.base.model

    def get(self, first: str, second: str) -> float:
        if first == second:
            return 0.0
        record = self.data.interaction(first, second)
        value = None if record is None else record.value(self.temperature_k)
        # DWSIM falls back to its ordinary cubic table when an advanced
        # expression is absent or evaluates to exactly zero.
        if value is not None and value != 0.0:
            return value
        try:
            return self.base.get(first, second)
        except ValidationError:
            # The DWSIM package's internal table uses zero for an unlisted pair.
            return 0.0


def load_advanced_cubic_data(path: str | Path) -> AdvancedCubicData:
    """Load and validate the complete DWSIM mercury interaction table."""
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if document["schema_version"] != "dwsim-prsrk-advanced-interactions-1":
            raise ValidationError("unsupported advanced cubic interaction schema")
        source_revision = document["source"]["revision"]
        records = tuple(
            AdvancedCubicInteraction(
                record["first"],
                record["second"],
                record["expression"],
                tuple(record["temperature_coefficients"]),
            )
            for record in document["pairs"]
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid advanced cubic interaction data: {error}") from error
    if (
        not isinstance(source_revision, str)
        or not source_revision
        or len(records) != 13
    ):
        raise ValidationError("advanced cubic source identity or record count is invalid")
    pairs = tuple(frozenset((record.first, record.second)) for record in records)
    if len(set(pairs)) != len(pairs):
        raise ValidationError("advanced cubic interaction pairs must be unique")
    for record in records:
        if (
            not isinstance(record.first, str)
            or not record.first
            or not isinstance(record.second, str)
            or not record.second
            or record.first == record.second
            or not isinstance(record.expression, str)
            or not record.expression
            or not 1 <= len(record.temperature_coefficients) <= 4
            or any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                for value in record.temperature_coefficients
            )
        ):
            raise ValidationError("advanced cubic interaction record is invalid")
    return AdvancedCubicData(source_revision, records)


__all__ = (
    "AdvancedCubicData",
    "AdvancedCubicInteraction",
    "EvaluatedAdvancedCubicInteractions",
    "load_advanced_cubic_data",
)
