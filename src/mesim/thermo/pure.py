"""Source-backed saturated pure-component property correlations."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from ..errors import OutOfRangeError, ValidationError
from .correlations import evaluate_temperature_equation
from .ideal import Provenance, Result


@dataclass(frozen=True, slots=True)
class TemperatureCorrelation:
    equation: int
    coefficients: tuple[float, float, float, float, float]
    minimum_k: float
    maximum_k: float
    unit: str

    def evaluate(
        self,
        temperature_k: float,
        *,
        critical_temperature_k: float,
        allow_extrapolation: bool = False,
    ) -> Result:
        if not isinstance(allow_extrapolation, bool):
            raise ValidationError("correlation extrapolation flag must be boolean")
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        warnings = ()
        if not self.minimum_k <= temperature_k <= self.maximum_k:
            message = (
                f"correlation extrapolated outside "
                f"{self.minimum_k:g}..{self.maximum_k:g} K"
            )
            if not allow_extrapolation:
                raise OutOfRangeError(message)
            warnings = (message,)
        value = evaluate_temperature_equation(
            self.equation,
            self.coefficients,
            temperature_k,
            critical_temperature_k,
        )
        if not math.isfinite(value) or value <= 0:
            raise ValidationError("pure-component correlation must be positive and finite")
        return Result(value, self.unit, warnings)


@dataclass(frozen=True, slots=True)
class SaturatedLiquidCorrelations:
    compound_id: str
    critical_temperature_k: float
    liquid_density_correlation: TemperatureCorrelation
    liquid_heat_capacity_correlation: TemperatureCorrelation
    heat_of_vaporization_correlation: TemperatureCorrelation
    surface_tension_correlation: TemperatureCorrelation
    provenance: Provenance

    def _evaluate(
        self,
        correlation: TemperatureCorrelation,
        temperature_k: float,
        allow_extrapolation: bool,
    ) -> Result:
        return correlation.evaluate(
            temperature_k,
            critical_temperature_k=self.critical_temperature_k,
            allow_extrapolation=allow_extrapolation,
        )

    def liquid_molar_density(
        self, temperature_k: float, allow_extrapolation: bool = False
    ) -> Result:
        return self._evaluate(
            self.liquid_density_correlation, temperature_k, allow_extrapolation
        )

    def liquid_heat_capacity(
        self, temperature_k: float, allow_extrapolation: bool = False
    ) -> Result:
        return self._evaluate(
            self.liquid_heat_capacity_correlation, temperature_k, allow_extrapolation
        )

    def heat_of_vaporization(
        self, temperature_k: float, allow_extrapolation: bool = False
    ) -> Result:
        return self._evaluate(
            self.heat_of_vaporization_correlation, temperature_k, allow_extrapolation
        )

    def surface_tension(
        self, temperature_k: float, allow_extrapolation: bool = False
    ) -> Result:
        return self._evaluate(
            self.surface_tension_correlation, temperature_k, allow_extrapolation
        )


def _correlation(data: dict, equations: tuple[int, ...], unit: str) -> TemperatureCorrelation:
    coefficients = tuple(data[key] for key in "ABCDE")
    minimum, maximum = data["minimum_k"], data["maximum_k"]
    values = (*coefficients, minimum, maximum)
    if (
        data["equation"] not in equations
        or data["unit"] != unit
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            for value in values
        )
        or minimum < 0
        or maximum <= minimum
    ):
        raise ValidationError(f"invalid saturated-property correlation in {unit}")
    return TemperatureCorrelation(
        data["equation"], coefficients, minimum, maximum, unit
    )


def load_saturated_liquid_correlations(
    path: str | Path,
) -> tuple[SaturatedLiquidCorrelations, ...]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if data["schema_version"] != "saturated-liquid-correlations-1":
            raise ValidationError("unsupported saturated-liquid correlation schema")
        provenance = Provenance(**data["provenance"])
        imported = datetime.fromisoformat(provenance.imported_utc)
        if (
            not all(vars(provenance).values())
            or not provenance.imported_utc.endswith("Z")
            or imported.utcoffset() != timedelta(0)
        ):
            raise ValidationError(
                "correlation provenance must be non-empty and use a UTC import timestamp"
            )
        records = tuple(
            SaturatedLiquidCorrelations(
                compound_id=record["compound_id"],
                critical_temperature_k=_critical_temperature(
                    record["critical_temperature"]
                ),
                liquid_density_correlation=_correlation(
                    record["liquid_density"], (105, 106), "kmol/m3"
                ),
                liquid_heat_capacity_correlation=_correlation(
                    record["liquid_heat_capacity"], (3, 4, 16, 100), "J/kmol/K"
                ),
                heat_of_vaporization_correlation=_correlation(
                    record["heat_of_vaporization"], (106,), "J/kmol"
                ),
                surface_tension_correlation=_correlation(
                    record["surface_tension"], (2, 16, 106, 116), "N/m"
                ),
                provenance=provenance,
            )
            for record in data["correlations"]
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid saturated-liquid correlation data: {error}") from error
    if (
        not records
        or len({record.compound_id for record in records}) != len(records)
        or any(not isinstance(record.compound_id, str) or not record.compound_id for record in records)
    ):
        raise ValidationError(
            "saturated-liquid correlation compound IDs must be non-empty and unique"
        )
    return records


def _critical_temperature(data: dict) -> float:
    value = data["value"]
    if (
        data["unit"] != "K"
        or isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ValidationError("critical temperature must be a positive finite K value")
    return value
