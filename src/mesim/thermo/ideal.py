import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from ..errors import OutOfRangeError, ValidationError


R = 8314.46261815324  # J/kmol/K


@dataclass(frozen=True)
class Result:
    value: float
    unit: str
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class Correlation:
    equation: int
    coefficients: tuple[float, float, float, float, float]
    minimum_k: float
    maximum_k: float
    unit: str


@dataclass(frozen=True)
class Provenance:
    source: str
    source_revision: str
    imported_utc: str


@dataclass(frozen=True)
class IdealCorrelations:
    compound_id: str
    heat_capacity_correlation: Correlation
    vapor_pressure_correlation: Correlation
    provenance: Provenance

    def _range(self, name: str, correlation: Correlation, temperatures: tuple[float, ...], allow: bool) -> tuple[str, ...]:
        if any(not math.isfinite(temperature) or temperature <= 0 for temperature in temperatures):
            raise ValidationError("absolute temperatures must be finite and positive")
        if all(correlation.minimum_k <= temperature <= correlation.maximum_k for temperature in temperatures):
            return ()
        message = f"{name} extrapolated outside {correlation.minimum_k:g}..{correlation.maximum_k:g} K"
        if not allow:
            raise OutOfRangeError(message)
        return (message,)

    def heat_capacity(self, temperature_k: float, allow_extrapolation: bool = False) -> Result:
        warnings = self._range("heat_capacity", self.heat_capacity_correlation, (temperature_k,), allow_extrapolation)
        a, b, c, d, e = self.heat_capacity_correlation.coefficients
        return Result(a + math.exp(b / temperature_k + c + d * temperature_k + e * temperature_k**2), "J/kmol/K", warnings)

    def vapor_pressure(self, temperature_k: float, allow_extrapolation: bool = False) -> Result:
        warnings = self._range("vapor_pressure", self.vapor_pressure_correlation, (temperature_k,), allow_extrapolation)
        a, b, c, d, e = self.vapor_pressure_correlation.coefficients
        return Result(math.exp(a + b / temperature_k + c * math.log(temperature_k) + d * temperature_k**e), "Pa", warnings)

    def enthalpy_change(self, temperature_k: float, reference_k: float = 298.15, allow_extrapolation: bool = False) -> Result:
        warnings = self._range("heat_capacity", self.heat_capacity_correlation, (temperature_k, reference_k), allow_extrapolation)
        return Result(_simpson(lambda value: self.heat_capacity(value, True).value, reference_k, temperature_k), "J/kmol", warnings)

    def entropy_change(self, temperature_k: float, pressure_pa: float, reference_k: float = 298.15, reference_pressure_pa: float = 101325, allow_extrapolation: bool = False) -> Result:
        if not all(math.isfinite(value) and value > 0 for value in (pressure_pa, reference_pressure_pa)):
            raise ValidationError("absolute pressures must be finite and positive")
        warnings = self._range("heat_capacity", self.heat_capacity_correlation, (temperature_k, reference_k), allow_extrapolation)
        thermal = _simpson(lambda value: self.heat_capacity(value, True).value / value, reference_k, temperature_k)
        return Result(thermal - R * math.log(pressure_pa / reference_pressure_pa), "J/kmol/K", warnings)


def _simpson(function, start: float, end: float) -> float:
    if start == end:
        return 0.0
    middle = (start + end) / 2
    first, center, last = function(start), function(middle), function(end)
    whole = (end - start) * (first + 4 * center + last) / 6

    def refine(left, right, f_left, f_middle, f_right, estimate, tolerance, depth):
        middle = (left + right) / 2
        left_middle, right_middle = (left + middle) / 2, (middle + right) / 2
        f_left_middle, f_right_middle = function(left_middle), function(right_middle)
        left_estimate = (middle - left) * (f_left + 4 * f_left_middle + f_middle) / 6
        right_estimate = (right - middle) * (f_middle + 4 * f_right_middle + f_right) / 6
        error = left_estimate + right_estimate - estimate
        if depth == 0:
            raise ValidationError("ideal-property integration did not converge")
        if abs(error) <= 15 * tolerance:
            return left_estimate + right_estimate + error / 15
        return refine(left, middle, f_left, f_left_middle, f_middle, left_estimate, tolerance / 2, depth - 1) + refine(middle, right, f_middle, f_right_middle, f_right, right_estimate, tolerance / 2, depth - 1)

    return refine(start, end, first, center, last, whole, max(1e-8, abs(whole) * 1e-12), 20)


def ideal_gas_density(molecular_weight_kg_per_kmol: float, temperature_k: float, pressure_pa: float) -> Result:
    if not all(math.isfinite(value) for value in (molecular_weight_kg_per_kmol, temperature_k, pressure_pa)):
        raise ValidationError("ideal-gas inputs must be finite")
    if molecular_weight_kg_per_kmol <= 0 or temperature_k <= 0 or pressure_pa <= 0:
        raise ValidationError("molecular weight, absolute temperature, and absolute pressure must be positive")
    return Result(pressure_pa * molecular_weight_kg_per_kmol / (R * temperature_k), "kg/m3")


def load_correlations(path: str | Path) -> tuple[IdealCorrelations, ...]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if data["schema_version"] != "ideal-correlations-1":
            raise ValidationError("unsupported ideal correlation schema")
        provenance = Provenance(**data["provenance"])
        imported = datetime.fromisoformat(provenance.imported_utc)
        if not all(vars(provenance).values()) or not provenance.imported_utc.endswith("Z") or imported.utcoffset() != timedelta(0):
            raise ValidationError("correlation provenance must be non-empty and use a UTC import timestamp")
        records = tuple(
            IdealCorrelations(
                compound_id=record["compound_id"],
                heat_capacity_correlation=_correlation(record["heat_capacity"], 16, "J/kmol/K"),
                vapor_pressure_correlation=_correlation(record["vapor_pressure"], 101, "Pa"),
                provenance=provenance,
            )
            for record in data["correlations"]
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid ideal correlation data: {error}") from error
    if len({record.compound_id for record in records}) != len(records):
        raise ValidationError("ideal correlation compound IDs must be unique")
    return records


def _correlation(data: dict, equation: int, unit: str) -> Correlation:
    if data["equation"] != equation or data["unit"] != unit:
        raise ValidationError(f"expected equation {equation} in {unit}")
    coefficients = tuple(data[key] for key in "ABCDE")
    values = coefficients + (data["minimum_k"], data["maximum_k"])
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in values):
        raise ValidationError("correlation values must be finite numbers")
    if data["minimum_k"] <= 0 or data["maximum_k"] <= data["minimum_k"]:
        raise ValidationError("invalid correlation temperature range")
    return Correlation(equation, coefficients, data["minimum_k"], data["maximum_k"], unit)
