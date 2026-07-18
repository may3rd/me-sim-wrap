"""Ideal liquid/vapor equilibrium using DWSIM's Raoult-law contract."""

import math
from dataclasses import dataclass

from ..errors import ValidationError
from .flash import SolverReport, rachford_rice
from .ideal import IdealCorrelations


@dataclass(frozen=True, slots=True)
class RaoultVLEResult:
    report: SolverReport
    kind: str
    temperature_k: float
    pressure_pa: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    equilibrium_ratios: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class RaoultTPFlashResult:
    report: SolverReport
    temperature_k: float
    pressure_pa: float
    phase: str
    vapor_fraction: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    equilibrium_ratios: tuple[float, ...]


def _validate_domain(
    correlations: tuple[IdealCorrelations, ...],
    composition: tuple[float, ...],
) -> None:
    if (
        not correlations
        or len(correlations) != len(composition)
        or any(not isinstance(record, IdealCorrelations) for record in correlations)
    ):
        raise ValidationError(
            "Raoult correlations and composition must have the same nonzero length"
        )
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
        for value in composition
    ):
        raise ValidationError("Raoult composition must be finite and non-negative")
    if not math.isclose(math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValidationError("Raoult composition must sum to one")


def _pressure(pressure_pa: float) -> float:
    if (
        isinstance(pressure_pa, bool)
        or not isinstance(pressure_pa, (int, float))
        or not math.isfinite(pressure_pa)
        or pressure_pa <= 0
    ):
        raise ValidationError("absolute pressure must be finite and positive")
    return float(pressure_pa)


def _vapor_pressures(
    correlations: tuple[IdealCorrelations, ...],
    temperature_k: float,
    *,
    allow_extrapolation: bool,
) -> tuple[tuple[float, ...], tuple[str, ...]]:
    if not isinstance(allow_extrapolation, bool):
        raise ValidationError("Raoult extrapolation flag must be boolean")
    results = tuple(
        record.vapor_pressure(temperature_k, allow_extrapolation)
        for record in correlations
    )
    warnings = tuple(dict.fromkeys(warning for result in results for warning in result.warnings))
    return tuple(result.value for result in results), warnings


def raoult_fugacity_coefficients(
    correlations: tuple[IdealCorrelations, ...],
    temperature_k: float,
    pressure_pa: float,
    phase: str,
    *,
    allow_extrapolation: bool = False,
) -> tuple[float, ...]:
    """Return Psat/P for liquid or unity for vapor, as DWSIM Ideal.vb does."""
    pressure = _pressure(pressure_pa)
    if phase == "vapor":
        return (1.0,) * len(correlations)
    if phase != "liquid":
        raise ValidationError("Raoult phase must be liquid or vapor")
    vapor_pressures, _ = _vapor_pressures(
        correlations, temperature_k, allow_extrapolation=allow_extrapolation
    )
    return tuple(value / pressure for value in vapor_pressures)


def raoult_equilibrium_ratios(
    correlations: tuple[IdealCorrelations, ...],
    temperature_k: float,
    pressure_pa: float,
    *,
    allow_extrapolation: bool = False,
) -> tuple[float, ...]:
    return raoult_fugacity_coefficients(
        correlations,
        temperature_k,
        pressure_pa,
        "liquid",
        allow_extrapolation=allow_extrapolation,
    )


def raoult_bubble_pressure(
    correlations: tuple[IdealCorrelations, ...],
    liquid_composition: tuple[float, ...],
    temperature_k: float,
    *,
    allow_extrapolation: bool = False,
) -> RaoultVLEResult:
    _validate_domain(correlations, liquid_composition)
    vapor_pressures, warnings = _vapor_pressures(
        correlations, temperature_k, allow_extrapolation=allow_extrapolation
    )
    pressure = math.fsum(
        fraction * vapor_pressure
        for fraction, vapor_pressure in zip(liquid_composition, vapor_pressures)
    )
    _pressure(pressure)
    vapor = tuple(
        fraction * vapor_pressure / pressure
        for fraction, vapor_pressure in zip(liquid_composition, vapor_pressures)
    )
    ratios = tuple(value / pressure for value in vapor_pressures)
    residual = abs(math.fsum(vapor) - 1.0)
    report = SolverReport(
        True, 0, residual, "Raoult bubble-pressure direct sum", warnings, None
    )
    return RaoultVLEResult(
        report,
        "bubble",
        temperature_k,
        pressure,
        liquid_composition,
        vapor,
        ratios,
    )


def raoult_dew_pressure(
    correlations: tuple[IdealCorrelations, ...],
    vapor_composition: tuple[float, ...],
    temperature_k: float,
    *,
    allow_extrapolation: bool = False,
) -> RaoultVLEResult:
    _validate_domain(correlations, vapor_composition)
    vapor_pressures, warnings = _vapor_pressures(
        correlations, temperature_k, allow_extrapolation=allow_extrapolation
    )
    denominator = math.fsum(
        fraction / vapor_pressure
        for fraction, vapor_pressure in zip(vapor_composition, vapor_pressures)
    )
    pressure = _pressure(1.0 / denominator)
    liquid = tuple(
        fraction * pressure / vapor_pressure
        for fraction, vapor_pressure in zip(vapor_composition, vapor_pressures)
    )
    ratios = tuple(value / pressure for value in vapor_pressures)
    residual = abs(math.fsum(liquid) - 1.0)
    report = SolverReport(
        True, 0, residual, "Raoult dew-pressure reciprocal sum", warnings, None
    )
    return RaoultVLEResult(
        report,
        "dew",
        temperature_k,
        pressure,
        liquid,
        vapor_composition,
        ratios,
    )


def raoult_tp_flash(
    correlations: tuple[IdealCorrelations, ...],
    composition: tuple[float, ...],
    temperature_k: float,
    pressure_pa: float,
    *,
    max_iterations: int = 100,
    tolerance: float = 1.0e-12,
    allow_extrapolation: bool = False,
) -> RaoultTPFlashResult:
    _validate_domain(correlations, composition)
    pressure = _pressure(pressure_pa)
    vapor_pressures, warnings = _vapor_pressures(
        correlations, temperature_k, allow_extrapolation=allow_extrapolation
    )
    ratios = tuple(value / pressure for value in vapor_pressures)
    flash = rachford_rice(
        composition,
        ratios,
        max_iterations=max_iterations,
        tolerance=tolerance,
    )
    report = SolverReport(
        flash.report.converged,
        flash.report.iterations,
        flash.report.residual,
        flash.report.algorithm,
        tuple(dict.fromkeys((*flash.report.warnings, *warnings))),
        flash.report.failure_reason,
    )
    return RaoultTPFlashResult(
        report,
        temperature_k,
        pressure,
        flash.phase,
        flash.vapor_fraction,
        flash.liquid_composition,
        flash.vapor_composition,
        ratios,
    )
