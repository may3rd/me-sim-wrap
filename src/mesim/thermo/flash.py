import math
from dataclasses import dataclass

from ..errors import ValidationError


@dataclass(frozen=True, slots=True)
class SolverReport:
    converged: bool
    iterations: int
    residual: float
    algorithm: str
    warnings: tuple[str, ...]
    failure_reason: str | None


@dataclass(frozen=True, slots=True)
class FlashResult:
    report: SolverReport
    phase: str
    vapor_fraction: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]


def _report(converged: bool, iterations: int, residual: float, failure_reason: str | None = None) -> SolverReport:
    return SolverReport(converged, iterations, residual, "Rachford-Rice bisection", (), failure_reason)


def _validate(composition: tuple[float, ...], k_values: tuple[float, ...], max_iterations: int, tolerance: float) -> None:
    if not composition or len(composition) != len(k_values):
        raise ValidationError("composition and K values must have the same nonzero length")
    if any(isinstance(value, bool) or not math.isfinite(value) or value < 0 for value in composition):
        raise ValidationError("composition must be finite and non-negative")
    if not math.isclose(math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("composition must sum to one")
    if any(isinstance(value, bool) or not math.isfinite(value) or value <= 0 for value in k_values):
        raise ValidationError("K values must be finite and positive")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValidationError("max_iterations must be a positive integer")
    if isinstance(tolerance, bool) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValidationError("tolerance must be finite and positive")


def _residual(beta: float, composition: tuple[float, ...], k_values: tuple[float, ...]) -> float:
    return math.fsum(
        fraction * (k_value - 1.0) / (1.0 + beta * (k_value - 1.0))
        for fraction, k_value in zip(composition, k_values)
    )


def _compositions(beta: float, composition: tuple[float, ...], k_values: tuple[float, ...]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    liquid = tuple(
        fraction / (1.0 + beta * (k_value - 1.0))
        for fraction, k_value in zip(composition, k_values)
    )
    return liquid, tuple(k_value * fraction for k_value, fraction in zip(k_values, liquid))


def rachford_rice(
    composition: tuple[float, ...],
    k_values: tuple[float, ...],
    *,
    max_iterations: int = 100,
    tolerance: float = 1e-12,
) -> FlashResult:
    _validate(composition, k_values, max_iterations, tolerance)
    if all(math.isclose(k_value, 1.0, rel_tol=0.0, abs_tol=tolerance) for k_value in k_values):
        reason = "Rachford-Rice is indeterminate when all K values equal 1"
        return FlashResult(_report(False, 0, 0.0, reason), "indeterminate", 0.5, composition, composition)

    at_liquid = _residual(0.0, composition, k_values)
    if at_liquid <= 0.0:
        return FlashResult(_report(True, 0, at_liquid), "liquid", 0.0, composition, ())
    at_vapor = _residual(1.0, composition, k_values)
    if at_vapor >= 0.0:
        return FlashResult(_report(True, 0, at_vapor), "vapor", 1.0, (), composition)

    low, high = 0.0, 1.0
    beta, residual = 0.5, _residual(0.5, composition, k_values)
    for iteration in range(1, max_iterations + 1):
        beta = (low + high) / 2.0
        residual = _residual(beta, composition, k_values)
        if residual > 0.0:
            low = beta
        else:
            high = beta
        if high - low <= tolerance and abs(residual) <= tolerance:
            liquid, vapor = _compositions(beta, composition, k_values)
            if not math.isclose(math.fsum(liquid), 1.0, rel_tol=0.0, abs_tol=1e-12) or not math.isclose(math.fsum(vapor), 1.0, rel_tol=0.0, abs_tol=1e-12):
                reason = "phase compositions do not sum to one"
                return FlashResult(_report(False, iteration, residual, reason), "two-phase", beta, liquid, vapor)
            return FlashResult(_report(True, iteration, residual), "two-phase", beta, liquid, vapor)

    liquid, vapor = _compositions(beta, composition, k_values)
    reason = "Rachford-Rice iteration limit reached"
    return FlashResult(_report(False, max_iterations, residual, reason), "two-phase", beta, liquid, vapor)
