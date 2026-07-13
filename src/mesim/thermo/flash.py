import math
from dataclasses import dataclass

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError
from .peng_robinson import PRMixtureState, PengRobinsonMixture


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


@dataclass(frozen=True, slots=True)
class StabilityTrial:
    report: SolverReport
    kind: str
    tangent_plane_distance: float | None
    composition: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class StabilityResult:
    report: SolverReport
    stable: bool
    feed_phase: str
    vapor_like: StabilityTrial
    liquid_like: StabilityTrial


@dataclass(frozen=True, slots=True)
class TPFlashResult:
    report: SolverReport
    stability: StabilityResult
    temperature_k: float
    pressure_pa: float
    phase: str
    vapor_fraction: float | None
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    liquid_state: PRMixtureState | None
    vapor_state: PRMixtureState | None


def _report(
    converged: bool,
    iterations: int,
    residual: float,
    failure_reason: str | None = None,
    algorithm: str = "Rachford-Rice bisection",
    warnings: tuple[str, ...] = (),
) -> SolverReport:
    return SolverReport(converged, iterations, residual, algorithm, warnings, failure_reason)


def _validate_solver_limits(max_iterations: int, tolerance: float) -> None:
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValidationError("max_iterations must be a positive integer")
    if isinstance(tolerance, bool) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValidationError("tolerance must be finite and positive")


def _validate(composition: tuple[float, ...], k_values: tuple[float, ...], max_iterations: int, tolerance: float) -> None:
    if not composition or len(composition) != len(k_values):
        raise ValidationError("composition and K values must have the same nonzero length")
    if any(isinstance(value, bool) or not math.isfinite(value) or value < 0 for value in composition):
        raise ValidationError("composition must be finite and non-negative")
    if not math.isclose(math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("composition must sum to one")
    if any(isinstance(value, bool) or not math.isfinite(value) or value <= 0 for value in k_values):
        raise ValidationError("K values must be finite and positive")
    _validate_solver_limits(max_iterations, tolerance)


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


def _wilson_log_k(compounds: tuple[Compound, ...], temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
    return tuple(
        math.log(compound.critical_pressure.value / pressure_pa)
        + 5.373
        * (1.0 + compound.acentric_factor.value)
        * (1.0 - compound.critical_temperature.value / temperature_k)
        for compound in compounds
    )


def _normalized_log_weights(values: tuple[float, ...]) -> tuple[float, ...]:
    maximum = max(values)
    weights = tuple(0.0 if value == -math.inf else math.exp(value - maximum) for value in values)
    total = math.fsum(weights)
    return tuple(value / total for value in weights)


def _stability_trial(
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    feed_fugacity: tuple[float, ...],
    kind: str,
    max_iterations: int,
    tolerance: float,
) -> StabilityTrial:
    log_k = _wilson_log_k(compounds, temperature_k, pressure_pa)
    sign = 1.0 if kind == "vapor-like" else -1.0
    trial = _normalized_log_weights(
        tuple(-math.inf if fraction == 0.0 else math.log(fraction) + sign * value for fraction, value in zip(composition, log_k))
    )
    phase = "vapor" if kind == "vapor-like" else "liquid"
    residual = math.inf

    for iteration in range(1, max_iterations + 1):
        trial_fugacity = PengRobinsonMixture(compounds, trial, interactions).state(
            temperature_k, pressure_pa, phase
        ).fugacity_coefficients
        updated = _normalized_log_weights(
            tuple(
                -math.inf
                if fraction == 0.0
                else math.log(fraction) + math.log(feed_phi) - math.log(trial_phi)
                for fraction, feed_phi, trial_phi in zip(composition, feed_fugacity, trial_fugacity)
            )
        )
        residual = max(abs(first - second) for first, second in zip(trial, updated))
        trial = updated
        if residual <= tolerance:
            trial_fugacity = PengRobinsonMixture(compounds, trial, interactions).state(
                temperature_k, pressure_pa, phase
            ).fugacity_coefficients
            tangent_plane_distance = math.fsum(
                fraction
                * (math.log(fraction) + math.log(trial_phi) - math.log(feed_fraction) - math.log(feed_phi))
                for fraction, trial_phi, feed_fraction, feed_phi in zip(
                    trial, trial_fugacity, composition, feed_fugacity
                )
                if fraction > 0.0 and feed_fraction > 0.0
            )
            report = _report(True, iteration, residual, algorithm="Michelsen tangent-plane iteration")
            return StabilityTrial(report, kind, tangent_plane_distance, trial)

    reason = f"{kind} stability iteration limit reached"
    report = _report(False, max_iterations, residual, reason, "Michelsen tangent-plane iteration")
    return StabilityTrial(report, kind, None, trial)


def pr_stability(
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    *,
    max_iterations: int = 100,
    tolerance: float = 1e-10,
) -> StabilityResult:
    _validate_solver_limits(max_iterations, tolerance)
    feed = PengRobinsonMixture(compounds, composition, interactions).stable_state(temperature_k, pressure_pa)
    vapor = _stability_trial(
        compounds,
        composition,
        interactions,
        temperature_k,
        pressure_pa,
        feed.fugacity_coefficients,
        "vapor-like",
        max_iterations,
        tolerance,
    )
    liquid = _stability_trial(
        compounds,
        composition,
        interactions,
        temperature_k,
        pressure_pa,
        feed.fugacity_coefficients,
        "liquid-like",
        max_iterations,
        tolerance,
    )
    converged = vapor.report.converged and liquid.report.converged
    failure_reason = None if converged else "one or more stability trials did not converge"
    residual = max(vapor.report.residual, liquid.report.residual)
    report = _report(
        converged,
        max(vapor.report.iterations, liquid.report.iterations),
        residual,
        failure_reason,
        "Michelsen tangent-plane distance",
    )
    stable = converged and min(vapor.tangent_plane_distance, liquid.tangent_plane_distance) >= -tolerance
    return StabilityResult(report, stable, feed.phase, vapor, liquid)


def _tp_failure(
    stability: StabilityResult,
    temperature_k: float,
    pressure_pa: float,
    iterations: int,
    residual: float,
    reason: str,
    vapor_fraction: float | None = None,
    liquid_composition: tuple[float, ...] = (),
    vapor_composition: tuple[float, ...] = (),
    liquid_state: PRMixtureState | None = None,
    vapor_state: PRMixtureState | None = None,
) -> TPFlashResult:
    report = _report(False, iterations, residual, reason, "PR TP successive substitution")
    return TPFlashResult(
        report,
        stability,
        temperature_k,
        pressure_pa,
        "unknown" if vapor_fraction is None else "two-phase",
        vapor_fraction,
        liquid_composition,
        vapor_composition,
        liquid_state,
        vapor_state,
    )


def tp_flash(
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    *,
    max_iterations: int = 100,
) -> TPFlashResult:
    _validate_solver_limits(max_iterations, 1e-10)
    stability = pr_stability(
        compounds,
        composition,
        interactions,
        temperature_k,
        pressure_pa,
        max_iterations=max_iterations,
    )
    if not stability.report.converged:
        return _tp_failure(
            stability,
            temperature_k,
            pressure_pa,
            stability.report.iterations,
            stability.report.residual,
            "phase-stability calculation did not converge",
        )

    feed_model = PengRobinsonMixture(compounds, composition, interactions)
    if stability.stable:
        state = feed_model.stable_state(temperature_k, pressure_pa)
        phase = stability.feed_phase
        report = _report(True, 0, stability.report.residual, algorithm="PR TP stability result")
        if phase == "liquid":
            return TPFlashResult(report, stability, temperature_k, pressure_pa, phase, 0.0, composition, (), state, None)
        if phase == "vapor":
            return TPFlashResult(report, stability, temperature_k, pressure_pa, phase, 1.0, (), composition, None, state)
        return TPFlashResult(report, stability, temperature_k, pressure_pa, "single", None, composition, composition, None, state)

    log_k = _wilson_log_k(compounds, temperature_k, pressure_pa)
    previous_residual = math.inf
    damping = 1.0
    last_rr = None
    liquid_state = None
    vapor_state = None
    fugacity_residual = math.inf
    material_residual = math.inf

    for iteration in range(1, max_iterations + 1):
        try:
            k_values = tuple(math.exp(value) for value in log_k)
        except OverflowError:
            return _tp_failure(stability, temperature_k, pressure_pa, iteration, math.inf, "K values are outside the representable range")
        last_rr = rachford_rice(composition, k_values)
        if not last_rr.report.converged or last_rr.phase != "two-phase":
            return _tp_failure(
                stability,
                temperature_k,
                pressure_pa,
                iteration,
                abs(last_rr.report.residual),
                "unstable feed did not produce a bounded two-phase split",
                last_rr.vapor_fraction,
                last_rr.liquid_composition,
                last_rr.vapor_composition,
            )
        liquid_state = PengRobinsonMixture(compounds, last_rr.liquid_composition, interactions).state(
            temperature_k, pressure_pa, "liquid"
        )
        vapor_state = PengRobinsonMixture(compounds, last_rr.vapor_composition, interactions).state(
            temperature_k, pressure_pa, "vapor"
        )
        fugacity_terms = tuple(
            math.log(liquid_fraction)
            + math.log(liquid_phi)
            - math.log(vapor_fraction)
            - math.log(vapor_phi)
            for overall, liquid_fraction, vapor_fraction, liquid_phi, vapor_phi in zip(
                composition,
                last_rr.liquid_composition,
                last_rr.vapor_composition,
                liquid_state.fugacity_coefficients,
                vapor_state.fugacity_coefficients,
            )
            if overall > 0.0
        )
        fugacity_residual = max(abs(value) for value in fugacity_terms)
        material_residual = max(
            abs(
                overall
                - (1.0 - last_rr.vapor_fraction) * liquid_fraction
                - last_rr.vapor_fraction * vapor_fraction
            )
            for overall, liquid_fraction, vapor_fraction in zip(
                composition, last_rr.liquid_composition, last_rr.vapor_composition
            )
        )
        residual = max(fugacity_residual, material_residual)
        if fugacity_residual <= 1e-8 and material_residual <= 1e-10:
            warnings = () if damping == 1.0 else (f"K-value damping reduced to {damping}",)
            report = _report(True, iteration, residual, algorithm="PR TP successive substitution", warnings=warnings)
            return TPFlashResult(
                report,
                stability,
                temperature_k,
                pressure_pa,
                "two-phase",
                last_rr.vapor_fraction,
                last_rr.liquid_composition,
                last_rr.vapor_composition,
                liquid_state,
                vapor_state,
            )
        if residual > previous_residual:
            damping = max(1.0 / 16.0, damping / 2.0)
        target_log_k = tuple(
            math.log(liquid_phi) - math.log(vapor_phi)
            for liquid_phi, vapor_phi in zip(
                liquid_state.fugacity_coefficients, vapor_state.fugacity_coefficients
            )
        )
        log_k = tuple(
            current + damping * (target - current) for current, target in zip(log_k, target_log_k)
        )
        previous_residual = residual

    return _tp_failure(
        stability,
        temperature_k,
        pressure_pa,
        max_iterations,
        max(fugacity_residual, material_residual),
        "PR TP iteration limit reached",
        last_rr.vapor_fraction,
        last_rr.liquid_composition,
        last_rr.vapor_composition,
        liquid_state,
        vapor_state,
    )
