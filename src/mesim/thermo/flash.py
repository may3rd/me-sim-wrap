import math
from dataclasses import dataclass

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError
from .ideal import IdealCorrelations
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


@dataclass(frozen=True, slots=True)
class PressureResult:
    report: SolverReport
    kind: str
    temperature_k: float
    pressure_pa: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    liquid_state: PRMixtureState | None
    vapor_state: PRMixtureState | None


@dataclass(frozen=True, slots=True)
class PHFlashResult:
    report: SolverReport
    temperature_k: float
    pressure_pa: float
    target_enthalpy_j_per_kmol: float
    enthalpy_j_per_kmol: float | None
    temperature_bracket_k: tuple[float, float]
    outer_iterations: int
    inner_iterations: int
    flash: TPFlashResult | None


@dataclass(frozen=True, slots=True)
class PSFlashResult:
    report: SolverReport
    temperature_k: float
    pressure_pa: float
    target_entropy_j_per_kmol_k: float
    entropy_j_per_kmol_k: float | None
    temperature_bracket_k: tuple[float, float]
    outer_iterations: int
    inner_iterations: int
    flash: TPFlashResult | None


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
    if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValidationError("tolerance must be finite and positive")


def _validate(composition: tuple[float, ...], k_values: tuple[float, ...], max_iterations: int, tolerance: float) -> None:
    if not composition or len(composition) != len(k_values):
        raise ValidationError("composition and K values must have the same nonzero length")
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 for value in composition):
        raise ValidationError("composition must be finite and non-negative")
    if not math.isclose(math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("composition must sum to one")
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0 for value in k_values):
        raise ValidationError("K values must be finite and positive")
    _validate_solver_limits(max_iterations, tolerance)


def _residual(beta: float, composition: tuple[float, ...], k_values: tuple[float, ...]) -> float:
    return math.fsum(
        fraction * (k_value - 1.0) / ((1.0 - beta) + beta * k_value)
        for fraction, k_value in zip(composition, k_values)
    )


def _compositions(beta: float, composition: tuple[float, ...], k_values: tuple[float, ...]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    liquid = tuple(
        fraction / ((1.0 - beta) + beta * k_value)
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
    _validate_absolute_pressure(pressure_pa)
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


def _validate_pressure_bracket(bracket_pa: tuple[float, float]) -> None:
    if len(bracket_pa) != 2 or any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0
        for value in bracket_pa
    ) or bracket_pa[0] >= bracket_pa[1]:
        raise ValidationError("pressure bracket must contain two increasing finite positive values")


def _equilibrium_composition(
    composition: tuple[float, ...], log_k: tuple[float, ...], sign: float
) -> tuple[tuple[float, ...], float]:
    terms = tuple(-math.inf if fraction == 0.0 else math.log(fraction) + sign * value for fraction, value in zip(composition, log_k))
    maximum = max(terms)
    weights = tuple(0.0 if value == -math.inf else math.exp(value - maximum) for value in terms)
    total = math.fsum(weights)
    return tuple(value / total for value in weights), maximum + math.log(total)


def _pressure_point(
    kind: str,
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    max_iterations: int,
) -> PressureResult:
    liquid_fixed = kind == "bubble"
    fixed_phase = "liquid" if liquid_fixed else "vapor"
    fixed_state = PengRobinsonMixture(compounds, composition, interactions).state(temperature_k, pressure_pa, fixed_phase)
    log_k = _wilson_log_k(compounds, temperature_k, pressure_pa)
    last_composition: tuple[float, ...] = ()
    last_state: PRMixtureState | None = None
    residual = math.inf

    for iteration in range(1, max_iterations + 1):
        try:
            last_composition, log_total = _equilibrium_composition(composition, log_k, 1.0 if liquid_fixed else -1.0)
            last_state = PengRobinsonMixture(compounds, last_composition, interactions).state(
                temperature_k, pressure_pa, "vapor" if liquid_fixed else "liquid"
            )
        except (OverflowError, ValueError, ValidationError) as error:
            report = _report(False, iteration, math.inf, f"{kind} pressure K iteration failed: {error}", f"PR {kind} K iteration")
            return PressureResult(report, kind, temperature_k, pressure_pa, composition if liquid_fixed else (), () if liquid_fixed else composition, fixed_state if liquid_fixed else last_state, last_state if liquid_fixed else fixed_state)
        target_log_k = tuple(
            math.log(liquid_phi) - math.log(vapor_phi)
            for liquid_phi, vapor_phi in zip(
                fixed_state.fugacity_coefficients if liquid_fixed else last_state.fugacity_coefficients,
                last_state.fugacity_coefficients if liquid_fixed else fixed_state.fugacity_coefficients,
            )
        )
        residual = max(abs(current - target) for current, target in zip(log_k, target_log_k))
        if residual <= 1e-10:
            if math.isclose(fixed_state.compressibility, last_state.compressibility, rel_tol=1e-9, abs_tol=1e-12):
                report = _report(False, iteration, 0.0, f"{kind} pressure requires distinct liquid and vapor roots", f"PR {kind} K iteration")
                return PressureResult(report, kind, temperature_k, pressure_pa, composition if liquid_fixed else last_composition, last_composition if liquid_fixed else composition, fixed_state if liquid_fixed else last_state, last_state if liquid_fixed else fixed_state)
            try:
                outer_residual = math.expm1(log_total)
            except OverflowError:
                outer_residual = math.inf
            report = _report(True, iteration, outer_residual, algorithm=f"PR {kind} K iteration")
            return PressureResult(
                report,
                kind,
                temperature_k,
                pressure_pa,
                composition if liquid_fixed else last_composition,
                last_composition if liquid_fixed else composition,
                fixed_state if liquid_fixed else last_state,
                last_state if liquid_fixed else fixed_state,
            )
        log_k = target_log_k

    report = _report(False, max_iterations, residual, f"{kind} pressure K iteration limit reached", f"PR {kind} K iteration")
    return PressureResult(report, kind, temperature_k, pressure_pa, composition if liquid_fixed else last_composition, last_composition if liquid_fixed else composition, fixed_state if liquid_fixed else last_state, last_state if liquid_fixed else fixed_state)


def _pressure_solve(
    kind: str,
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    bracket_pa: tuple[float, float],
    max_iterations: int,
) -> PressureResult:
    _validate_solver_limits(max_iterations, 1e-10)
    _validate_pressure_bracket(bracket_pa)
    PengRobinsonMixture(compounds, composition, interactions)
    low, high = bracket_pa
    low_result = _pressure_point(kind, compounds, composition, interactions, temperature_k, low, max_iterations)
    high_result = _pressure_point(kind, compounds, composition, interactions, temperature_k, high, max_iterations)
    if not low_result.report.converged and not high_result.report.converged:
        failed = low_result if not low_result.report.converged else high_result
        report = _report(False, failed.report.iterations, failed.report.residual, failed.report.failure_reason, f"PR {kind}-pressure bisection")
        return PressureResult(report, kind, temperature_k, failed.pressure_pa, failed.liquid_composition, failed.vapor_composition, failed.liquid_state, failed.vapor_state)
    if not low_result.report.converged or not high_result.report.converged:
        raise ValidationError("pressure bracket endpoints must admit distinct liquid and vapor roots")
    if low_result.report.residual == 0.0:
        return low_result
    if high_result.report.residual == 0.0:
        return high_result
    if low_result.report.residual * high_result.report.residual > 0.0:
        raise ValidationError("pressure bracket does not enclose a bubble/dew residual root")

    latest = low_result
    for iteration in range(1, max_iterations + 1):
        pressure_pa = (low + high) / 2.0
        latest = _pressure_point(kind, compounds, composition, interactions, temperature_k, pressure_pa, max_iterations)
        if not latest.report.converged:
            report = _report(False, iteration, latest.report.residual, latest.report.failure_reason, f"PR {kind}-pressure bisection")
            return PressureResult(report, kind, temperature_k, pressure_pa, latest.liquid_composition, latest.vapor_composition, latest.liquid_state, latest.vapor_state)
        if abs(latest.report.residual) <= 1e-8 and high - low <= 1e-8 * pressure_pa:
            report = _report(True, iteration, latest.report.residual, algorithm=f"PR {kind}-pressure bisection")
            return PressureResult(report, kind, temperature_k, pressure_pa, latest.liquid_composition, latest.vapor_composition, latest.liquid_state, latest.vapor_state)
        if low_result.report.residual * latest.report.residual <= 0.0:
            high, high_result = pressure_pa, latest
        else:
            low, low_result = pressure_pa, latest

    report = _report(False, max_iterations, latest.report.residual, f"{kind} pressure iteration limit reached", f"PR {kind}-pressure bisection")
    return PressureResult(report, kind, temperature_k, latest.pressure_pa, latest.liquid_composition, latest.vapor_composition, latest.liquid_state, latest.vapor_state)


def bubble_pressure(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], interactions: PRInteractions,
    temperature_k: float, bracket_pa: tuple[float, float], *, max_iterations: int = 100,
) -> PressureResult:
    return _pressure_solve("bubble", compounds, composition, interactions, temperature_k, bracket_pa, max_iterations)


def dew_pressure(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], interactions: PRInteractions,
    temperature_k: float, bracket_pa: tuple[float, float], *, max_iterations: int = 100,
) -> PressureResult:
    return _pressure_solve("dew", compounds, composition, interactions, temperature_k, bracket_pa, max_iterations)


def _validate_composition(compounds: tuple[Compound, ...], composition: tuple[float, ...]) -> None:
    if len(compounds) != len(composition) or not composition or any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0
        for value in composition
    ) or not math.isclose(math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("compound composition must be finite, non-negative, and sum to one")


def _validate_absolute_pressure(pressure_pa: float) -> None:
    if isinstance(pressure_pa, bool) or not isinstance(pressure_pa, (int, float)) or not math.isfinite(pressure_pa) or pressure_pa <= 0:
        raise ValidationError("absolute pressure must be finite and positive")


def ideal_mixture_enthalpy(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], correlations: tuple[IdealCorrelations, ...], temperature_k: float
) -> float:
    _validate_composition(compounds, composition)
    records = {record.compound_id: record for record in correlations}
    if any(compound.id not in records for compound in compounds):
        raise ValidationError("missing ideal heat-capacity correlation for a flash compound")
    return math.fsum(
        fraction * records[compound.id].enthalpy_change(temperature_k, 298.15).value
        for compound, fraction in zip(compounds, composition)
    )


def phase_enthalpy(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], correlations: tuple[IdealCorrelations, ...],
    temperature_k: float, state: PRMixtureState | None,
) -> float:
    if state is None or not math.isfinite(state.departure_enthalpy_j_per_kmol):
        raise ValidationError("a finite PR phase state is required for phase enthalpy")
    return ideal_mixture_enthalpy(compounds, composition, correlations, temperature_k) + state.departure_enthalpy_j_per_kmol


def phase_entropy(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], correlations: tuple[IdealCorrelations, ...],
    temperature_k: float, pressure_pa: float, state: PRMixtureState | None,
) -> float:
    _validate_composition(compounds, composition)
    if state is None or not math.isfinite(state.departure_entropy_j_per_kmol_k):
        raise ValidationError("a finite PR phase state is required for phase entropy")
    records = {record.compound_id: record for record in correlations}
    if any(compound.id not in records for compound in compounds):
        raise ValidationError("missing ideal heat-capacity correlation for a flash compound")
    ideal = math.fsum(
        fraction * records[compound.id].entropy_change(temperature_k, pressure_pa).value
        for compound, fraction in zip(compounds, composition)
    )
    mixing = -8314.46261815324 * math.fsum(fraction * math.log(fraction) for fraction in composition if fraction > 0.0)
    return ideal + mixing + state.departure_entropy_j_per_kmol_k


def flash_enthalpy(compounds: tuple[Compound, ...], correlations: tuple[IdealCorrelations, ...], flash: TPFlashResult) -> float:
    if not flash.report.converged:
        raise ValidationError("a converged TP flash is required for total enthalpy")
    if flash.phase == "two-phase":
        liquid = phase_enthalpy(compounds, flash.liquid_composition, correlations, flash.temperature_k, flash.liquid_state)
        vapor = phase_enthalpy(compounds, flash.vapor_composition, correlations, flash.temperature_k, flash.vapor_state)
        return (1.0 - flash.vapor_fraction) * liquid + flash.vapor_fraction * vapor
    if flash.phase == "liquid":
        return phase_enthalpy(compounds, flash.liquid_composition, correlations, flash.temperature_k, flash.liquid_state)
    return phase_enthalpy(compounds, flash.vapor_composition, correlations, flash.temperature_k, flash.vapor_state)


def flash_entropy(compounds: tuple[Compound, ...], correlations: tuple[IdealCorrelations, ...], flash: TPFlashResult) -> float:
    if not flash.report.converged:
        raise ValidationError("a converged TP flash is required for total entropy")
    if flash.phase == "two-phase":
        liquid = phase_entropy(compounds, flash.liquid_composition, correlations, flash.temperature_k, flash.pressure_pa, flash.liquid_state)
        vapor = phase_entropy(compounds, flash.vapor_composition, correlations, flash.temperature_k, flash.pressure_pa, flash.vapor_state)
        return (1.0 - flash.vapor_fraction) * liquid + flash.vapor_fraction * vapor
    if flash.phase == "liquid":
        return phase_entropy(compounds, flash.liquid_composition, correlations, flash.temperature_k, flash.pressure_pa, flash.liquid_state)
    return phase_entropy(compounds, flash.vapor_composition, correlations, flash.temperature_k, flash.pressure_pa, flash.vapor_state)


def _validate_temperature_bracket(bracket_k: tuple[float, float]) -> None:
    if len(bracket_k) != 2 or any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0
        for value in bracket_k
    ) or bracket_k[0] >= bracket_k[1]:
        raise ValidationError("temperature bracket must contain two increasing finite positive values")


def _ph_result(
    converged: bool, iterations: int, residual: float, reason: str | None, temperature_k: float,
    pressure_pa: float, target: float, bracket: tuple[float, float], inner_iterations: int,
    flash: TPFlashResult | None, enthalpy: float | None,
) -> PHFlashResult:
    report = _report(converged, iterations, residual, reason, "PR PH bisection")
    return PHFlashResult(report, temperature_k, pressure_pa, target, enthalpy, bracket, iterations, inner_iterations, flash)


def ph_flash(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...], pressure_pa: float, target_enthalpy_j_per_kmol: float,
    temperature_bracket_k: tuple[float, float], *, max_iterations: int = 100,
) -> PHFlashResult:
    _validate_solver_limits(max_iterations, 1e-10)
    _validate_temperature_bracket(temperature_bracket_k)
    _validate_absolute_pressure(pressure_pa)
    if isinstance(target_enthalpy_j_per_kmol, bool) or not isinstance(target_enthalpy_j_per_kmol, (int, float)) or not math.isfinite(target_enthalpy_j_per_kmol):
        raise ValidationError("target enthalpy must be finite")
    low, high = temperature_bracket_k
    inner_iterations = 0

    def trial(temperature_k: float) -> TPFlashResult:
        nonlocal inner_iterations
        flash = tp_flash(compounds, composition, interactions, temperature_k, pressure_pa, max_iterations=max_iterations)
        inner_iterations += flash.report.iterations
        return flash

    try:
        low_flash = trial(low)
        high_flash = trial(high)
    except ValidationError as error:
        return _ph_result(False, 0, math.inf, str(error), low, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, None, None)
    if not low_flash.report.converged or not high_flash.report.converged:
        failed = low_flash if not low_flash.report.converged else high_flash
        return _ph_result(False, 0, failed.report.residual, failed.report.failure_reason, low, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, failed, None)
    try:
        low_enthalpy = flash_enthalpy(compounds, correlations, low_flash)
        high_enthalpy = flash_enthalpy(compounds, correlations, high_flash)
    except ValidationError as error:
        return _ph_result(False, 0, math.inf, str(error), low, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, low_flash, None)
    low_residual = low_enthalpy - target_enthalpy_j_per_kmol
    high_residual = high_enthalpy - target_enthalpy_j_per_kmol
    if low_residual == 0.0:
        return _ph_result(True, 0, 0.0, None, low, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, low_flash, low_enthalpy)
    if high_residual == 0.0:
        return _ph_result(True, 0, 0.0, high, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, high_flash, high_enthalpy)
    if low_residual * high_residual > 0.0:
        return _ph_result(False, 0, min(abs(low_residual), abs(high_residual)), "temperature bracket does not enclose target enthalpy", low, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, low_flash, low_enthalpy)

    latest_flash, latest_enthalpy = low_flash, low_enthalpy
    tolerance = max(1e-6 * abs(target_enthalpy_j_per_kmol), 1e-3)
    for iteration in range(1, max_iterations + 1):
        temperature_k = (low + high) / 2.0
        try:
            latest_flash = trial(temperature_k)
        except ValidationError as error:
            return _ph_result(False, iteration, math.inf, str(error), temperature_k, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, None, None)
        if not latest_flash.report.converged:
            return _ph_result(False, iteration, latest_flash.report.residual, latest_flash.report.failure_reason, temperature_k, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, latest_flash, None)
        try:
            latest_enthalpy = flash_enthalpy(compounds, correlations, latest_flash)
        except ValidationError as error:
            return _ph_result(False, iteration, math.inf, str(error), temperature_k, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, latest_flash, None)
        residual = latest_enthalpy - target_enthalpy_j_per_kmol
        if abs(residual) <= tolerance and high - low <= 1e-8 * temperature_k:
            return _ph_result(True, iteration, residual, None, temperature_k, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, latest_flash, latest_enthalpy)
        if low_residual * residual <= 0.0:
            high, high_residual = temperature_k, residual
        else:
            low, low_residual = temperature_k, residual

    return _ph_result(False, max_iterations, latest_enthalpy - target_enthalpy_j_per_kmol, "PH iteration limit reached", temperature_k, pressure_pa, target_enthalpy_j_per_kmol, (low, high), inner_iterations, latest_flash, latest_enthalpy)


def ps_flash(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...], pressure_pa: float, target_entropy_j_per_kmol_k: float,
    temperature_bracket_k: tuple[float, float], *, max_iterations: int = 100,
) -> PSFlashResult:
    _validate_solver_limits(max_iterations, 1e-10)
    _validate_temperature_bracket(temperature_bracket_k)
    _validate_absolute_pressure(pressure_pa)
    if isinstance(target_entropy_j_per_kmol_k, bool) or not isinstance(target_entropy_j_per_kmol_k, (int, float)) or not math.isfinite(target_entropy_j_per_kmol_k):
        raise ValidationError("target entropy must be finite")
    low, high = temperature_bracket_k
    inner_iterations = 0

    def trial(temperature_k: float) -> TPFlashResult:
        nonlocal inner_iterations
        flash = tp_flash(compounds, composition, interactions, temperature_k, pressure_pa, max_iterations=max_iterations)
        inner_iterations += flash.report.iterations
        return flash

    try:
        low_flash, high_flash = trial(low), trial(high)
        low_entropy = flash_entropy(compounds, correlations, low_flash)
        high_entropy = flash_entropy(compounds, correlations, high_flash)
    except ValidationError as error:
        return PSFlashResult(_report(False, 0, math.inf, str(error), "PR PS bisection"), low, pressure_pa, target_entropy_j_per_kmol_k, None, (low, high), 0, inner_iterations, None)
    if not low_flash.report.converged or not high_flash.report.converged:
        failed = low_flash if not low_flash.report.converged else high_flash
        return PSFlashResult(_report(False, 0, failed.report.residual, failed.report.failure_reason, "PR PS bisection"), low, pressure_pa, target_entropy_j_per_kmol_k, None, (low, high), 0, inner_iterations, failed)
    low_residual, high_residual = low_entropy - target_entropy_j_per_kmol_k, high_entropy - target_entropy_j_per_kmol_k
    if low_residual * high_residual > 0.0:
        return PSFlashResult(_report(False, 0, min(abs(low_residual), abs(high_residual)), "temperature bracket does not enclose target entropy", "PR PS bisection"), low, pressure_pa, target_entropy_j_per_kmol_k, low_entropy, (low, high), 0, inner_iterations, low_flash)

    latest_flash, latest_entropy = low_flash, low_entropy
    tolerance = max(1e-6 * abs(target_entropy_j_per_kmol_k), 1e-3)
    for iteration in range(1, max_iterations + 1):
        temperature_k = (low + high) / 2.0
        try:
            latest_flash = trial(temperature_k)
            latest_entropy = flash_entropy(compounds, correlations, latest_flash)
        except ValidationError as error:
            return PSFlashResult(_report(False, iteration, math.inf, str(error), "PR PS bisection"), temperature_k, pressure_pa, target_entropy_j_per_kmol_k, None, (low, high), iteration, inner_iterations, None)
        if not latest_flash.report.converged:
            return PSFlashResult(_report(False, iteration, latest_flash.report.residual, latest_flash.report.failure_reason, "PR PS bisection"), temperature_k, pressure_pa, target_entropy_j_per_kmol_k, None, (low, high), iteration, inner_iterations, latest_flash)
        residual = latest_entropy - target_entropy_j_per_kmol_k
        if abs(residual) <= tolerance and high - low <= 1e-8 * temperature_k:
            return PSFlashResult(_report(True, iteration, residual, algorithm="PR PS bisection"), temperature_k, pressure_pa, target_entropy_j_per_kmol_k, latest_entropy, (low, high), iteration, inner_iterations, latest_flash)
        if low_residual * residual <= 0.0:
            high, high_residual = temperature_k, residual
        else:
            low, low_residual = temperature_k, residual
    return PSFlashResult(_report(False, max_iterations, latest_entropy - target_entropy_j_per_kmol_k, "PS iteration limit reached", "PR PS bisection"), temperature_k, pressure_pa, target_entropy_j_per_kmol_k, latest_entropy, (low, high), max_iterations, inner_iterations, latest_flash)
