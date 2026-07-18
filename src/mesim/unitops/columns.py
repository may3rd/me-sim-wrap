from dataclasses import dataclass
import math

from scipy.optimize import least_squares

from ..errors import ValidationError
from ..thermo.activity import NRTLVLEResult
from ..thermo.systems import NRTLSystem


def _finite_number(value: object) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
    )


@dataclass(frozen=True, slots=True)
class NRTLColumnEquilibriumProfile:
    equilibrium_ratios: tuple[tuple[float, ...], ...]
    bubble_pressures_pa: tuple[float, ...]
    relative_pressure_residuals: tuple[float, ...]


def nrtl_column_equilibrium_profile(
    data: NRTLSystem,
    compound_ids: tuple[str, ...],
    temperatures_k: tuple[float, ...],
    pressures_pa: tuple[float, ...],
    liquid_mole_fractions: tuple[tuple[float, ...], ...],
) -> NRTLColumnEquilibriumProfile:
    try:
        temperatures = tuple(temperatures_k)
        pressures = tuple(pressures_pa)
        liquid = tuple(tuple(row) for row in liquid_mole_fractions)
    except TypeError as error:
        raise ValidationError("NRTL column profile inputs must be finite sequences") from error
    if (
        not isinstance(data, NRTLSystem)
        or tuple(compound_ids) != data.compound_ids
        or not temperatures
        or len(temperatures) != len(pressures)
        or len(temperatures) != len(liquid)
    ):
        raise ValidationError("NRTL column profile stage arrays must be non-empty and aligned")
    ratios = []
    bubble_pressures = []
    residuals = []
    for temperature_k, pressure_pa, composition in zip(temperatures, pressures, liquid):
        stage_ratios = data.equilibrium_ratios(composition, temperature_k, pressure_pa)
        bubble = data.bubble_pressure(composition, temperature_k)
        ratios.append(stage_ratios)
        bubble_pressures.append(bubble.pressure_pa)
        residuals.append((bubble.pressure_pa - pressure_pa) / pressure_pa)
    return NRTLColumnEquilibriumProfile(tuple(ratios), tuple(bubble_pressures), tuple(residuals))


def nrtl_column_bubble_temperature_profile(
    data: NRTLSystem,
    compound_ids: tuple[str, ...],
    pressures_pa: tuple[float, ...],
    liquid_mole_fractions: tuple[tuple[float, ...], ...],
    bracket_k: tuple[float, float],
    max_iterations: int = 100,
    tolerance: float = 1e-10,
) -> tuple[NRTLVLEResult, ...]:
    try:
        pressures = tuple(pressures_pa)
        liquid = tuple(tuple(row) for row in liquid_mole_fractions)
    except TypeError as error:
        raise ValidationError("NRTL column bubble-point inputs must be finite sequences") from error
    if (
        not isinstance(data, NRTLSystem)
        or tuple(compound_ids) != data.compound_ids
        or not pressures
        or len(pressures) != len(liquid)
    ):
        raise ValidationError("NRTL column pressure and composition arrays must be non-empty and aligned")
    return tuple(
        data.bubble_temperature(
            composition,
            pressure_pa,
            bracket_k,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )
        for pressure_pa, composition in zip(pressures, liquid)
    )


@dataclass(frozen=True, slots=True)
class NRTLColumnEnthalpyProfile:
    liquid_enthalpies_j_per_kmol: tuple[float, ...]
    vapor_enthalpies_j_per_kmol: tuple[float, ...]
    liquid_densities_kg_per_m3: tuple[float, ...]
    excess_enthalpies_j_per_kmol: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class NRTLRigorousColumnIteration:
    residual_evaluation: int
    scaled_residual_norm: float


@dataclass(frozen=True, slots=True)
class NRTLRigorousColumnResult:
    temperatures_k: tuple[float, ...]
    pressures_pa: tuple[float, ...]
    liquid_flows_kmol_s: tuple[float, ...]
    vapor_flows_kmol_s: tuple[float, ...]
    liquid_mole_fractions: tuple[tuple[float, ...], ...]
    vapor_mole_fractions: tuple[tuple[float, ...], ...]
    equilibrium_ratios: tuple[tuple[float, ...], ...]
    liquid_enthalpies_j_per_kmol: tuple[float, ...]
    vapor_enthalpies_j_per_kmol: tuple[float, ...]
    distillate_flow_kmol_s: float
    bottoms_flow_kmol_s: float
    condenser_duty_w: float
    reboiler_duty_w: float
    scaled_residual_norm: float
    solver_evaluations: int
    residual_evaluations: int
    history: tuple[NRTLRigorousColumnIteration, ...]


@dataclass(frozen=True, slots=True)
class NRTLRigorousReboiledAbsorberResult:
    temperatures_k: tuple[float, ...]
    pressures_pa: tuple[float, ...]
    liquid_flows_kmol_s: tuple[float, ...]
    vapor_flows_kmol_s: tuple[float, ...]
    liquid_mole_fractions: tuple[tuple[float, ...], ...]
    vapor_mole_fractions: tuple[tuple[float, ...], ...]
    equilibrium_ratios: tuple[tuple[float, ...], ...]
    liquid_enthalpies_j_per_kmol: tuple[float, ...]
    vapor_enthalpies_j_per_kmol: tuple[float, ...]
    overhead_vapor_flow_kmol_s: float
    bottoms_flow_kmol_s: float
    reboiler_duty_w: float
    scaled_residual_norm: float
    solver_evaluations: int
    residual_evaluations: int
    history: tuple[NRTLRigorousColumnIteration, ...]


class NRTLRigorousColumnConvergenceError(RuntimeError):
    def __init__(self, message: str, history: tuple[NRTLRigorousColumnIteration, ...]):
        super().__init__(message)
        self.history = history


def nrtl_column_enthalpy_profile(
    data: NRTLSystem,
    compound_ids: tuple[str, ...],
    temperatures_k: tuple[float, ...],
    pressures_pa: tuple[float, ...],
    liquid_mole_fractions: tuple[tuple[float, ...], ...],
    vapor_mole_fractions: tuple[tuple[float, ...], ...],
) -> NRTLColumnEnthalpyProfile:
    try:
        temperatures = tuple(temperatures_k)
        pressures = tuple(pressures_pa)
        liquid = tuple(tuple(row) for row in liquid_mole_fractions)
        vapor = tuple(tuple(row) for row in vapor_mole_fractions)
    except TypeError as error:
        raise ValidationError("NRTL column enthalpy inputs must be finite sequences") from error
    stage_count = len(temperatures)
    if (
        not isinstance(data, NRTLSystem)
        or tuple(compound_ids) != data.compound_ids
        or not temperatures
        or len(pressures) != stage_count
        or len(liquid) != stage_count
        or len(vapor) != stage_count
    ):
        raise ValidationError("NRTL column enthalpy stage arrays must be non-empty and aligned")
    results = tuple(
        data.phase_enthalpies(
            liquid_row,
            vapor_row,
            temperature_k,
            pressure_pa,
        )
        for temperature_k, pressure_pa, liquid_row, vapor_row in zip(
            temperatures, pressures, liquid, vapor
        )
    )
    return NRTLColumnEnthalpyProfile(
        tuple(result.liquid_j_per_kmol for result in results),
        tuple(result.vapor_j_per_kmol for result in results),
        tuple(result.liquid_density_kg_per_m3 for result in results),
        tuple(result.excess_enthalpy_j_per_kmol for result in results),
    )


def nrtl_rigorous_total_condenser_column(
    data: NRTLSystem,
    compound_ids: tuple[str, ...],
    feed_component_flows_by_stage_kmol_s: tuple[tuple[float, ...], ...],
    feed_energy_flows_by_stage_w: tuple[float, ...],
    pressures_pa: tuple[float, ...],
    initial_temperatures_k: tuple[float, ...],
    initial_liquid_flows_kmol_s: tuple[float, ...],
    initial_vapor_flows_kmol_s: tuple[float, ...],
    initial_liquid_mole_fractions: tuple[tuple[float, ...], ...],
    reflux_ratio: float,
    bottoms_flow_kmol_s: float,
    temperature_bounds_k: tuple[float, float],
    residual_tolerance: float = 1.0e-8,
    maximum_solver_evaluations: int = 100,
) -> NRTLRigorousColumnResult:
    """Solve one total-condenser NRTL column with simultaneous MESH equations.

    The current predictive slice fixes a constant supplied pressure profile,
    one liquid distillate at stage zero, one liquid bottoms product at the last
    stage, a reflux-ratio specification, and a bottoms-flow specification.
    Condenser and reboiler duties are calculated outputs. Feed energy arrives
    on the same watt basis as stage energy balances.
    """
    try:
        ids = tuple(compound_ids)
        feeds = tuple(tuple(float(value) for value in row) for row in feed_component_flows_by_stage_kmol_s)
        feed_energy = tuple(float(value) for value in feed_energy_flows_by_stage_w)
        pressures = tuple(float(value) for value in pressures_pa)
        initial_temperatures = tuple(float(value) for value in initial_temperatures_k)
        initial_liquid = tuple(float(value) for value in initial_liquid_flows_kmol_s)
        initial_vapor = tuple(float(value) for value in initial_vapor_flows_kmol_s)
        initial_fractions = tuple(
            tuple(float(value) for value in row) for row in initial_liquid_mole_fractions
        )
        temperature_bounds = tuple(float(value) for value in temperature_bounds_k)
    except (TypeError, ValueError) as error:
        raise ValidationError("NRTL rigorous-column inputs must be finite sequences") from error
    stage_count = len(feeds)
    component_count = len(ids)
    scalar_values = (reflux_ratio, bottoms_flow_kmol_s, residual_tolerance)
    if (
        not isinstance(data, NRTLSystem)
        or ids != data.compound_ids
        or stage_count < 2
        or component_count < 2
        or len(set(ids)) != component_count
        or any(not isinstance(value, str) or not value for value in ids)
        or len(feed_energy) != stage_count
        or len(pressures) != stage_count
        or len(initial_temperatures) != stage_count
        or len(initial_liquid) != stage_count
        or len(initial_vapor) != stage_count
        or len(initial_fractions) != stage_count
        or len(temperature_bounds) != 2
        or any(len(row) != component_count for row in feeds + initial_fractions)
        or any(not _finite_number(value) or value < 0.0 for row in feeds for value in row)
        or any(not _finite_number(value) for value in feed_energy)
        or any(not _finite_number(value) or value <= 0.0 for value in pressures)
        or any(not _finite_number(value) or value <= 0.0 for value in initial_temperatures)
        or any(not _finite_number(value) or value <= 0.0 for value in initial_liquid)
        or initial_vapor[0] != 0.0
        or any(not _finite_number(value) or value <= 0.0 for value in initial_vapor[1:])
        or any(
            not _finite_number(value) or value <= 0.0
            for row in initial_fractions for value in row
        )
        or any(abs(math.fsum(row) - 1.0) > 1.0e-8 for row in initial_fractions)
        or any(not _finite_number(value) or value <= 0.0 for value in scalar_values)
        or not all(_finite_number(value) and value > 0.0 for value in temperature_bounds)
        or temperature_bounds[0] >= temperature_bounds[1]
        or any(not temperature_bounds[0] <= value <= temperature_bounds[1] for value in initial_temperatures)
        or isinstance(maximum_solver_evaluations, bool)
        or not isinstance(maximum_solver_evaluations, int)
        or maximum_solver_evaluations <= 0
    ):
        raise ValidationError("NRTL rigorous total-condenser column inputs are invalid")
    flow_scale = math.fsum(value for row in feeds for value in row)
    if flow_scale <= 0.0 or bottoms_flow_kmol_s >= flow_scale:
        raise ValidationError("NRTL rigorous column requires feed flow above the bottoms specification")
    distillate_initial = flow_scale - bottoms_flow_kmol_s

    composition_offset = 2 * stage_count - 1
    temperature_offset = composition_offset + stage_count * (component_count - 1)
    distillate_index = temperature_offset + stage_count
    condenser_index = distillate_index + 1
    reboiler_index = condenser_index + 1

    initial_logits = [
        math.log(row[component]) - math.log(row[-1])
        for row in initial_fractions
        for component in range(component_count - 1)
    ]
    variables = (
        [math.log(value) for value in initial_liquid]
        + [math.log(value) for value in initial_vapor[1:]]
        + initial_logits
        + list(initial_temperatures)
        + [math.log(distillate_initial)]
    )

    def decode(values):
        liquid = tuple(math.exp(float(value)) for value in values[:stage_count])
        vapor = (0.0,) + tuple(
            math.exp(float(value))
            for value in values[stage_count:composition_offset]
        )
        fractions = []
        reduced_count = component_count - 1
        for stage in range(stage_count):
            start = composition_offset + stage * reduced_count
            logits = [float(value) for value in values[start:start + reduced_count]] + [0.0]
            maximum = max(logits)
            exponentials = [math.exp(value - maximum) for value in logits]
            total = math.fsum(exponentials)
            fractions.append(tuple(value / total for value in exponentials))
        temperatures = tuple(
            float(value) for value in values[temperature_offset:distillate_index]
        )
        distillate = math.exp(float(values[distillate_index]))
        condenser_heat_input = float(values[condenser_index])
        reboiler_heat_input = float(values[reboiler_index])
        return (
            liquid,
            vapor,
            tuple(fractions),
            temperatures,
            distillate,
            condenser_heat_input,
            reboiler_heat_input,
        )

    def thermodynamics(temperatures, fractions):
        ratios = []
        vapor_fractions = []
        liquid_enthalpies = []
        vapor_enthalpies = []
        for temperature_k, pressure_pa, liquid_row in zip(
            temperatures, pressures, fractions
        ):
            ratio_row = data.equilibrium_ratios(
                liquid_row, temperature_k, pressure_pa
            )
            vapor_row = tuple(
                ratio * fraction for ratio, fraction in zip(ratio_row, liquid_row)
            )
            vapor_total = math.fsum(vapor_row)
            normalized_vapor = tuple(value / vapor_total for value in vapor_row)
            enthalpies = data.phase_enthalpies(
                liquid_row,
                normalized_vapor,
                temperature_k,
                pressure_pa,
            )
            ratios.append(ratio_row)
            vapor_fractions.append(vapor_row)
            liquid_enthalpies.append(enthalpies.liquid_j_per_kmol)
            vapor_enthalpies.append(enthalpies.vapor_j_per_kmol)
        return (
            tuple(ratios),
            tuple(vapor_fractions),
            tuple(liquid_enthalpies),
            tuple(vapor_enthalpies),
        )

    _, _, initial_hl, initial_hv = thermodynamics(
        initial_temperatures, initial_fractions
    )
    condenser_initial_w = (
        (initial_liquid[0] + distillate_initial) * initial_hl[0]
        - initial_vapor[1] * initial_hv[1]
        - feed_energy[0]
    )
    reboiler_initial_w = (
        initial_liquid[-1] * initial_hl[-1]
        + initial_vapor[-1] * initial_hv[-1]
        - initial_liquid[-2] * initial_hl[-2]
        - feed_energy[-1]
    )
    if condenser_initial_w >= 0.0 or reboiler_initial_w <= 0.0:
        raise ValidationError("initial column profile does not imply condenser removal and reboiler input")
    energy_scale = max(abs(condenser_initial_w), reboiler_initial_w, 1.0)
    variables.extend((condenser_initial_w / energy_scale, reboiler_initial_w / energy_scale))

    lower_flow = max(flow_scale * 1.0e-12, 1.0e-300)
    upper_flow = max(
        flow_scale * 1.0e4,
        max(initial_liquid + initial_vapor[1:] + (distillate_initial,)) * 10.0,
    )
    lower_logit = min(-50.0, min(initial_logits) - 1.0)
    upper_logit = max(50.0, max(initial_logits) + 1.0)
    lower_bounds = (
        [math.log(lower_flow)] * stage_count
        + [math.log(lower_flow)] * (stage_count - 1)
        + [lower_logit] * (stage_count * (component_count - 1))
        + [temperature_bounds[0]] * stage_count
        + [math.log(lower_flow), -100.0, 0.0]
    )
    upper_bounds = (
        [math.log(upper_flow)] * stage_count
        + [math.log(upper_flow)] * (stage_count - 1)
        + [upper_logit] * (stage_count * (component_count - 1))
        + [temperature_bounds[1]] * stage_count
        + [math.log(upper_flow), 0.0, 100.0]
    )

    history: list[NRTLRigorousColumnIteration] = []
    residual_evaluation = 0
    best_norm = math.inf

    def residuals(values):
        nonlocal residual_evaluation, best_norm
        residual_evaluation += 1
        liquid, vapor, fractions, temperatures, distillate, condenser, reboiler = decode(values)
        ratios, vapor_fractions, liquid_enthalpies, vapor_enthalpies = thermodynamics(
            temperatures, fractions
        )
        result = []
        for stage in range(stage_count):
            for component in range(component_count):
                residual = feeds[stage][component]
                if stage > 0:
                    residual += liquid[stage - 1] * fractions[stage - 1][component]
                if stage < stage_count - 1:
                    residual += vapor[stage + 1] * vapor_fractions[stage + 1][component]
                residual -= liquid[stage] * fractions[stage][component]
                residual -= vapor[stage] * vapor_fractions[stage][component]
                if stage == 0:
                    residual -= distillate * fractions[stage][component]
                result.append(residual / flow_scale)
        result.extend(math.fsum(row) - 1.0 for row in vapor_fractions)
        for stage in range(stage_count):
            residual = feed_energy[stage]
            if stage > 0:
                residual += liquid[stage - 1] * liquid_enthalpies[stage - 1]
            if stage < stage_count - 1:
                residual += vapor[stage + 1] * vapor_enthalpies[stage + 1]
            residual -= liquid[stage] * liquid_enthalpies[stage]
            residual -= vapor[stage] * vapor_enthalpies[stage]
            if stage == 0:
                residual -= distillate * liquid_enthalpies[stage]
                residual += condenser * energy_scale
            if stage == stage_count - 1:
                residual += reboiler * energy_scale
            result.append(residual / energy_scale)
        result.append((liquid[0] - reflux_ratio * distillate) / flow_scale)
        result.append((liquid[-1] - bottoms_flow_kmol_s) / flow_scale)
        norm = max(abs(value) for value in result)
        if math.isfinite(norm) and norm < best_norm:
            best_norm = norm
            history.append(NRTLRigorousColumnIteration(residual_evaluation, norm))
        return tuple(result)

    try:
        solved = least_squares(
            residuals,
            variables,
            bounds=(lower_bounds, upper_bounds),
            x_scale="jac",
            ftol=1.0e-11,
            xtol=1.0e-11,
            gtol=1.0e-11,
            max_nfev=maximum_solver_evaluations,
        )
    except (ValidationError, ValueError, OverflowError) as error:
        raise NRTLRigorousColumnConvergenceError(str(error), tuple(history)) from error
    final_residuals = residuals(solved.x)
    final_norm = max(abs(value) for value in final_residuals)
    if not math.isfinite(final_norm) or final_norm > residual_tolerance:
        raise NRTLRigorousColumnConvergenceError(
            "NRTL rigorous total-condenser column did not converge", tuple(history)
        )
    liquid, vapor, fractions, temperatures, distillate, condenser, reboiler = decode(solved.x)
    ratios, vapor_fractions, liquid_enthalpies, vapor_enthalpies = thermodynamics(
        temperatures, fractions
    )
    return NRTLRigorousColumnResult(
        temperatures,
        pressures,
        liquid,
        vapor,
        fractions,
        vapor_fractions,
        ratios,
        liquid_enthalpies,
        vapor_enthalpies,
        distillate,
        liquid[-1],
        -condenser * energy_scale,
        reboiler * energy_scale,
        final_norm,
        int(solved.nfev),
        residual_evaluation,
        tuple(history),
    )


def nrtl_rigorous_reboiled_absorber(
    data: NRTLSystem,
    compound_ids: tuple[str, ...],
    feed_component_flows_by_stage_kmol_s: tuple[tuple[float, ...], ...],
    feed_energy_flows_by_stage_w: tuple[float, ...],
    pressures_pa: tuple[float, ...],
    initial_temperatures_k: tuple[float, ...],
    initial_liquid_flows_kmol_s: tuple[float, ...],
    initial_vapor_flows_kmol_s: tuple[float, ...],
    initial_liquid_mole_fractions: tuple[tuple[float, ...], ...],
    bottoms_flow_kmol_s: float,
    temperature_bounds_k: tuple[float, float],
    residual_tolerance: float = 1.0e-8,
    maximum_solver_evaluations: int = 100,
) -> NRTLRigorousReboiledAbsorberResult:
    """Solve one NRTL reboiled absorber with simultaneous MESH equations.

    Vapor leaves stage zero as the overhead product, liquid leaves the last
    stage as bottoms, and the positive reboiler input is calculated. The
    current predictive slice uses a supplied pressure profile, one feed, and
    a specified bottoms flow; it has no condenser, reflux, or liquid overhead.
    """
    try:
        ids = tuple(compound_ids)
        feeds = tuple(
            tuple(float(value) for value in row)
            for row in feed_component_flows_by_stage_kmol_s
        )
        feed_energy = tuple(float(value) for value in feed_energy_flows_by_stage_w)
        pressures = tuple(float(value) for value in pressures_pa)
        initial_temperatures = tuple(float(value) for value in initial_temperatures_k)
        initial_liquid = tuple(float(value) for value in initial_liquid_flows_kmol_s)
        initial_vapor = tuple(float(value) for value in initial_vapor_flows_kmol_s)
        initial_fractions = tuple(
            tuple(float(value) for value in row)
            for row in initial_liquid_mole_fractions
        )
        temperature_bounds = tuple(float(value) for value in temperature_bounds_k)
    except (TypeError, ValueError) as error:
        raise ValidationError("NRTL reboiled-absorber inputs must be finite sequences") from error

    stage_count = len(feeds)
    component_count = len(ids)
    scalar_values = (bottoms_flow_kmol_s, residual_tolerance)
    if (
        not isinstance(data, NRTLSystem)
        or ids != data.compound_ids
        or stage_count < 2
        or component_count < 2
        or len(set(ids)) != component_count
        or any(not isinstance(value, str) or not value for value in ids)
        or len(feed_energy) != stage_count
        or len(pressures) != stage_count
        or len(initial_temperatures) != stage_count
        or len(initial_liquid) != stage_count
        or len(initial_vapor) != stage_count
        or len(initial_fractions) != stage_count
        or len(temperature_bounds) != 2
        or any(len(row) != component_count for row in feeds + initial_fractions)
        or any(not _finite_number(value) or value < 0.0 for row in feeds for value in row)
        or any(not _finite_number(value) for value in feed_energy)
        or any(not _finite_number(value) or value <= 0.0 for value in pressures)
        or any(not _finite_number(value) or value <= 0.0 for value in initial_temperatures)
        or any(not _finite_number(value) or value <= 0.0 for value in initial_liquid)
        or any(not _finite_number(value) or value <= 0.0 for value in initial_vapor)
        or any(
            not _finite_number(value) or value <= 0.0
            for row in initial_fractions for value in row
        )
        or any(abs(math.fsum(row) - 1.0) > 1.0e-8 for row in initial_fractions)
        or any(not _finite_number(value) or value <= 0.0 for value in scalar_values)
        or not all(_finite_number(value) and value > 0.0 for value in temperature_bounds)
        or temperature_bounds[0] >= temperature_bounds[1]
        or any(
            not temperature_bounds[0] <= value <= temperature_bounds[1]
            for value in initial_temperatures
        )
        or isinstance(maximum_solver_evaluations, bool)
        or not isinstance(maximum_solver_evaluations, int)
        or maximum_solver_evaluations <= 0
    ):
        raise ValidationError("NRTL rigorous reboiled-absorber inputs are invalid")

    flow_scale = math.fsum(value for row in feeds for value in row)
    if flow_scale <= 0.0 or bottoms_flow_kmol_s >= flow_scale:
        raise ValidationError("NRTL reboiled absorber requires feed flow above the bottoms specification")

    composition_offset = 2 * stage_count
    temperature_offset = composition_offset + stage_count * (component_count - 1)
    reboiler_index = temperature_offset + stage_count
    initial_logits = [
        math.log(row[component]) - math.log(row[-1])
        for row in initial_fractions
        for component in range(component_count - 1)
    ]
    variables = (
        [math.log(value) for value in initial_liquid]
        + [math.log(value) for value in initial_vapor]
        + initial_logits
        + list(initial_temperatures)
    )

    def decode(values):
        liquid = tuple(math.exp(float(value)) for value in values[:stage_count])
        vapor = tuple(
            math.exp(float(value))
            for value in values[stage_count:composition_offset]
        )
        fractions = []
        reduced_count = component_count - 1
        for stage in range(stage_count):
            start = composition_offset + stage * reduced_count
            logits = [float(value) for value in values[start:start + reduced_count]] + [0.0]
            maximum = max(logits)
            exponentials = [math.exp(value - maximum) for value in logits]
            total = math.fsum(exponentials)
            fractions.append(tuple(value / total for value in exponentials))
        temperatures = tuple(float(value) for value in values[temperature_offset:reboiler_index])
        reboiler_heat_input = float(values[reboiler_index])
        return liquid, vapor, tuple(fractions), temperatures, reboiler_heat_input

    def thermodynamics(temperatures, fractions):
        ratios = []
        vapor_fractions = []
        liquid_enthalpies = []
        vapor_enthalpies = []
        for temperature_k, pressure_pa, liquid_row in zip(
            temperatures, pressures, fractions
        ):
            ratio_row = data.equilibrium_ratios(
                liquid_row, temperature_k, pressure_pa
            )
            vapor_row = tuple(
                ratio * fraction for ratio, fraction in zip(ratio_row, liquid_row)
            )
            vapor_total = math.fsum(vapor_row)
            normalized_vapor = tuple(value / vapor_total for value in vapor_row)
            enthalpies = data.phase_enthalpies(
                liquid_row,
                normalized_vapor,
                temperature_k,
                pressure_pa,
            )
            ratios.append(ratio_row)
            vapor_fractions.append(vapor_row)
            liquid_enthalpies.append(enthalpies.liquid_j_per_kmol)
            vapor_enthalpies.append(enthalpies.vapor_j_per_kmol)
        return (
            tuple(ratios),
            tuple(vapor_fractions),
            tuple(liquid_enthalpies),
            tuple(vapor_enthalpies),
        )

    _, _, initial_hl, initial_hv = thermodynamics(
        initial_temperatures, initial_fractions
    )
    reboiler_initial_w = (
        initial_liquid[-1] * initial_hl[-1]
        + initial_vapor[-1] * initial_hv[-1]
        - initial_liquid[-2] * initial_hl[-2]
        - feed_energy[-1]
    )
    if not _finite_number(reboiler_initial_w) or reboiler_initial_w <= 0.0:
        raise ValidationError("initial reboiled-absorber profile does not imply reboiler input")
    energy_scale = max(reboiler_initial_w, 1.0)
    variables.append(reboiler_initial_w / energy_scale)

    lower_flow = max(flow_scale * 1.0e-12, 1.0e-300)
    upper_flow = max(
        flow_scale * 1.0e4,
        max(initial_liquid + initial_vapor) * 10.0,
    )
    lower_logit = min(-50.0, min(initial_logits) - 1.0)
    upper_logit = max(50.0, max(initial_logits) + 1.0)
    lower_bounds = (
        [math.log(lower_flow)] * (2 * stage_count)
        + [lower_logit] * (stage_count * (component_count - 1))
        + [temperature_bounds[0]] * stage_count
        + [0.0]
    )
    upper_bounds = (
        [math.log(upper_flow)] * (2 * stage_count)
        + [upper_logit] * (stage_count * (component_count - 1))
        + [temperature_bounds[1]] * stage_count
        + [100.0]
    )

    history: list[NRTLRigorousColumnIteration] = []
    residual_evaluation = 0
    best_norm = math.inf

    def residuals(values):
        nonlocal residual_evaluation, best_norm
        residual_evaluation += 1
        liquid, vapor, fractions, temperatures, reboiler = decode(values)
        ratios, vapor_fractions, liquid_enthalpies, vapor_enthalpies = thermodynamics(
            temperatures, fractions
        )
        result = []
        for stage in range(stage_count):
            for component in range(component_count):
                residual = feeds[stage][component]
                if stage > 0:
                    residual += liquid[stage - 1] * fractions[stage - 1][component]
                if stage < stage_count - 1:
                    residual += vapor[stage + 1] * vapor_fractions[stage + 1][component]
                residual -= liquid[stage] * fractions[stage][component]
                residual -= vapor[stage] * vapor_fractions[stage][component]
                result.append(residual / flow_scale)
        result.extend(math.fsum(row) - 1.0 for row in vapor_fractions)
        for stage in range(stage_count):
            residual = feed_energy[stage]
            if stage > 0:
                residual += liquid[stage - 1] * liquid_enthalpies[stage - 1]
            if stage < stage_count - 1:
                residual += vapor[stage + 1] * vapor_enthalpies[stage + 1]
            residual -= liquid[stage] * liquid_enthalpies[stage]
            residual -= vapor[stage] * vapor_enthalpies[stage]
            if stage == stage_count - 1:
                residual += reboiler * energy_scale
            result.append(residual / energy_scale)
        result.append((liquid[-1] - bottoms_flow_kmol_s) / flow_scale)
        norm = max(abs(value) for value in result)
        if math.isfinite(norm) and norm < best_norm:
            best_norm = norm
            history.append(NRTLRigorousColumnIteration(residual_evaluation, norm))
        return tuple(result)

    try:
        solved = least_squares(
            residuals,
            variables,
            bounds=(lower_bounds, upper_bounds),
            x_scale="jac",
            ftol=1.0e-11,
            xtol=1.0e-11,
            gtol=1.0e-11,
            max_nfev=maximum_solver_evaluations,
        )
    except (ValidationError, ValueError, OverflowError) as error:
        raise NRTLRigorousColumnConvergenceError(str(error), tuple(history)) from error
    final_residuals = residuals(solved.x)
    final_norm = max(abs(value) for value in final_residuals)
    if not math.isfinite(final_norm) or final_norm > residual_tolerance:
        raise NRTLRigorousColumnConvergenceError(
            "NRTL rigorous reboiled absorber did not converge", tuple(history)
        )
    liquid, vapor, fractions, temperatures, reboiler = decode(solved.x)
    ratios, vapor_fractions, liquid_enthalpies, vapor_enthalpies = thermodynamics(
        temperatures, fractions
    )
    return NRTLRigorousReboiledAbsorberResult(
        temperatures,
        pressures,
        liquid,
        vapor,
        fractions,
        vapor_fractions,
        ratios,
        liquid_enthalpies,
        vapor_enthalpies,
        vapor[0],
        liquid[-1],
        reboiler * energy_scale,
        final_norm,
        int(solved.nfev),
        residual_evaluation,
        tuple(history),
    )


@dataclass(frozen=True, slots=True)
class StageStream:
    """One stage-boundary stream on a coherent kmol/s and J/kmol basis."""

    molar_flow_kmol_s: float
    mole_fractions: tuple[float, ...]
    molar_enthalpy_j_per_kmol: float

    def __post_init__(self) -> None:
        if (
            not _finite_number(self.molar_flow_kmol_s)
            or self.molar_flow_kmol_s < 0.0
            or not _finite_number(self.molar_enthalpy_j_per_kmol)
        ):
            raise ValidationError("stage stream flow and enthalpy are invalid")
        if not self.mole_fractions or any(
            not _finite_number(value) or value < 0.0 for value in self.mole_fractions
        ):
            raise ValidationError("stage stream trial mole fractions are invalid")


@dataclass(frozen=True, slots=True)
class EquilibriumStageState:
    temperature_k: float
    pressure_pa: float
    liquid: StageStream
    vapor: StageStream

    def __post_init__(self) -> None:
        if (
            not _finite_number(self.temperature_k)
            or not _finite_number(self.pressure_pa)
            or self.temperature_k <= 0.0
            or self.pressure_pa <= 0.0
        ):
            raise ValidationError("stage temperature and pressure must be positive finite values")
        if not isinstance(self.liquid, StageStream) or not isinstance(self.vapor, StageStream):
            raise ValidationError("stage liquid and vapor states are required")
        if len(self.liquid.mole_fractions) != len(self.vapor.mole_fractions):
            raise ValidationError("stage liquid and vapor component counts must match")


@dataclass(frozen=True, slots=True)
class StageResiduals:
    component_material_kmol_s: tuple[float, ...]
    phase_equilibrium: tuple[float, ...]
    liquid_summation: float
    vapor_summation: float
    energy_w: float

    def is_closed(
        self,
        component_tolerance_kmol_s: float,
        equilibrium_tolerance: float,
        summation_tolerance: float,
        energy_tolerance_w: float,
    ) -> bool:
        tolerances = (
            component_tolerance_kmol_s,
            equilibrium_tolerance,
            summation_tolerance,
            energy_tolerance_w,
        )
        if any(not _finite_number(value) or value <= 0.0 for value in tolerances):
            raise ValidationError("stage closure tolerances must be positive finite values")
        return (
            max(abs(value) for value in self.component_material_kmol_s)
            <= component_tolerance_kmol_s
            and max(abs(value) for value in self.phase_equilibrium)
            <= equilibrium_tolerance
            and abs(self.liquid_summation) <= summation_tolerance
            and abs(self.vapor_summation) <= summation_tolerance
            and abs(self.energy_w) <= energy_tolerance_w
        )


def equilibrium_stage_residuals(
    state: EquilibriumStageState,
    incoming_streams: tuple[StageStream, ...],
    equilibrium_ratios: tuple[float, ...],
    heat_duty_w: float = 0.0,
) -> StageResiduals:
    """Return component M, phase E, summation S, and energy H residuals."""
    if not isinstance(state, EquilibriumStageState):
        raise ValidationError("an equilibrium stage state is required")
    try:
        incoming = tuple(incoming_streams)
        ratios = tuple(equilibrium_ratios)
    except TypeError as error:
        raise ValidationError("stage incoming streams and K values must be finite sequences") from error
    component_count = len(state.liquid.mole_fractions)
    if (
        not incoming
        or any(not isinstance(stream, StageStream) for stream in incoming)
        or any(len(stream.mole_fractions) != component_count for stream in incoming)
        or len(ratios) != component_count
        or any(not _finite_number(value) or value <= 0.0 for value in ratios)
        or not _finite_number(heat_duty_w)
    ):
        raise ValidationError("stage inputs, component counts, K values, or duty are invalid")

    liquid = state.liquid
    vapor = state.vapor
    component_material = tuple(
        math.fsum(
            stream.molar_flow_kmol_s * stream.mole_fractions[index]
            for stream in incoming
        )
        - liquid.molar_flow_kmol_s * liquid.mole_fractions[index]
        - vapor.molar_flow_kmol_s * vapor.mole_fractions[index]
        for index in range(component_count)
    )
    phase_equilibrium = tuple(
        vapor.mole_fractions[index] - ratios[index] * liquid.mole_fractions[index]
        for index in range(component_count)
    )
    energy = (
        math.fsum(
            stream.molar_flow_kmol_s * stream.molar_enthalpy_j_per_kmol
            for stream in incoming
        )
        + float(heat_duty_w)
        - liquid.molar_flow_kmol_s * liquid.molar_enthalpy_j_per_kmol
        - vapor.molar_flow_kmol_s * vapor.molar_enthalpy_j_per_kmol
    )
    return StageResiduals(
        component_material_kmol_s=component_material,
        phase_equilibrium=phase_equilibrium,
        liquid_summation=math.fsum(liquid.mole_fractions) - 1.0,
        vapor_summation=math.fsum(vapor.mole_fractions) - 1.0,
        energy_w=energy,
    )


@dataclass(frozen=True, slots=True)
class ShortcutColumnResult:
    """Fenske-Underwood-Gilliland result on a coherent kmol/s basis."""

    distillate_flow_kmol_s: float
    bottoms_flow_kmol_s: float
    distillate_mole_fractions: tuple[float, ...]
    bottoms_mole_fractions: tuple[float, ...]
    minimum_stages: float
    minimum_reflux_ratio: float
    actual_stages: float
    feed_stage: float
    rectifying_liquid_flow_kmol_s: float
    stripping_liquid_flow_kmol_s: float
    rectifying_vapor_flow_kmol_s: float
    stripping_vapor_flow_kmol_s: float
    iterations: int


@dataclass(frozen=True, slots=True)
class ColumnBalanceResiduals:
    component_kmol_s: tuple[float, ...]
    total_kmol_s: float

    def is_closed(self, tolerance_kmol_s: float) -> bool:
        if not _finite_number(tolerance_kmol_s) or tolerance_kmol_s <= 0.0:
            raise ValidationError("column balance tolerance must be positive and finite")
        return (
            max(abs(value) for value in self.component_kmol_s) <= tolerance_kmol_s
            and abs(self.total_kmol_s) <= tolerance_kmol_s
        )


def column_balance_residuals(
    inlet_component_flows_kmol_s: tuple[tuple[float, ...], ...],
    outlet_component_flows_kmol_s: tuple[tuple[float, ...], ...],
) -> ColumnBalanceResiduals:
    """Return steady column component and total material residuals."""
    try:
        inlets = tuple(tuple(stream) for stream in inlet_component_flows_kmol_s)
        outlets = tuple(tuple(stream) for stream in outlet_component_flows_kmol_s)
    except TypeError as error:
        raise ValidationError("column balance streams must be finite sequences") from error
    if not inlets or not outlets or not inlets[0]:
        raise ValidationError("column balance requires inlet and outlet component flows")
    component_count = len(inlets[0])
    streams = inlets + outlets
    if any(
        len(stream) != component_count
        or any(not _finite_number(value) or value < 0.0 for value in stream)
        for stream in streams
    ):
        raise ValidationError("column balance component flows are invalid")
    residuals = tuple(
        math.fsum(stream[index] for stream in inlets)
        - math.fsum(stream[index] for stream in outlets)
        for index in range(component_count)
    )
    return ColumnBalanceResiduals(residuals, math.fsum(residuals))


@dataclass(frozen=True, slots=True)
class ColumnProfileResiduals:
    component_material_kmol_s: tuple[tuple[float, ...], ...]
    phase_equilibrium: tuple[tuple[float, ...], ...]
    liquid_summation: tuple[float, ...]
    vapor_summation: tuple[float, ...]
    flow_scale_kmol_s: float

    @property
    def maximum_scaled_component_residual(self) -> float:
        return max(
            abs(value) / self.flow_scale_kmol_s
            for row in self.component_material_kmol_s
            for value in row
        )

    def is_closed(
        self,
        scaled_component_tolerance: float,
        equilibrium_tolerance: float,
        summation_tolerance: float,
    ) -> bool:
        tolerances = (
            scaled_component_tolerance,
            equilibrium_tolerance,
            summation_tolerance,
        )
        if any(not _finite_number(value) or value <= 0.0 for value in tolerances):
            raise ValidationError("column profile tolerances must be positive and finite")
        return (
            self.maximum_scaled_component_residual <= scaled_component_tolerance
            and max(
                abs(value) for row in self.phase_equilibrium for value in row
            )
            <= equilibrium_tolerance
            and max(abs(value) for value in self.liquid_summation)
            <= summation_tolerance
            and max(abs(value) for value in self.vapor_summation)
            <= summation_tolerance
        )


def fixed_k_column_profile_residuals(
    feed_component_flows_by_stage_kmol_s: tuple[tuple[float, ...], ...],
    liquid_flows_kmol_s: tuple[float, ...],
    vapor_flows_kmol_s: tuple[float, ...],
    liquid_mole_fractions: tuple[tuple[float, ...], ...],
    vapor_mole_fractions: tuple[tuple[float, ...], ...],
    equilibrium_ratios_by_stage: tuple[tuple[float, ...], ...],
    liquid_product_flows_by_stage_kmol_s: tuple[float, ...] | None = None,
    vapor_product_flows_by_stage_kmol_s: tuple[float, ...] | None = None,
) -> ColumnProfileResiduals:
    """Evaluate fixed-K material, equilibrium, and summation profile closure.

    Internal liquid flows point toward increasing stage index and internal
    vapor flows point toward decreasing stage index. End-stage internal flows
    therefore represent the ordinary vapor-overhead and liquid-bottoms
    products. Explicit product arrays cover side draws and a total-condenser
    liquid distillate without changing that convention.
    """
    try:
        feeds = tuple(tuple(float(value) for value in row) for row in feed_component_flows_by_stage_kmol_s)
        liquid = tuple(float(value) for value in liquid_flows_kmol_s)
        vapor = tuple(float(value) for value in vapor_flows_kmol_s)
        liquid_fractions = tuple(tuple(float(value) for value in row) for row in liquid_mole_fractions)
        vapor_fractions = tuple(tuple(float(value) for value in row) for row in vapor_mole_fractions)
        ratios = tuple(tuple(float(value) for value in row) for row in equilibrium_ratios_by_stage)
        liquid_products = (
            (0.0,) * len(feeds)
            if liquid_product_flows_by_stage_kmol_s is None
            else tuple(float(value) for value in liquid_product_flows_by_stage_kmol_s)
        )
        vapor_products = (
            (0.0,) * len(feeds)
            if vapor_product_flows_by_stage_kmol_s is None
            else tuple(float(value) for value in vapor_product_flows_by_stage_kmol_s)
        )
    except (TypeError, ValueError) as error:
        raise ValidationError("column profile inputs must be finite sequences") from error
    stage_count = len(feeds)
    component_count = len(feeds[0]) if feeds else 0
    matrices = feeds + liquid_fractions + vapor_fractions + ratios
    if (
        stage_count < 2
        or component_count < 2
        or any(len(values) != stage_count for values in (
            liquid, vapor, liquid_products, vapor_products,
        ))
        or any(len(matrix) != stage_count for matrix in (
            liquid_fractions, vapor_fractions, ratios,
        ))
        or any(len(row) != component_count for row in matrices)
        or any(not _finite_number(value) or value < 0.0 for row in feeds for value in row)
        or any(not _finite_number(value) or value < 0.0 for value in (
            liquid + vapor + liquid_products + vapor_products
        ))
        or any(not _finite_number(value) or value < 0.0 for row in (
            liquid_fractions + vapor_fractions
        ) for value in row)
        or any(not _finite_number(value) or value <= 0.0 for row in ratios for value in row)
    ):
        raise ValidationError("column profile inputs are invalid")
    flow_scale = math.fsum(value for row in feeds for value in row)
    if flow_scale <= 0.0:
        raise ValidationError("column profile requires at least one feed")

    material: list[tuple[float, ...]] = []
    for stage in range(stage_count):
        stage_residuals = []
        for component in range(component_count):
            residual = feeds[stage][component]
            if stage > 0:
                residual += liquid[stage - 1] * liquid_fractions[stage - 1][component]
            if stage < stage_count - 1:
                residual += vapor[stage + 1] * vapor_fractions[stage + 1][component]
            residual -= (
                liquid[stage] + liquid_products[stage]
            ) * liquid_fractions[stage][component]
            residual -= (
                vapor[stage] + vapor_products[stage]
            ) * vapor_fractions[stage][component]
            stage_residuals.append(residual)
        material.append(tuple(stage_residuals))
    equilibrium = tuple(
        tuple(
            vapor_fraction - ratio * liquid_fraction
            for vapor_fraction, ratio, liquid_fraction in zip(
                vapor_row, ratio_row, liquid_row
            )
        )
        for vapor_row, ratio_row, liquid_row in zip(
            vapor_fractions, ratios, liquid_fractions
        )
    )
    return ColumnProfileResiduals(
        tuple(material),
        equilibrium,
        tuple(math.fsum(row) - 1.0 for row in liquid_fractions),
        tuple(math.fsum(row) - 1.0 for row in vapor_fractions),
        flow_scale,
    )


def column_energy_residual_w(
    inlet_mass_enthalpy_kj_s: tuple[tuple[float, float], ...],
    outlet_mass_enthalpy_kj_s: tuple[tuple[float, float], ...],
    heat_input_w: float = 0.0,
    heat_output_w: float = 0.0,
) -> float:
    """Return steady column energy residual from kg/s and kJ/kg streams."""
    try:
        inlets = tuple((float(flow), float(enthalpy)) for flow, enthalpy in inlet_mass_enthalpy_kj_s)
        outlets = tuple((float(flow), float(enthalpy)) for flow, enthalpy in outlet_mass_enthalpy_kj_s)
    except (TypeError, ValueError) as error:
        raise ValidationError("column energy streams must be flow-enthalpy pairs") from error
    if (
        not inlets
        or not outlets
        or any(
            not _finite_number(value)
            for pair in inlets + outlets
            for value in pair
        )
        or any(flow < 0.0 for flow, _ in inlets + outlets)
        or not _finite_number(heat_input_w)
        or not _finite_number(heat_output_w)
    ):
        raise ValidationError("column energy inputs are invalid")
    return (
        1000.0 * math.fsum(flow * enthalpy for flow, enthalpy in inlets)
        + float(heat_input_w)
        - 1000.0 * math.fsum(flow * enthalpy for flow, enthalpy in outlets)
        - float(heat_output_w)
    )


def column_profile_energy_residuals(
    feed_energy_by_stage_w: tuple[float, ...],
    liquid_flows_kmol_s: tuple[float, ...],
    vapor_flows_kmol_s: tuple[float, ...],
    liquid_molar_enthalpies_j_kmol: tuple[float, ...],
    vapor_molar_enthalpies_j_kmol: tuple[float, ...],
    heat_duties_by_stage_w: tuple[float, ...] | None = None,
    liquid_product_flows_by_stage_kmol_s: tuple[float, ...] | None = None,
    vapor_product_flows_by_stage_kmol_s: tuple[float, ...] | None = None,
) -> tuple[float, ...]:
    """Evaluate the steady energy balance on every column stage.

    Internal liquid flows point toward increasing stage index and internal
    vapor flows point toward decreasing stage index, matching
    :func:`fixed_k_column_profile_residuals`. Positive stage duties add heat.
    Molar flows and molar enthalpies use one coherent kmol basis, so their
    products are watts.
    """
    try:
        feeds = tuple(float(value) for value in feed_energy_by_stage_w)
        liquid = tuple(float(value) for value in liquid_flows_kmol_s)
        vapor = tuple(float(value) for value in vapor_flows_kmol_s)
        liquid_enthalpy = tuple(float(value) for value in liquid_molar_enthalpies_j_kmol)
        vapor_enthalpy = tuple(float(value) for value in vapor_molar_enthalpies_j_kmol)
        stage_count = len(feeds)
        duties = (
            (0.0,) * stage_count
            if heat_duties_by_stage_w is None
            else tuple(float(value) for value in heat_duties_by_stage_w)
        )
        liquid_products = (
            (0.0,) * stage_count
            if liquid_product_flows_by_stage_kmol_s is None
            else tuple(float(value) for value in liquid_product_flows_by_stage_kmol_s)
        )
        vapor_products = (
            (0.0,) * stage_count
            if vapor_product_flows_by_stage_kmol_s is None
            else tuple(float(value) for value in vapor_product_flows_by_stage_kmol_s)
        )
    except (TypeError, ValueError) as error:
        raise ValidationError("column stage energy inputs must be finite sequences") from error
    vectors = (
        feeds,
        liquid,
        vapor,
        liquid_enthalpy,
        vapor_enthalpy,
        duties,
        liquid_products,
        vapor_products,
    )
    if (
        stage_count < 2
        or any(len(values) != stage_count for values in vectors)
        or any(not _finite_number(value) for values in vectors for value in values)
        or any(
            value < 0.0
            for values in (liquid, vapor, liquid_products, vapor_products)
            for value in values
        )
    ):
        raise ValidationError("column stage energy inputs are invalid")

    residuals = []
    for stage in range(stage_count):
        residual = feeds[stage] + duties[stage]
        if stage > 0:
            residual += liquid[stage - 1] * liquid_enthalpy[stage - 1]
        if stage < stage_count - 1:
            residual += vapor[stage + 1] * vapor_enthalpy[stage + 1]
        residual -= (
            liquid[stage] + liquid_products[stage]
        ) * liquid_enthalpy[stage]
        residual -= (
            vapor[stage] + vapor_products[stage]
        ) * vapor_enthalpy[stage]
        residuals.append(residual)
    return tuple(residuals)


@dataclass(frozen=True, slots=True)
class SumRatesIteration:
    iteration: int
    maximum_vapor_flow_change_kmol_s: float


@dataclass(frozen=True, slots=True)
class SumRatesProfile:
    liquid_flows_kmol_s: tuple[float, ...]
    vapor_flows_kmol_s: tuple[float, ...]
    liquid_mole_fractions: tuple[tuple[float, ...], ...]
    vapor_mole_fractions: tuple[tuple[float, ...], ...]
    history: tuple[SumRatesIteration, ...]


class ColumnConvergenceError(RuntimeError):
    def __init__(self, message: str, history: tuple[SumRatesIteration, ...]):
        super().__init__(message)
        self.history = history


@dataclass(frozen=True, slots=True)
class ColumnNewtonIteration:
    iteration: int
    scaled_residual_norm: float
    damping: float


@dataclass(frozen=True, slots=True)
class FixedKColumnProfile:
    liquid_flows_kmol_s: tuple[float, ...]
    vapor_flows_kmol_s: tuple[float, ...]
    liquid_mole_fractions: tuple[tuple[float, ...], ...]
    vapor_mole_fractions: tuple[tuple[float, ...], ...]
    history: tuple[ColumnNewtonIteration, ...]


class ColumnNewtonConvergenceError(RuntimeError):
    def __init__(self, message: str, history: tuple[ColumnNewtonIteration, ...]):
        super().__init__(message)
        self.history = history


def _dense_linear_solve(matrix: list[list[float]], right: list[float]) -> list[float]:
    """Solve a dense square system with scaled partial pivoting."""
    count = len(right)
    augmented = [row[:] + [value] for row, value in zip(matrix, right)]
    scales = [max(abs(value) for value in row) for row in matrix]
    if any(scale == 0.0 or not math.isfinite(scale) for scale in scales):
        raise ValidationError("column Newton Jacobian is singular")
    for column in range(count):
        pivot = max(
            range(column, count),
            key=lambda row: abs(augmented[row][column]) / scales[row],
        )
        if abs(augmented[pivot][column]) <= 1.0e-14 * scales[pivot]:
            raise ValidationError("column Newton Jacobian is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scales[column], scales[pivot] = scales[pivot], scales[column]
        divisor = augmented[column][column]
        for row in range(column + 1, count):
            factor = augmented[row][column] / divisor
            if factor == 0.0:
                continue
            augmented[row][column] = 0.0
            for index in range(column + 1, count + 1):
                augmented[row][index] -= factor * augmented[column][index]
    solution = [0.0] * count
    for row in range(count - 1, -1, -1):
        value = augmented[row][-1] - math.fsum(
            augmented[row][column] * solution[column]
            for column in range(row + 1, count)
        )
        solution[row] = value / augmented[row][row]
    if any(not math.isfinite(value) for value in solution):
        raise ValidationError("column Newton direction is not finite")
    return solution


def fixed_k_material_column(
    feed_component_flows_by_stage_kmol_s: tuple[tuple[float, ...], ...],
    equilibrium_ratios_by_stage: tuple[tuple[float, ...], ...],
    initial_liquid_flows_kmol_s: tuple[float, ...],
    initial_vapor_flows_kmol_s: tuple[float, ...],
    initial_liquid_mole_fractions: tuple[tuple[float, ...], ...],
    residual_tolerance: float = 1.0e-10,
    maximum_iterations: int = 30,
    jacobian_step: float = 1.0e-6,
    minimum_damping: float = 2.0**-24,
    liquid_product_flows_by_stage_kmol_s: tuple[float, ...] | None = None,
    vapor_product_flows_by_stage_kmol_s: tuple[float, ...] | None = None,
) -> FixedKColumnProfile:
    """Close fixed-K stage material balances and both phase summations.

    Positive phase flows are represented logarithmically and each liquid
    composition uses reduced softmax coordinates. A zero initial vapor flow
    fixes that phase flow at zero, as required by a total condenser, and omits
    its otherwise redundant vapor-summation equation. Vapor compositions
    follow ``y = K*x`` on active vapor stages.
    """
    try:
        feeds = tuple(
            tuple(float(value) for value in row)
            for row in feed_component_flows_by_stage_kmol_s
        )
        ratios = tuple(
            tuple(float(value) for value in row)
            for row in equilibrium_ratios_by_stage
        )
        initial_liquid = tuple(float(value) for value in initial_liquid_flows_kmol_s)
        initial_vapor = tuple(float(value) for value in initial_vapor_flows_kmol_s)
        initial_fractions = tuple(
            tuple(float(value) for value in row) for row in initial_liquid_mole_fractions
        )
        liquid_products = (
            (0.0,) * len(feeds)
            if liquid_product_flows_by_stage_kmol_s is None
            else tuple(float(value) for value in liquid_product_flows_by_stage_kmol_s)
        )
        vapor_products = (
            (0.0,) * len(feeds)
            if vapor_product_flows_by_stage_kmol_s is None
            else tuple(float(value) for value in vapor_product_flows_by_stage_kmol_s)
        )
    except (TypeError, ValueError) as error:
        raise ValidationError("fixed-K column inputs must be finite sequences") from error
    stage_count = len(feeds)
    component_count = len(feeds[0]) if feeds else 0
    scalar_controls = (residual_tolerance, jacobian_step, minimum_damping)
    if (
        stage_count < 2
        or component_count < 2
        or len(ratios) != stage_count
        or len(initial_liquid) != stage_count
        or len(initial_vapor) != stage_count
        or len(initial_fractions) != stage_count
        or len(liquid_products) != stage_count
        or len(vapor_products) != stage_count
        or any(
            len(row) != component_count
            for row in feeds + ratios + initial_fractions
        )
        or any(not _finite_number(value) or value < 0.0 for row in feeds for value in row)
        or any(not _finite_number(value) or value <= 0.0 for row in ratios for value in row)
        or any(not _finite_number(value) or value <= 0.0 for value in initial_liquid)
        or any(not _finite_number(value) or value < 0.0 for value in initial_vapor)
        or any(
            not _finite_number(value) or value < 0.0
            for value in liquid_products + vapor_products
        )
        or any(
            not _finite_number(value) or value <= 0.0
            for row in initial_fractions
            for value in row
        )
        or any(abs(math.fsum(row) - 1.0) > 1.0e-8 for row in initial_fractions)
        or any(not _finite_number(value) or value <= 0.0 for value in scalar_controls)
        or minimum_damping > 1.0
        or isinstance(maximum_iterations, bool)
        or not isinstance(maximum_iterations, int)
        or maximum_iterations <= 0
    ):
        raise ValidationError("fixed-K column inputs are invalid")
    flow_scale = math.fsum(value for row in feeds for value in row)
    if flow_scale <= 0.0:
        raise ValidationError("fixed-K column requires at least one feed")

    active_vapor_stages = tuple(
        stage for stage, value in enumerate(initial_vapor) if value > 0.0
    )
    flow_variable_count = stage_count + len(active_vapor_stages)
    variables = (
        [math.log(value) for value in initial_liquid]
        + [math.log(initial_vapor[stage]) for stage in active_vapor_stages]
        + [
            math.log(row[component] / row[-1])
            for row in initial_fractions
            for component in range(component_count - 1)
        ]
    )

    def decode(
        values: list[float],
    ) -> tuple[list[float], list[float], list[list[float]]]:
        if any(value < -700.0 or value > 700.0 for value in values[:flow_variable_count]):
            raise ValidationError(
                "fixed-K column trial flows are outside the representable range"
            )
        liquid = [math.exp(value) for value in values[:stage_count]]
        vapor = [0.0] * stage_count
        for index, stage in enumerate(active_vapor_stages, start=stage_count):
            vapor[stage] = math.exp(values[index])
        liquid_fractions: list[list[float]] = []
        offset = flow_variable_count
        reduced_count = component_count - 1
        for stage in range(stage_count):
            start = offset + stage * reduced_count
            logits = values[start : start + reduced_count] + [0.0]
            maximum = max(logits)
            exponentials = [math.exp(value - maximum) for value in logits]
            total = math.fsum(exponentials)
            liquid_fractions.append([value / total for value in exponentials])
        return liquid, vapor, liquid_fractions

    def evaluate(
        values: list[float],
    ) -> tuple[list[float], list[float], list[list[float]], list[list[float]]]:
        liquid, vapor, liquid_fractions = decode(values)
        vapor_fractions = [
            [ratio * fraction for ratio, fraction in zip(ratio_row, fraction_row)]
            for ratio_row, fraction_row in zip(ratios, liquid_fractions)
        ]
        residuals: list[float] = []
        for stage in range(stage_count):
            for component in range(component_count):
                residual = feeds[stage][component]
                if stage > 0:
                    residual += liquid[stage - 1] * liquid_fractions[stage - 1][component]
                if stage < stage_count - 1:
                    residual += vapor[stage + 1] * vapor_fractions[stage + 1][component]
                residual -= (
                    liquid[stage] + liquid_products[stage]
                ) * liquid_fractions[stage][component]
                residual -= (
                    vapor[stage] + vapor_products[stage]
                ) * vapor_fractions[stage][component]
                residuals.append(residual / flow_scale)
        residuals.extend(
            math.fsum(vapor_fractions[stage]) - 1.0
            for stage in active_vapor_stages
        )
        return residuals, liquid, liquid_fractions, vapor_fractions

    history: list[ColumnNewtonIteration] = []
    for iteration in range(1, maximum_iterations + 1):
        current, _, _, _ = evaluate(variables)
        norm = max(abs(value) for value in current)
        if norm <= residual_tolerance:
            liquid, vapor, liquid_fractions = decode(variables)
            vapor_fractions = [
                [ratio * fraction for ratio, fraction in zip(ratio_row, fraction_row)]
                for ratio_row, fraction_row in zip(ratios, liquid_fractions)
            ]
            return FixedKColumnProfile(
                tuple(liquid),
                tuple(vapor),
                tuple(tuple(row) for row in liquid_fractions),
                tuple(tuple(row) for row in vapor_fractions),
                tuple(history),
            )
        size = len(variables)
        jacobian = [[0.0] * size for _ in range(size)]
        for column in range(size):
            step = jacobian_step * max(1.0, abs(variables[column]))
            trial = variables[:]
            trial[column] += step
            shifted, _, _, _ = evaluate(trial)
            for row in range(size):
                jacobian[row][column] = (shifted[row] - current[row]) / step
        try:
            direction = _dense_linear_solve(jacobian, [-value for value in current])
        except ValidationError as error:
            raise ColumnNewtonConvergenceError(str(error), tuple(history)) from error
        damping = 1.0
        accepted = False
        while damping >= minimum_damping:
            trial = [value + damping * delta for value, delta in zip(variables, direction)]
            try:
                shifted, _, _, _ = evaluate(trial)
            except ValidationError:
                damping *= 0.5
                continue
            shifted_norm = max(abs(value) for value in shifted)
            if shifted_norm < norm:
                variables = trial
                history.append(ColumnNewtonIteration(iteration, shifted_norm, damping))
                accepted = True
                break
            damping *= 0.5
        if not accepted:
            raise ColumnNewtonConvergenceError(
                "fixed-K column line search did not reduce the residual", tuple(history)
            )
    final_residuals, final_liquid, final_liquid_fractions, final_vapor_fractions = evaluate(
        variables
    )
    if max(abs(value) for value in final_residuals) <= residual_tolerance:
        _, final_vapor, _ = decode(variables)
        return FixedKColumnProfile(
            tuple(final_liquid),
            tuple(final_vapor),
            tuple(tuple(row) for row in final_liquid_fractions),
            tuple(tuple(row) for row in final_vapor_fractions),
            tuple(history),
        )
    raise ColumnNewtonConvergenceError(
        "fixed-K column did not converge", tuple(history)
    )


def _tridiagonal_solve(
    lower: list[float], diagonal: list[float], upper: list[float], right: list[float],
) -> list[float]:
    count = len(right)
    upper_work = [0.0] * count
    right_work = [0.0] * count
    if abs(diagonal[0]) <= 1.0e-300:
        raise ValidationError("column tridiagonal matrix is singular")
    upper_work[0] = upper[0] / diagonal[0]
    right_work[0] = right[0] / diagonal[0]
    for index in range(1, count):
        denominator = diagonal[index] - lower[index] * upper_work[index - 1]
        if abs(denominator) <= 1.0e-300 or not math.isfinite(denominator):
            raise ValidationError("column tridiagonal matrix is singular")
        if index < count - 1:
            upper_work[index] = upper[index] / denominator
        right_work[index] = (
            right[index] - lower[index] * right_work[index - 1]
        ) / denominator
    result = [0.0] * count
    result[-1] = right_work[-1]
    for index in range(count - 2, -1, -1):
        result[index] = right_work[index] - upper_work[index] * result[index + 1]
    return result


def fixed_k_sum_rates_absorber(
    feed_component_flows_by_stage_kmol_s: tuple[tuple[float, ...], ...],
    equilibrium_ratios_by_stage: tuple[tuple[float, ...], ...],
    initial_liquid_flows_kmol_s: tuple[float, ...],
    initial_vapor_flows_kmol_s: tuple[float, ...],
    flow_tolerance_kmol_s: float = 1.0e-12,
    maximum_iterations: int = 200,
    component_floor_kmol_s: float = 1.0e-10,
) -> SumRatesProfile:
    """Solve the fixed-K material loop of DWSIM's Burningham–Otto method."""
    try:
        feeds = tuple(tuple(row) for row in feed_component_flows_by_stage_kmol_s)
        ratios = tuple(tuple(row) for row in equilibrium_ratios_by_stage)
        liquid = [float(value) for value in initial_liquid_flows_kmol_s]
        vapor = [float(value) for value in initial_vapor_flows_kmol_s]
    except (TypeError, ValueError) as error:
        raise ValidationError("sum-rates inputs must be finite sequences") from error
    stage_count = len(feeds)
    component_count = len(feeds[0]) if feeds else 0
    if (
        stage_count < 2
        or component_count < 2
        or len(ratios) != stage_count
        or len(liquid) != stage_count
        or len(vapor) != stage_count
        or any(len(row) != component_count for row in feeds + ratios)
        or any(not _finite_number(value) or value < 0.0 for row in feeds for value in row)
        or any(not _finite_number(value) or value <= 0.0 for row in ratios for value in row)
        or any(not _finite_number(value) or value <= 0.0 for value in liquid + vapor)
        or not _finite_number(flow_tolerance_kmol_s)
        or flow_tolerance_kmol_s <= 0.0
        or not _finite_number(component_floor_kmol_s)
        or component_floor_kmol_s <= 0.0
        or isinstance(maximum_iterations, bool)
        or not isinstance(maximum_iterations, int)
        or maximum_iterations <= 0
    ):
        raise ValidationError("sum-rates absorber inputs are invalid")

    feed_flows = [math.fsum(row) for row in feeds]
    if math.fsum(feed_flows) <= 0.0:
        raise ValidationError("sum-rates absorber requires at least one feed")
    feed_fractions = [
        tuple(value / flow for value in row) if flow > 0.0 else (0.0,) * component_count
        for row, flow in zip(feeds, feed_flows)
    ]
    history: list[SumRatesIteration] = []
    liquid_fractions: list[list[float]] = []
    vapor_fractions: list[list[float]] = []

    for iteration in range(1, maximum_iterations + 1):
        cumulative_feed = [math.fsum(feed_flows[: index + 1]) for index in range(stage_count)]
        preceding_feed = [math.fsum(feed_flows[:index]) for index in range(stage_count)]
        liquid_components = [[0.0] * component_count for _ in range(stage_count)]
        for component in range(component_count):
            lower = [0.0] * stage_count
            diagonal = [0.0] * stage_count
            upper = [0.0] * stage_count
            right = [
                -feed_flows[stage] * feed_fractions[stage][component]
                for stage in range(stage_count)
            ]
            for stage in range(stage_count):
                next_vapor = vapor[stage + 1] if stage < stage_count - 1 else 0.0
                diagonal[stage] = -(
                    next_vapor
                    + cumulative_feed[stage]
                    - vapor[0]
                    + vapor[stage] * ratios[stage][component]
                )
                if stage < stage_count - 1:
                    upper[stage] = vapor[stage + 1] * ratios[stage + 1][component]
                if stage > 0:
                    lower[stage] = vapor[stage] + preceding_feed[stage] - vapor[0]
            solution = _tridiagonal_solve(lower, diagonal, upper, right)
            for stage, value in enumerate(solution):
                liquid_components[stage][component] = max(value, component_floor_kmol_s)

        liquid_sums = [math.fsum(row) for row in liquid_components]
        liquid = [flow * total for flow, total in zip(liquid, liquid_sums)]
        liquid_fractions = [
            [value / total for value in row]
            for row, total in zip(liquid_components, liquid_sums)
        ]
        vapor_fractions = []
        for stage in range(stage_count):
            trial = [
                liquid_fractions[stage][component] * ratios[stage][component]
                for component in range(component_count)
            ]
            total = math.fsum(trial)
            vapor_fractions.append([value / total for value in trial])

        following_feed = [math.fsum(feed_flows[index:]) for index in range(stage_count)]
        previous_vapor = vapor
        vapor = [
            abs(
                (liquid[stage - 1] if stage > 0 else 0.0)
                - liquid[-1]
                + following_feed[stage]
            )
            for stage in range(stage_count)
        ]
        maximum_change = max(
            abs(value - previous) for value, previous in zip(vapor, previous_vapor)
        )
        history.append(SumRatesIteration(iteration, maximum_change))
        if maximum_change <= flow_tolerance_kmol_s:
            return SumRatesProfile(
                tuple(liquid),
                tuple(vapor),
                tuple(tuple(row) for row in liquid_fractions),
                tuple(tuple(row) for row in vapor_fractions),
                tuple(history),
            )

    raise ColumnConvergenceError(
        "fixed-K sum-rates absorber did not converge", tuple(history)
    )


def _underwood_value(
    theta: float,
    relative_volatilities: tuple[float, ...],
    feed_mole_fractions: tuple[float, ...],
    feed_liquid_fraction: float,
) -> float:
    return (
        math.fsum(
            alpha * fraction / (alpha - theta)
            for alpha, fraction in zip(relative_volatilities, feed_mole_fractions)
            if fraction != 0.0
        )
        - 1.0
        + feed_liquid_fraction
    )


def _bounded_underwood_minimum(
    lower: float,
    upper: float,
    relative_volatilities: tuple[float, ...],
    feed_mole_fractions: tuple[float, ...],
    feed_liquid_fraction: float,
) -> float:
    """Minimize DWSIM's squared Underwood residual on one pole-free interval."""
    upper = math.nextafter(upper, lower)
    if not lower < upper:
        raise ValidationError("shortcut-column Underwood interval is invalid")

    def objective(theta: float) -> float:
        value = _underwood_value(
            theta, relative_volatilities, feed_mole_fractions, feed_liquid_fraction,
        )
        return value * value

    # A bounded golden-section search also handles DWSIM cases where the
    # constrained minimum lies at the 1.01 * heavy-key endpoint.
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = lower
    right = upper
    first = right - ratio * (right - left)
    second = left + ratio * (right - left)
    first_value = objective(first)
    second_value = objective(second)
    for _ in range(160):
        if right - left <= 1.0e-12 * max(1.0, abs(left), abs(right)):
            break
        if first_value <= second_value:
            right = second
            second = first
            second_value = first_value
            first = right - ratio * (right - left)
            first_value = objective(first)
        else:
            left = first
            first = second
            first_value = second_value
            second = left + ratio * (right - left)
            second_value = objective(second)
    candidates = (lower, first, second, right)
    return min(candidates, key=objective)


def shortcut_column(
    feed_flow_kmol_s: float,
    feed_mole_fractions: tuple[float, ...],
    relative_volatilities: tuple[float, ...],
    light_key_index: int,
    heavy_key_index: int,
    heavy_key_distillate_fraction: float,
    light_key_bottoms_fraction: float,
    reflux_ratio: float,
    feed_liquid_fraction: float,
    convergence_tolerance: float = 1.0e-4,
    maximum_iterations: int = 100,
) -> ShortcutColumnResult:
    """Solve DWSIM-compatible Fenske-Underwood-Gilliland shortcut balances.

    Relative volatilities are referenced internally to the heavy key. The
    current slice supports the ordinary Underwood root with no distributed
    non-key roots; those are detected explicitly instead of silently applying
    the wrong minimum-reflux equation.
    """
    try:
        composition = tuple(feed_mole_fractions)
        alpha_input = tuple(relative_volatilities)
    except TypeError as error:
        raise ValidationError("shortcut-column composition and volatility must be sequences") from error
    component_count = len(composition)
    scalars = (
        feed_flow_kmol_s,
        heavy_key_distillate_fraction,
        light_key_bottoms_fraction,
        reflux_ratio,
        feed_liquid_fraction,
        convergence_tolerance,
    )
    if (
        component_count < 2
        or len(alpha_input) != component_count
        or any(not _finite_number(value) or value < 0.0 for value in composition)
        or math.fsum(composition) <= 0.0
        or any(not _finite_number(value) or value <= 0.0 for value in alpha_input)
        or any(not _finite_number(value) for value in scalars)
        or feed_flow_kmol_s <= 0.0
        or not 0 <= light_key_index < component_count
        or not 0 <= heavy_key_index < component_count
        or light_key_index == heavy_key_index
        or not 0.0 < heavy_key_distillate_fraction < 1.0
        or not 0.0 < light_key_bottoms_fraction < 1.0
        or reflux_ratio < 0.0
        or convergence_tolerance <= 0.0
        or isinstance(maximum_iterations, bool)
        or not isinstance(maximum_iterations, int)
        or maximum_iterations <= 0
    ):
        raise ValidationError("shortcut-column inputs are invalid")

    heavy_alpha = alpha_input[heavy_key_index]
    alpha = tuple(value / heavy_alpha for value in alpha_input)
    if alpha[light_key_index] <= 1.0:
        raise ValidationError("shortcut-column light key must be more volatile than the heavy key")

    distillate_flow = feed_flow_kmol_s * composition[light_key_index]
    distillate_flow += math.fsum(
        feed_flow_kmol_s * composition[index]
        for index, value in enumerate(alpha)
        if value > alpha[light_key_index]
    )
    distillate = [0.0] * component_count
    bottoms = [0.0] * component_count
    minimum_stages = math.nan
    bottoms_flow = math.nan

    for iteration in range(1, maximum_iterations + 1):
        bottoms_flow = feed_flow_kmol_s - distillate_flow
        if bottoms_flow <= 0.0 or distillate_flow <= 0.0:
            raise ValidationError("shortcut-column product flow estimate is invalid")
        distillate[heavy_key_index] = heavy_key_distillate_fraction
        bottoms[light_key_index] = light_key_bottoms_fraction
        bottoms[heavy_key_index] = (
            feed_flow_kmol_s * composition[heavy_key_index]
            - distillate_flow * distillate[heavy_key_index]
        ) / bottoms_flow
        distillate[light_key_index] = (
            feed_flow_kmol_s * composition[light_key_index]
            - bottoms_flow * bottoms[light_key_index]
        ) / distillate_flow
        key_values = (
            bottoms[heavy_key_index], distillate[light_key_index],
        )
        if any(not _finite_number(value) or value <= 0.0 for value in key_values):
            raise ValidationError("shortcut-column key specifications are infeasible")

        separation = (
            distillate[light_key_index]
            / distillate[heavy_key_index]
            * bottoms[heavy_key_index]
            / bottoms[light_key_index]
        )
        minimum_stages = math.log(separation) / math.log(alpha[light_key_index])
        constant = math.log10(
            distillate[heavy_key_index] / bottoms[heavy_key_index]
        )
        for index in range(component_count):
            if index in (light_key_index, heavy_key_index):
                continue
            distribution = 10.0 ** (
                minimum_stages * math.log10(alpha[index]) + constant
            )
            bottoms[index] = (
                feed_flow_kmol_s * composition[index]
                / (bottoms_flow + distillate_flow * distribution)
            )
            distillate[index] = bottoms[index] * distribution

        previous_distillate_flow = distillate_flow
        distillate_flow = previous_distillate_flow * math.fsum(
            fraction
            for fraction, feed_fraction in zip(distillate, composition)
            if feed_fraction != 0.0
        )
        if not _finite_number(distillate_flow) or distillate_flow <= 0.0:
            raise ValidationError("shortcut-column distillate iteration failed")
        if abs((distillate_flow - previous_distillate_flow) / distillate_flow) < convergence_tolerance:
            break
    else:
        raise ValidationError("shortcut-column product split did not converge")

    distributed_nonkeys = []
    for index, value in enumerate(alpha):
        distribution_fraction = (
            (value - 1.0)
            / (alpha[light_key_index] - 1.0)
            * distillate_flow
            * distillate[light_key_index]
            / (feed_flow_kmol_s * composition[light_key_index])
            + (alpha[light_key_index] - value)
            / (alpha[light_key_index] - 1.0)
            * distillate_flow
            * distillate[heavy_key_index]
            / (feed_flow_kmol_s * composition[heavy_key_index])
        )
        if (
            index not in (light_key_index, heavy_key_index)
            and composition[index] != 0.0
            and 0.0 < distribution_fraction < 1.0
        ):
            distributed_nonkeys.append(index)
    if distributed_nonkeys:
        raise ValidationError("distributed non-key Underwood roots are not implemented")

    theta = _bounded_underwood_minimum(
        1.01,
        alpha[light_key_index],
        alpha,
        composition,
        feed_liquid_fraction,
    )
    minimum_reflux = math.fsum(
        value * fraction / (value - theta)
        for value, fraction, feed_fraction in zip(alpha, distillate, composition)
        if feed_fraction != 0.0
    ) - 1.0
    if minimum_reflux > reflux_ratio:
        raise ValidationError("shortcut-column reflux ratio is below the minimum")

    gilliland_x = (reflux_ratio - minimum_reflux) / (reflux_ratio + 1.0)
    gilliland_y = 0.75 * (1.0 - gilliland_x ** 0.5668)
    actual_stages = (gilliland_y + minimum_stages) / (1.0 - gilliland_y)
    stripping_separation = (
        composition[light_key_index]
        / composition[heavy_key_index]
        * bottoms[heavy_key_index]
        / bottoms[light_key_index]
    )
    stripping_minimum_stages = math.log(stripping_separation) / math.log(
        alpha[light_key_index]
    )
    feed_stage = stripping_minimum_stages * actual_stages / minimum_stages

    rectifying_liquid = reflux_ratio * distillate_flow
    stripping_liquid = rectifying_liquid + feed_liquid_fraction * feed_flow_kmol_s
    stripping_vapor = stripping_liquid - bottoms_flow
    rectifying_vapor = distillate_flow + rectifying_liquid
    if stripping_liquid < 0.0 or stripping_vapor < 0.0:
        raise ValidationError("shortcut-column internal flow is invalid")

    return ShortcutColumnResult(
        distillate_flow_kmol_s=distillate_flow,
        bottoms_flow_kmol_s=bottoms_flow,
        distillate_mole_fractions=tuple(distillate),
        bottoms_mole_fractions=tuple(bottoms),
        minimum_stages=minimum_stages,
        minimum_reflux_ratio=minimum_reflux,
        actual_stages=actual_stages,
        feed_stage=feed_stage,
        rectifying_liquid_flow_kmol_s=rectifying_liquid,
        stripping_liquid_flow_kmol_s=stripping_liquid,
        rectifying_vapor_flow_kmol_s=rectifying_vapor,
        stripping_vapor_flow_kmol_s=stripping_vapor,
        iterations=iteration,
    )
