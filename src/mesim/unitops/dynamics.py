"""Deterministic dynamic-model primitives on explicit physical bases.

The adaptive integration entry point intentionally accepts only an explicit
ODE.  Algebraic constraints can be described and checked here, but a model
that has not been reduced to an ODE requires an IDA-capable DAE solver.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

from scipy.integrate import solve_ivp

from ..errors import ConvergenceError, ValidationError


Vector = tuple[float, ...]
ExplicitRhs = Callable[[float, Vector], Sequence[float]]


def _finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValidationError(f"{name} must be finite")
    return float(value)


def _positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise ValidationError(f"{name} must be positive")
    return result


def _extended_bound(value: float, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or math.isnan(value)
    ):
        raise ValidationError(f"{name} must be numeric and not NaN")
    return float(value)


def _finite_vector(values: Sequence[float], name: str, *, non_negative: bool = False) -> Vector:
    try:
        result = tuple(_finite(value, name) for value in values)
    except TypeError as error:
        raise ValidationError(f"{name} must be a finite sequence") from error
    if not result:
        raise ValidationError(f"{name} must not be empty")
    if non_negative and any(value < 0.0 for value in result):
        raise ValidationError(f"{name} must be non-negative")
    return result


@dataclass(frozen=True, slots=True)
class HoldupState:
    """Component amounts in kmol and total internal energy in J."""

    component_amounts_kmol: Vector
    internal_energy_j: float


@dataclass(frozen=True, slots=True)
class HoldupRates:
    """Rates on the same kmol/s and W bases as :class:`HoldupState`."""

    component_in_kmol_s: Vector
    component_out_kmol_s: Vector
    component_generation_kmol_s: Vector
    energy_in_w: float
    energy_out_w: float
    heat_input_w: float = 0.0
    energy_generation_w: float = 0.0


@dataclass(frozen=True, slots=True)
class HoldupStepResult:
    state: HoldupState
    component_balance_residuals_kmol: Vector
    energy_balance_residual_j: float


def advance_holdup(state: HoldupState, rates: HoldupRates, step_s: float) -> HoldupStepResult:
    """Advance a lumped holdup with explicit Euler and audit the exact step balance."""
    amounts = _finite_vector(state.component_amounts_kmol, "component amounts", non_negative=True)
    internal_energy = _finite(state.internal_energy_j, "internal energy")
    inlet = _finite_vector(rates.component_in_kmol_s, "component inlet rates", non_negative=True)
    outlet = _finite_vector(rates.component_out_kmol_s, "component outlet rates", non_negative=True)
    generation = _finite_vector(rates.component_generation_kmol_s, "component generation rates")
    if not (len(amounts) == len(inlet) == len(outlet) == len(generation)):
        raise ValidationError("holdup component vectors must have the same length")
    duration = _positive(step_s, "holdup step")
    energy_rate = (
        _finite(rates.energy_in_w, "inlet energy rate")
        - _finite(rates.energy_out_w, "outlet energy rate")
        + _finite(rates.heat_input_w, "heat input")
        + _finite(rates.energy_generation_w, "energy generation")
    )
    component_rates = tuple(
        entering - leaving + produced
        for entering, leaving, produced in zip(inlet, outlet, generation, strict=True)
    )
    next_amounts = tuple(amount + duration * rate for amount, rate in zip(amounts, component_rates, strict=True))
    if any(amount < -1.0e-12 for amount in next_amounts):
        raise ValidationError("holdup step produced a negative component amount")
    next_amounts = tuple(max(amount, 0.0) for amount in next_amounts)
    next_energy = internal_energy + duration * energy_rate
    if not math.isfinite(next_energy):
        raise ValidationError("holdup step produced an unrepresentable energy state")
    component_residuals = tuple(
        new - old - duration * rate
        for new, old, rate in zip(next_amounts, amounts, component_rates, strict=True)
    )
    energy_residual = next_energy - internal_energy - duration * energy_rate
    return HoldupStepResult(
        HoldupState(next_amounts, next_energy), component_residuals, energy_residual,
    )


@dataclass(frozen=True, slots=True)
class AlgebraicConstraint:
    name: str
    residual: float
    scale: float = 1.0


def validate_algebraic_constraints(
    constraints: Sequence[AlgebraicConstraint], tolerance: float,
) -> tuple[AlgebraicConstraint, ...]:
    """Validate that an initialized state is consistent with its algebraic equations."""
    limit = _positive(tolerance, "algebraic tolerance")
    try:
        records = tuple(constraints)
    except TypeError as error:
        raise ValidationError("algebraic constraints must be a finite sequence") from error
    if not records:
        raise ValidationError("at least one algebraic constraint is required")
    for record in records:
        if not isinstance(record, AlgebraicConstraint) or not record.name:
            raise ValidationError("algebraic constraints require non-empty names")
        residual = _finite(record.residual, f"{record.name} residual")
        scale = _positive(record.scale, f"{record.name} scale")
        if abs(residual) / scale > limit:
            raise ValidationError(f"inconsistent algebraic initialization: {record.name}")
    return records


@dataclass(frozen=True, slots=True)
class StateEvent:
    time_s: float
    state_index: int
    value: float


@dataclass(frozen=True, slots=True)
class DynamicTrajectory:
    times_s: Vector
    states: tuple[Vector, ...]


def fixed_step_explicit_euler(
    rhs: ExplicitRhs,
    initial_state: Sequence[float],
    duration_s: float,
    step_s: float,
    events: Sequence[StateEvent] = (),
) -> DynamicTrajectory:
    """Integrate a reduced ODE reproducibly; events apply at recorded step boundaries."""
    if not callable(rhs):
        raise ValidationError("fixed-step ODE right-hand side must be callable")
    state = _finite_vector(initial_state, "initial ODE state")
    duration = _positive(duration_s, "ODE duration")
    step = _positive(step_s, "ODE step")
    step_count_float = duration / step
    step_count = round(step_count_float)
    if not math.isclose(step_count_float, step_count, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValidationError("fixed-step duration must be an integer number of steps")
    event_map: dict[int, list[StateEvent]] = {}
    for event in tuple(events):
        if not isinstance(event, StateEvent):
            raise ValidationError("fixed-step events must be StateEvent records")
        event_time = _finite(event.time_s, "event time")
        event_step_float = event_time / step
        event_step = round(event_step_float)
        if (
            event_time < 0.0
            or event_time > duration
            or not math.isclose(event_step_float, event_step, rel_tol=0.0, abs_tol=1.0e-12)
        ):
            raise ValidationError("fixed-step event times must align with the integration grid")
        if isinstance(event.state_index, bool) or not isinstance(event.state_index, int):
            raise ValidationError("event state index must be an integer")
        if not 0 <= event.state_index < len(state):
            raise ValidationError("event state index is outside the state vector")
        _finite(event.value, "event value")
        event_map.setdefault(event_step, []).append(event)

    times: list[float] = []
    states: list[Vector] = []
    for index in range(step_count + 1):
        if index in event_map:
            mutable = list(state)
            for event in event_map[index]:
                mutable[event.state_index] = float(event.value)
            state = tuple(mutable)
        time = index * step
        times.append(time)
        states.append(state)
        if index == step_count:
            break
        derivative = _finite_vector(rhs(time, state), "ODE derivative")
        if len(derivative) != len(state):
            raise ValidationError("ODE derivative length does not match the state")
        state = tuple(value + step * rate for value, rate in zip(state, derivative, strict=True))
        if not all(math.isfinite(value) for value in state):
            raise ValidationError("fixed-step ODE produced an unrepresentable state")
    return DynamicTrajectory(tuple(times), tuple(states))


def adaptive_explicit_ode(
    rhs: ExplicitRhs,
    initial_state: Sequence[float],
    duration_s: float,
    *,
    output_times_s: Sequence[float] | None = None,
    relative_tolerance: float = 1.0e-7,
    absolute_tolerance: float = 1.0e-9,
    algebraic_equations_present: bool = False,
) -> DynamicTrajectory:
    """Integrate an explicitly reduced ODE with SciPy's adaptive solver."""
    if algebraic_equations_present:
        raise ValidationError("DAE models require an IDA-capable solver; solve_ivp accepts only reduced ODEs")
    if not callable(rhs):
        raise ValidationError("adaptive ODE right-hand side must be callable")
    state = _finite_vector(initial_state, "initial ODE state")
    duration = _positive(duration_s, "ODE duration")
    rtol = _positive(relative_tolerance, "relative tolerance")
    atol = _positive(absolute_tolerance, "absolute tolerance")
    output_times = None
    if output_times_s is not None:
        output_times = _finite_vector(output_times_s, "ODE output times")
        if any(time < 0.0 or time > duration for time in output_times):
            raise ValidationError("ODE output times must lie in the integration interval")
        if any(right <= left for left, right in zip(output_times, output_times[1:])):
            raise ValidationError("ODE output times must be strictly increasing")

    def scipy_rhs(time: float, values: Sequence[float]) -> Vector:
        derivative = _finite_vector(rhs(time, tuple(float(value) for value in values)), "ODE derivative")
        if len(derivative) != len(state):
            raise ValidationError("ODE derivative length does not match the state")
        return derivative

    solution = solve_ivp(
        scipy_rhs, (0.0, duration), state, t_eval=output_times, rtol=rtol, atol=atol,
    )
    if not solution.success:
        raise ConvergenceError(f"adaptive ODE integration failed: {solution.message}")
    return DynamicTrajectory(
        tuple(float(value) for value in solution.t),
        tuple(tuple(float(value) for value in column) for column in solution.y.T),
    )


@dataclass(frozen=True, slots=True)
class PIDConfig:
    proportional_gain: float
    integral_gain: float
    derivative_gain: float
    setpoint: float
    output_scale: float
    integral_guard: float
    reverse_acting: bool = False
    output_min: float = -math.inf
    output_max: float = math.inf


@dataclass(frozen=True, slots=True)
class PIDState:
    integral: float = 0.0
    previous_error: float = 0.0


@dataclass(frozen=True, slots=True)
class PIDStepResult:
    state: PIDState
    normalized_error: float
    output: float


def pid_step(config: PIDConfig, state: PIDState, process_value: float, step_s: float) -> PIDStepResult:
    """Apply the discrete DWSIM PID equation used by the tank reference."""
    kp = _finite(config.proportional_gain, "PID proportional gain")
    ki = _finite(config.integral_gain, "PID integral gain")
    kd = _finite(config.derivative_gain, "PID derivative gain")
    setpoint = _finite(config.setpoint, "PID setpoint")
    if setpoint == 0.0:
        raise ValidationError("PID setpoint must be non-zero for normalized error")
    scale = _finite(config.output_scale, "PID output scale")
    guard = _positive(config.integral_guard, "PID integral guard")
    output_min = _extended_bound(config.output_min, "PID output minimum")
    output_max = _extended_bound(config.output_max, "PID output maximum")
    if output_min > output_max:
        raise ValidationError("PID output bounds are invalid")
    duration = _positive(step_s, "PID step")
    value = _finite(process_value, "PID process value")
    integral = _finite(state.integral, "PID integral")
    previous_error = _finite(state.previous_error, "PID previous error")
    error = (value - setpoint) / abs(setpoint)
    integral = min(max(integral + error * duration, -guard), guard)
    derivative = (error - previous_error) / duration if previous_error != 0.0 else 0.0
    relative_output = kp * error + ki * integral + kd * derivative
    if config.reverse_acting:
        output = (1.0 + relative_output) * scale
    else:
        output = (1.0 - relative_output) * scale
    output = min(max(output, output_min), output_max)
    return PIDStepResult(PIDState(integral, error), error, output)


@dataclass(frozen=True, slots=True)
class TankLevelControlConfig:
    duration_s: float
    step_s: float
    tank_volume_m3: float
    tank_height_m: float
    liquid_density_kg_m3: float
    inlet_mass_flow_kg_s: float
    initial_contents_volume_m3: float
    initial_outlet_opening_percent: float
    inlet_opening_percent: float
    tank_base_pressure_pa: float
    downstream_pressure_pa: float
    gravity_m_s2: float
    valve_kv: float
    pid: PIDConfig


@dataclass(frozen=True, slots=True)
class TankControlPoint:
    time_s: float
    liquid_level_m: float
    inlet_opening_percent: float
    outlet_opening_percent: float
    outlet_mass_flow_kg_s: float


@dataclass(frozen=True, slots=True)
class TankControlTrajectory:
    points: tuple[TankControlPoint, ...]
    mass_balance_residuals_kg: Vector


def simulate_dwsim_tank_level_control(config: TankLevelControlConfig) -> TankControlTrajectory:
    """Reproduce DWSIM's explicit tank/valve/PID schedule and one-step flow lag."""
    duration = _positive(config.duration_s, "tank simulation duration")
    step = _positive(config.step_s, "tank simulation step")
    step_count_float = duration / step
    step_count = round(step_count_float)
    if not math.isclose(step_count_float, step_count, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValidationError("tank duration must be an integer number of steps")
    volume = _positive(config.tank_volume_m3, "tank volume")
    height = _positive(config.tank_height_m, "tank height")
    area = volume / height
    density = _positive(config.liquid_density_kg_m3, "liquid density")
    inlet_mass_flow = _finite(config.inlet_mass_flow_kg_s, "inlet mass flow")
    if inlet_mass_flow < 0.0:
        raise ValidationError("inlet mass flow must be non-negative")
    initial_volume = _finite(config.initial_contents_volume_m3, "initial contents volume")
    if not 0.0 <= initial_volume <= volume:
        raise ValidationError("initial contents volume must lie within the vessel")
    opening = _finite(config.initial_outlet_opening_percent, "initial outlet opening")
    inlet_opening = _finite(config.inlet_opening_percent, "inlet opening")
    if not 0.0 <= opening <= 100.0 or not 0.0 <= inlet_opening <= 100.0:
        raise ValidationError("valve openings must be between zero and 100 percent")
    base_pressure = _positive(config.tank_base_pressure_pa, "tank base pressure")
    downstream_pressure = _positive(config.downstream_pressure_pa, "downstream pressure")
    gravity = _positive(config.gravity_m_s2, "gravity")
    valve_kv = _positive(config.valve_kv, "valve Kv")

    inlet_volume_flow = inlet_mass_flow / density
    level = initial_volume / area
    outlet_mass_flow = 0.0
    pid_state = PIDState()
    points: list[TankControlPoint] = []
    residuals: list[float] = []
    for index in range(step_count + 1):
        points.append(TankControlPoint(index * step, level, inlet_opening, opening, outlet_mass_flow))
        if index == step_count:
            break
        controller = pid_step(config.pid, pid_state, level, step)
        pid_state = controller.state
        next_opening = min(max(controller.output, 0.0), 100.0)

        old_mass = density * area * level
        next_level = level + step * (inlet_volume_flow - outlet_mass_flow / density) / area
        # DWSIM's dynamic tank reports overfill instead of clipping the state at
        # the nominal vessel height; negative inventory is still non-physical.
        if next_level < 0.0 or not math.isfinite(next_level):
            raise ValidationError("tank level became negative or unrepresentable")
        next_mass = density * area * next_level
        residuals.append(next_mass - old_mass - step * (inlet_mass_flow - outlet_mass_flow))

        upstream_pressure = base_pressure + density * gravity * next_level
        pressure_drop = upstream_pressure - downstream_pressure
        if pressure_drop <= 0.0 or next_opening <= 0.0:
            next_outlet_mass_flow = 0.0
        else:
            effective_kv = valve_kv * next_opening / 100.0
            next_outlet_mass_flow = effective_kv * math.sqrt(
                1000.0 * density * pressure_drop / 100000.0
            ) / 3600.0
        level = next_level
        opening = next_opening
        outlet_mass_flow = next_outlet_mass_flow
    return TankControlTrajectory(tuple(points), tuple(residuals))


@dataclass(frozen=True, slots=True)
class HeatExchangerState:
    hot_temperature_k: float
    cold_temperature_k: float


@dataclass(frozen=True, slots=True)
class HeatExchangerStepResult:
    state: HeatExchangerState
    transferred_heat_w: float
    total_energy_balance_residual_j: float


def lumped_heat_exchanger_step(
    state: HeatExchangerState,
    *,
    hot_inlet_temperature_k: float,
    cold_inlet_temperature_k: float,
    hot_mass_flow_kg_s: float,
    cold_mass_flow_kg_s: float,
    hot_holdup_kg: float,
    cold_holdup_kg: float,
    hot_heat_capacity_j_kg_k: float,
    cold_heat_capacity_j_kg_k: float,
    ua_w_k: float,
    step_s: float,
) -> HeatExchangerStepResult:
    """Advance two well-mixed sides with equal-and-opposite UA heat transfer."""
    hot_temperature = _positive(state.hot_temperature_k, "hot temperature")
    cold_temperature = _positive(state.cold_temperature_k, "cold temperature")
    hot_inlet = _positive(hot_inlet_temperature_k, "hot inlet temperature")
    cold_inlet = _positive(cold_inlet_temperature_k, "cold inlet temperature")
    hot_flow = _finite(hot_mass_flow_kg_s, "hot mass flow")
    cold_flow = _finite(cold_mass_flow_kg_s, "cold mass flow")
    if hot_flow < 0.0 or cold_flow < 0.0:
        raise ValidationError("heat-exchanger mass flows must be non-negative")
    hot_holdup = _positive(hot_holdup_kg, "hot holdup")
    cold_holdup = _positive(cold_holdup_kg, "cold holdup")
    hot_cp = _positive(hot_heat_capacity_j_kg_k, "hot heat capacity")
    cold_cp = _positive(cold_heat_capacity_j_kg_k, "cold heat capacity")
    ua = _finite(ua_w_k, "heat-exchanger UA")
    if ua < 0.0:
        raise ValidationError("heat-exchanger UA must be non-negative")
    duration = _positive(step_s, "heat-exchanger step")
    transferred_heat = ua * (hot_temperature - cold_temperature)
    hot_external_rate = hot_flow * hot_cp * (hot_inlet - hot_temperature)
    cold_external_rate = cold_flow * cold_cp * (cold_inlet - cold_temperature)
    hot_energy_change = duration * (hot_external_rate - transferred_heat)
    cold_energy_change = duration * (cold_external_rate + transferred_heat)
    next_hot = hot_temperature + hot_energy_change / (hot_holdup * hot_cp)
    next_cold = cold_temperature + cold_energy_change / (cold_holdup * cold_cp)
    if next_hot <= 0.0 or next_cold <= 0.0 or not math.isfinite(next_hot + next_cold):
        raise ValidationError("heat-exchanger step produced an invalid temperature")
    residual = hot_energy_change + cold_energy_change - duration * (
        hot_external_rate + cold_external_rate
    )
    return HeatExchangerStepResult(
        HeatExchangerState(next_hot, next_cold), transferred_heat, residual,
    )


def dynamic_cstr_step(
    state: HoldupState,
    *,
    component_in_kmol_s: Sequence[float],
    component_out_kmol_s: Sequence[float],
    stoichiometry: Sequence[Sequence[float]],
    reaction_extents_kmol_s: Sequence[float],
    energy_in_w: float,
    energy_out_w: float,
    reaction_enthalpies_j_kmol: Sequence[float],
    heat_input_w: float,
    step_s: float,
) -> HoldupStepResult:
    """Advance a dynamic CSTR from explicit reaction extents and stoichiometry."""
    extents = _finite_vector(reaction_extents_kmol_s, "reaction extents")
    enthalpies = _finite_vector(reaction_enthalpies_j_kmol, "reaction enthalpies")
    try:
        coefficients = tuple(tuple(_finite(value, "stoichiometric coefficient") for value in row) for row in stoichiometry)
    except TypeError as error:
        raise ValidationError("CSTR stoichiometry must be a finite matrix") from error
    component_count = len(state.component_amounts_kmol)
    if len(coefficients) != component_count or any(len(row) != len(extents) for row in coefficients):
        raise ValidationError("CSTR stoichiometry dimensions do not match components and reactions")
    if len(enthalpies) != len(extents):
        raise ValidationError("CSTR reaction enthalpies must match reaction extents")
    generation = tuple(
        math.fsum(coefficient * extent for coefficient, extent in zip(row, extents, strict=True))
        for row in coefficients
    )
    reaction_energy = -math.fsum(
        extent * enthalpy for extent, enthalpy in zip(extents, enthalpies, strict=True)
    )
    try:
        inlet = tuple(component_in_kmol_s)
        outlet = tuple(component_out_kmol_s)
    except TypeError as error:
        raise ValidationError("CSTR component rates must be finite sequences") from error
    return advance_holdup(
        state,
        HoldupRates(
            inlet, outlet, generation,
            energy_in_w, energy_out_w, heat_input_w, reaction_energy,
        ),
        step_s,
    )
