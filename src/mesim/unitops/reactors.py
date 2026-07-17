"""Reaction unit operations."""
import math
import sys
from dataclasses import dataclass

from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from ..compounds import Compound, PRInteractions
from ..errors import ConvergenceError, ValidationError
from ..reactions import CompoundThermochemistry, ReactionDefinition
from ..thermo.ideal import IdealCorrelations
from ..thermo.peng_robinson import PengRobinsonMixture


DWSIM_EQUILIBRIUM_R = 8314.0  # J/kmol/K, as used by DWSIM AUX_DELGF_T
DWSIM_REFERENCE_PRESSURE_PA = 101325.0
DWSIM_KINETIC_R = 8.31446261815324  # J/mol/K


@dataclass(frozen=True, slots=True)
class ConversionReactorResult:
    outlet_component_flows_kmol_s: tuple[tuple[str, float], ...]
    extent_kmol_s: float
    conversion_fraction: float
    reaction_heat_w: float

    @property
    def total_molar_flow_kmol_s(self) -> float:
        return math.fsum(flow for _, flow in self.outlet_component_flows_kmol_s)


@dataclass(frozen=True, slots=True)
class EquilibriumReactorResult:
    outlet_component_flows_kmol_s: tuple[tuple[str, float], ...]
    extents_kmol_s: tuple[tuple[str, float], ...]
    component_conversions: tuple[tuple[str, float], ...]
    equilibrium_log_residuals: tuple[tuple[str, float], ...]
    reference_reaction_heat_w: float
    isothermal_duty_w: float
    iterations: int

    @property
    def total_molar_flow_kmol_s(self) -> float:
        return math.fsum(flow for _, flow in self.outlet_component_flows_kmol_s)


@dataclass(frozen=True, slots=True)
class GibbsReactorResult:
    outlet_component_flows_kmol_s: tuple[tuple[str, float], ...]
    component_conversions: tuple[tuple[str, float], ...]
    element_balance_residuals_kmol_s: tuple[tuple[str, float], ...]
    stationarity_residual: float
    initial_gibbs_energy_w: float
    final_gibbs_energy_w: float
    isothermal_duty_w: float
    iterations: int

    @property
    def total_molar_flow_kmol_s(self) -> float:
        return math.fsum(flow for _, flow in self.outlet_component_flows_kmol_s)


@dataclass(frozen=True, slots=True)
class CSTRResult:
    outlet_component_flows_kmol_s: tuple[tuple[str, float], ...]
    extent_kmol_s: float
    reaction_rate_kmol_m3_s: float
    component_conversions: tuple[tuple[str, float], ...]
    reference_reaction_heat_w: float
    material_rate_residual_kmol_s: float
    evaluations: int

    @property
    def total_molar_flow_kmol_s(self) -> float:
        return math.fsum(flow for _, flow in self.outlet_component_flows_kmol_s)


@dataclass(frozen=True, slots=True)
class PFRResult:
    outlet_component_flows_kmol_s: tuple[tuple[str, float], ...]
    extent_kmol_s: float
    average_reaction_rate_kmol_m3_s: float
    component_conversions: tuple[tuple[str, float], ...]
    reference_reaction_heat_w: float
    material_balance_residual_kmol_s: float
    integration_steps: int

    @property
    def total_molar_flow_kmol_s(self) -> float:
        return math.fsum(flow for _, flow in self.outlet_component_flows_kmol_s)


def conversion_reactor(
    inlet_component_flows_kmol_s: tuple[tuple[str, float], ...],
    reaction: ReactionDefinition,
    conversion_fraction: float | None = None,
) -> ConversionReactorResult:
    """Apply one balanced DWSIM-style conversion reaction on a kmol/s basis."""
    try:
        inlet_items = tuple(inlet_component_flows_kmol_s)
    except TypeError as error:
        raise ValidationError("conversion reactor inlet flows must be a finite sequence") from error
    if not inlet_items or len({compound for compound, _ in inlet_items}) != len(inlet_items):
        raise ValidationError("conversion reactor inlet compound IDs must be non-empty and unique")
    inlet = {}
    for compound, flow in inlet_items:
        if not isinstance(compound, str) or not compound:
            raise ValidationError("conversion reactor compound IDs must be non-empty strings")
        if isinstance(flow, bool) or not isinstance(flow, (int, float)) or not math.isfinite(flow) or flow < 0.0:
            raise ValidationError("conversion reactor component flows must be finite and non-negative")
        inlet[compound] = float(flow)

    if reaction.reaction_type != "conversion" or reaction.conversion_fraction is None:
        raise ValidationError("conversion reactor requires a conversion reaction")
    conversion = reaction.conversion_fraction if conversion_fraction is None else conversion_fraction
    if isinstance(conversion, bool) or not isinstance(conversion, (int, float)) or not math.isfinite(conversion) or not 0.0 <= conversion <= 1.0:
        raise ValidationError("conversion reactor conversion must be between zero and one")
    stoichiometry = dict(reaction.stoichiometry)
    base_coefficient = stoichiometry.get(reaction.base_reactant)
    base_flow = inlet.get(reaction.base_reactant, 0.0)
    if base_coefficient is None or base_coefficient >= 0.0 or base_flow <= 0.0:
        raise ValidationError("conversion reactor requires a flowing base reactant with a negative coefficient")

    extent = conversion * base_flow / -base_coefficient
    ordered_compounds = tuple(inlet) + tuple(compound for compound in stoichiometry if compound not in inlet)
    outlet = []
    for compound in ordered_compounds:
        flow = inlet.get(compound, 0.0) + stoichiometry.get(compound, 0.0) * extent
        if flow < -1.0e-12:
            raise ValidationError(f"conversion reactor produced a negative {compound} flow")
        outlet.append((compound, max(flow, 0.0)))
    return ConversionReactorResult(
        tuple(outlet), extent, float(conversion), extent * reaction.reaction_heat_j_per_kmol,
    )


def _arrhenius_rate(
    pre_exponential_factor: float,
    activation_energy_j_per_mol: float,
    orders: tuple[tuple[str, float], ...],
    concentrations_kmol_m3: dict[str, float],
    temperature_k: float,
) -> float:
    if pre_exponential_factor == 0.0:
        return 0.0
    logarithm = math.log(pre_exponential_factor) - activation_energy_j_per_mol / (
        DWSIM_KINETIC_R * temperature_k
    )
    for compound, order in orders:
        if order == 0.0:
            continue
        concentration = concentrations_kmol_m3.get(compound, 0.0)
        if concentration <= 0.0:
            return 0.0
        logarithm += order * math.log(concentration)
    if logarithm > math.log(sys.float_info.max):
        raise ValidationError("kinetic rate is outside the representable range")
    if logarithm < math.log(sys.float_info.min):
        return 0.0
    return math.exp(logarithm)


def _kinetic_rate_kmol_m3_s(
    reaction: ReactionDefinition,
    concentrations_kmol_m3: dict[str, float],
    temperature_k: float,
) -> float:
    kinetics = reaction.kinetics
    if kinetics is None:
        raise ValidationError("kinetic reaction is missing its rate definition")
    forward = _arrhenius_rate(
        kinetics.forward.pre_exponential_factor,
        kinetics.forward.activation_energy_j_per_mol,
        kinetics.forward.orders,
        concentrations_kmol_m3,
        temperature_k,
    )
    reverse = _arrhenius_rate(
        kinetics.reverse.pre_exponential_factor,
        kinetics.reverse.activation_energy_j_per_mol,
        kinetics.reverse.orders,
        concentrations_kmol_m3,
        temperature_k,
    )
    rate = forward - reverse
    if kinetics.rate_unit == "kmol/[m3.h]":
        rate /= 3600.0
    return rate


def continuous_stirred_tank_reactor(
    inlet_component_flows_kmol_s: tuple[tuple[str, float], ...],
    reaction: ReactionDefinition,
    temperature_k: float,
    reactor_volume_m3: float,
    outlet_volumetric_flow_m3_s: float,
    tolerance: float = 1.0e-10,
    max_evaluations: int = 100,
) -> CSTRResult:
    """Solve one steady, ideally mixed kinetic reaction using its original DWSIM units."""
    inlet_items = tuple(inlet_component_flows_kmol_s)
    if not inlet_items or len({compound for compound, _ in inlet_items}) != len(inlet_items):
        raise ValidationError("CSTR inlet compound IDs must be non-empty and unique")
    inlet: dict[str, float] = {}
    for compound, flow in inlet_items:
        if not isinstance(compound, str) or not compound:
            raise ValidationError("CSTR compound IDs must be non-empty strings")
        if isinstance(flow, bool) or not isinstance(flow, (int, float)) or not math.isfinite(flow) or flow < 0.0:
            raise ValidationError("CSTR component flows must be finite and non-negative")
        inlet[compound] = float(flow)
    if reaction.reaction_type != "kinetic" or reaction.kinetics is None:
        raise ValidationError("CSTR requires a kinetic reaction")
    if reaction.phase not in {"mixture", "liquid"}:
        raise ValidationError("CSTR currently supports mixture or liquid kinetic reactions")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value) or value <= 0.0
        for value in (temperature_k, reactor_volume_m3, outlet_volumetric_flow_m3_s, tolerance)
    ) or isinstance(max_evaluations, bool) or not isinstance(max_evaluations, int) or max_evaluations <= 0:
        raise ValidationError("CSTR state and solver controls must be finite and positive")

    stoichiometry = dict(reaction.stoichiometry)
    ordered_compounds = tuple(inlet) + tuple(compound for compound in stoichiometry if compound not in inlet)
    limits = [
        inlet.get(compound, 0.0) / -coefficient
        for compound, coefficient in stoichiometry.items()
        if coefficient < 0.0
    ]
    if not limits or min(limits) <= 0.0:
        raise ValidationError("CSTR reaction requires every consumed compound to have positive inlet flow")
    limiting_extent = min(limits)

    def outlet_flows(extent: float) -> tuple[float, ...]:
        return tuple(
            inlet.get(compound, 0.0) + stoichiometry.get(compound, 0.0) * extent
            for compound in ordered_compounds
        )

    def rate_at(extent: float) -> float:
        flows = outlet_flows(extent)
        if any(flow < 0.0 for flow in flows):
            raise ValidationError("CSTR trial produced a negative component flow")
        concentrations = {
            compound: flow / outlet_volumetric_flow_m3_s
            for compound, flow in zip(ordered_compounds, flows)
        }
        return _kinetic_rate_kmol_m3_s(reaction, concentrations, float(temperature_k))

    initial_rate = rate_at(0.0)
    if initial_rate < 0.0:
        raise ValidationError("CSTR reverse-rate-dominated solutions are outside the current domain")
    scale = max(limiting_extent, 1.0e-12)
    initial_extent = min(max(initial_rate * reactor_volume_m3, limiting_extent * 1.0e-12), limiting_extent * 0.5)

    def residual(values) -> tuple[float]:
        extent = float(values[0])
        return ((extent - rate_at(extent) * reactor_volume_m3) / scale,)

    solved = least_squares(
        residual, (initial_extent,), bounds=(0.0, limiting_extent),
        ftol=1.0e-14, xtol=1.0e-14, gtol=1.0e-14, max_nfev=max_evaluations,
    )
    extent = float(solved.x[0])
    reaction_rate = rate_at(extent)
    material_residual = extent - reaction_rate * reactor_volume_m3
    if not solved.success or abs(material_residual) > tolerance * scale:
        raise ConvergenceError("CSTR material-rate balance did not converge")

    outlet_values = outlet_flows(extent)
    outlet = tuple(zip(ordered_compounds, outlet_values))
    conversions = tuple(
        (compound, (flow - outlet_values[index]) / flow)
        for index, (compound, flow) in enumerate(inlet.items())
        if flow > 0.0 and outlet_values[index] < flow
    )
    return CSTRResult(
        outlet,
        extent,
        reaction_rate,
        conversions,
        extent * reaction.reaction_heat_j_per_kmol,
        material_residual,
        solved.nfev,
    )


def plug_flow_reactor(
    inlet_component_flows_kmol_s: tuple[tuple[str, float], ...],
    reaction: ReactionDefinition,
    inlet_temperature_k: float,
    outlet_temperature_k: float,
    reactor_volume_m3: float,
    inlet_volumetric_flow_m3_s: float,
    outlet_volumetric_flow_m3_s: float,
    relative_tolerance: float = 1.0e-10,
    absolute_tolerance: float = 1.0e-14,
    max_step_m3: float | None = None,
) -> PFRResult:
    """Integrate one kinetic reaction over reactor volume with supplied T/Q endpoints."""
    inlet_items = tuple(inlet_component_flows_kmol_s)
    if not inlet_items or len({compound for compound, _ in inlet_items}) != len(inlet_items):
        raise ValidationError("PFR inlet compound IDs must be non-empty and unique")
    inlet: dict[str, float] = {}
    for compound, flow in inlet_items:
        if not isinstance(compound, str) or not compound:
            raise ValidationError("PFR compound IDs must be non-empty strings")
        if isinstance(flow, bool) or not isinstance(flow, (int, float)) or not math.isfinite(flow) or flow < 0.0:
            raise ValidationError("PFR component flows must be finite and non-negative")
        inlet[compound] = float(flow)
    if reaction.reaction_type != "kinetic" or reaction.kinetics is None:
        raise ValidationError("PFR requires a kinetic reaction")
    if reaction.phase not in {"mixture", "liquid"}:
        raise ValidationError("PFR currently supports mixture or liquid kinetic reactions")
    scalars = (
        inlet_temperature_k, outlet_temperature_k, reactor_volume_m3,
        inlet_volumetric_flow_m3_s, outlet_volumetric_flow_m3_s,
        relative_tolerance, absolute_tolerance,
    )
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value) or value <= 0.0
        for value in scalars
    ):
        raise ValidationError("PFR state and integration tolerances must be finite and positive")
    if max_step_m3 is not None and (
        isinstance(max_step_m3, bool) or not isinstance(max_step_m3, (int, float))
        or not math.isfinite(max_step_m3) or max_step_m3 <= 0.0
    ):
        raise ValidationError("PFR maximum step must be finite and positive")

    stoichiometry = dict(reaction.stoichiometry)
    if any(inlet.get(compound, 0.0) <= 0.0 for compound, coefficient in stoichiometry.items() if coefficient < 0.0):
        raise ValidationError("PFR reaction requires every consumed compound to have positive inlet flow")
    ordered_compounds = tuple(inlet) + tuple(compound for compound in stoichiometry if compound not in inlet)
    initial = tuple(inlet.get(compound, 0.0) for compound in ordered_compounds)

    def derivative(volume: float, flows) -> tuple[float, ...]:
        fraction = volume / reactor_volume_m3
        temperature = inlet_temperature_k + fraction * (outlet_temperature_k - inlet_temperature_k)
        volumetric_flow = inlet_volumetric_flow_m3_s + fraction * (
            outlet_volumetric_flow_m3_s - inlet_volumetric_flow_m3_s
        )
        concentrations = {
            compound: max(float(flow), 0.0) / volumetric_flow
            for compound, flow in zip(ordered_compounds, flows)
        }
        rate = _kinetic_rate_kmol_m3_s(reaction, concentrations, temperature)
        if rate < 0.0:
            raise ValidationError("PFR reverse-rate-dominated solutions are outside the current domain")
        return tuple(stoichiometry.get(compound, 0.0) * rate for compound in ordered_compounds)

    integration = solve_ivp(
        derivative,
        (0.0, reactor_volume_m3),
        initial,
        method="RK45",
        rtol=relative_tolerance,
        atol=absolute_tolerance,
        max_step=reactor_volume_m3 if max_step_m3 is None else max_step_m3,
    )
    if not integration.success:
        raise ConvergenceError(f"PFR integration failed: {integration.message}")
    outlet_values = tuple(float(value) for value in integration.y[:, -1])
    if any(value < -absolute_tolerance for value in outlet_values):
        raise ConvergenceError("PFR integration produced a negative component flow")
    outlet_values = tuple(max(value, 0.0) for value in outlet_values)
    base_coefficient = stoichiometry[reaction.base_reactant]
    base_index = ordered_compounds.index(reaction.base_reactant)
    extent = (outlet_values[base_index] - initial[base_index]) / base_coefficient
    residual = max(
        abs(
            outlet_values[index]
            - (initial[index] + stoichiometry.get(compound, 0.0) * extent)
        )
        for index, compound in enumerate(ordered_compounds)
    )
    conversions = tuple(
        (compound, (flow - outlet_values[index]) / flow)
        for index, (compound, flow) in enumerate(inlet.items())
        if flow > 0.0 and outlet_values[index] < flow
    )
    return PFRResult(
        tuple(zip(ordered_compounds, outlet_values)),
        extent,
        extent / reactor_volume_m3,
        conversions,
        extent * reaction.reaction_heat_j_per_kmol,
        residual,
        len(integration.t) - 1,
    )


def _dwsim_midpoint_integrals(
    correlation: IdealCorrelations, start_k: float, end_k: float,
) -> tuple[float, float]:
    """Reproduce DWSIM's AUX_INT_CPDTi/AUX_INT_CPDT_Ti midpoint quadrature."""
    difference = abs(end_k - start_k)
    base_steps = int(round(difference / 10.0))
    if base_steps < 10:
        if difference < 1.0:
            base_steps = 2
        elif difference < 3.0:
            base_steps = 4
        elif difference < 5.0:
            base_steps = 6
        else:
            base_steps = 10

    def integrate(steps: int, divide_by_temperature: bool) -> float:
        delta = (end_k - start_k) / steps
        values = tuple(start_k + delta / 2.0 + index * delta for index in range(steps))
        heat_capacities = tuple(
            correlation.heat_capacity(value, allow_extrapolation=True).value for value in values
        )
        if divide_by_temperature:
            return math.fsum(
                heat_capacity / temperature
                for heat_capacity, temperature in zip(heat_capacities, values)
            ) * delta
        return math.fsum(heat_capacities) * delta

    return integrate(min(base_steps, 100), False), integrate(base_steps, True)


def _standard_gibbs_rt(
    thermochemistry: CompoundThermochemistry,
    correlation: IdealCorrelations,
    temperature_k: float,
) -> float:
    enthalpy_change, entropy_change = _dwsim_midpoint_integrals(
        correlation, thermochemistry.formation_temperature_k, temperature_k,
    )
    formation_entropy = (
        thermochemistry.ideal_gas_formation_enthalpy_j_per_kmol
        - thermochemistry.ideal_gas_formation_gibbs_energy_j_per_kmol
    ) / thermochemistry.formation_temperature_k
    standard_gibbs = (
        thermochemistry.ideal_gas_formation_enthalpy_j_per_kmol
        + enthalpy_change
        - temperature_k * (formation_entropy + entropy_change)
    )
    return standard_gibbs / (DWSIM_EQUILIBRIUM_R * temperature_k)


def _solve_linear(matrix: list[list[float]], right: list[float]) -> list[float]:
    size = len(right)
    augmented = [row[:] + [value] for row, value in zip(matrix, right)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1.0e-14:
            raise ConvergenceError("equilibrium reactor Jacobian is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][-1] for row in range(size)]


def equilibrium_reactor(
    inlet_component_flows_kmol_s: tuple[tuple[str, float], ...],
    reactions: tuple[ReactionDefinition, ...],
    thermochemistry: tuple[CompoundThermochemistry, ...],
    compounds: tuple[Compound, ...],
    correlations: tuple[IdealCorrelations, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    initial_extents_kmol_s: tuple[float, ...] | None = None,
    tolerance: float = 1.0e-10,
    max_iterations: int = 60,
) -> EquilibriumReactorResult:
    """Solve simultaneous vapor reactions on the same kmol/s and PR fugacity basis as DWSIM."""
    inlet_items = tuple(inlet_component_flows_kmol_s)
    if not inlet_items or len({compound for compound, _ in inlet_items}) != len(inlet_items):
        raise ValidationError("equilibrium reactor inlet compound IDs must be non-empty and unique")
    if not reactions or any(
        reaction.reaction_type != "equilibrium"
        or reaction.equilibrium_constant_model != "gibbs"
        or reaction.reaction_basis != "fugacity"
        or reaction.phase != "vapor"
        for reaction in reactions
    ):
        raise ValidationError("equilibrium reactor requires vapor Gibbs-fugacity reactions")
    scalar_values = (temperature_k, pressure_pa, tolerance)
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value) or value <= 0.0
        for value in scalar_values
    ) or isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValidationError("equilibrium reactor state and solver controls must be finite and positive")

    inlet = {}
    for compound_id, flow in inlet_items:
        if not isinstance(compound_id, str) or not compound_id:
            raise ValidationError("equilibrium reactor compound IDs must be non-empty strings")
        if isinstance(flow, bool) or not isinstance(flow, (int, float)) or not math.isfinite(flow) or flow < 0.0:
            raise ValidationError("equilibrium reactor component flows must be finite and non-negative")
        inlet[compound_id] = float(flow)
    compound_ids = tuple(inlet)
    if tuple(compound.id for compound in compounds) != compound_ids:
        raise ValidationError("equilibrium reactor compound order must match the inlet flow order")
    thermo_by_id = {record.compound_id: record for record in thermochemistry}
    correlation_by_id = {record.compound_id: record for record in correlations}
    if any(compound_id not in thermo_by_id or compound_id not in correlation_by_id for compound_id in compound_ids):
        raise ValidationError("equilibrium reactor is missing thermochemistry or heat-capacity data")
    if any(compound not in inlet for reaction in reactions for compound, _ in reaction.stoichiometry):
        raise ValidationError("equilibrium reactor inlet must explicitly include every reaction compound")

    stoichiometry = tuple(
        tuple(dict(reaction.stoichiometry).get(compound_id, 0.0) for compound_id in compound_ids)
        for reaction in reactions
    )
    standard_gibbs_rt = tuple(
        _standard_gibbs_rt(thermo_by_id[compound_id], correlation_by_id[compound_id], temperature_k)
        for compound_id in compound_ids
    )
    inlet_vector = tuple(inlet.values())
    total_inlet = math.fsum(inlet_vector)
    positive_floor = max(total_inlet * 1.0e-14, 1.0e-18)

    def flows(extents: tuple[float, ...] | list[float]) -> tuple[float, ...]:
        return tuple(
            inlet_vector[index] + math.fsum(
                stoichiometry[reaction_index][index] * extents[reaction_index]
                for reaction_index in range(len(reactions))
            )
            for index in range(len(compound_ids))
        )

    if initial_extents_kmol_s is None:
        guesses = []
        for reaction, coefficients in zip(reactions, stoichiometry):
            limits = [
                inlet_vector[index] / -coefficient
                for index, coefficient in enumerate(coefficients)
                if coefficient < 0.0 and inlet_vector[index] > 0.0
            ]
            if not limits:
                raise ValidationError(f"equilibrium reaction {reaction.id} has no flowing reactant")
            guesses.append(0.1 * min(limits) / len(reactions))
        extents = guesses
    else:
        if len(initial_extents_kmol_s) != len(reactions) or any(
            isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            for value in initial_extents_kmol_s
        ):
            raise ValidationError("equilibrium reactor initial extents must match the reaction count")
        extents = [float(value) for value in initial_extents_kmol_s]

    active_compounds = {
        compound_id for reaction in reactions for compound_id, _ in reaction.stoichiometry
    }

    def residual(extent_values: tuple[float, ...] | list[float]) -> tuple[float, ...]:
        component_flows = flows(extent_values)
        if any(
            flow <= positive_floor for compound_id, flow in zip(compound_ids, component_flows)
            if compound_id in active_compounds
        ) or any(flow < 0.0 for flow in component_flows):
            raise ValidationError("equilibrium trial produced a non-positive reacting component flow")
        total = math.fsum(component_flows)
        composition = tuple(flow / total for flow in component_flows)
        state = PengRobinsonMixture(compounds, composition, interactions).state(
            temperature_k, pressure_pa, "vapor",
        )
        log_activities = tuple(
            math.log(fraction) + math.log(coefficient) + math.log(pressure_pa / DWSIM_REFERENCE_PRESSURE_PA)
            for fraction, coefficient in zip(composition, state.fugacity_coefficients)
        )
        return tuple(
            math.fsum(
                coefficient * (standard_gibbs_rt[index] + log_activities[index])
                for index, coefficient in enumerate(reaction_coefficients)
            )
            for reaction_coefficients in stoichiometry
        )

    last_residual: tuple[float, ...] | None = None
    for iteration in range(1, max_iterations + 1):
        current = residual(extents)
        last_residual = current
        norm = max(abs(value) for value in current)
        if norm <= tolerance:
            break
        jacobian = [[0.0 for _ in reactions] for _ in reactions]
        for column in range(len(reactions)):
            step = max(total_inlet * 1.0e-7, abs(extents[column]) * 1.0e-5, 1.0e-12)
            trial = extents[:]
            trial[column] += step
            try:
                shifted = residual(trial)
            except ValidationError:
                trial[column] -= 2.0 * step
                shifted = residual(trial)
                step = -step
            for row in range(len(reactions)):
                jacobian[row][column] = (shifted[row] - current[row]) / step
        direction = _solve_linear(jacobian, [-value for value in current])
        scale = 1.0
        accepted = False
        while scale >= 2.0**-24:
            trial = [value + scale * delta for value, delta in zip(extents, direction)]
            try:
                shifted = residual(trial)
            except ValidationError:
                scale *= 0.5
                continue
            if max(abs(value) for value in shifted) < norm:
                extents = trial
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            raise ConvergenceError("equilibrium reactor line search did not reduce the residual")
    else:
        raise ConvergenceError("equilibrium reactor did not converge")

    outlet_vector = flows(extents)
    outlet = tuple(zip(compound_ids, outlet_vector))
    conversions = tuple(
        (compound_id, (inlet[compound_id] - outlet_vector[index]) / inlet[compound_id])
        for index, compound_id in enumerate(compound_ids)
        if inlet[compound_id] > 0.0 and outlet_vector[index] < inlet[compound_id]
    )
    reference_heat = math.fsum(
        extent * reaction.reaction_heat_j_per_kmol for extent, reaction in zip(extents, reactions)
    )

    def enthalpy_flow(component_flows: tuple[float, ...]) -> float:
        total = math.fsum(component_flows)
        composition = tuple(flow / total for flow in component_flows)
        state = PengRobinsonMixture(compounds, composition, interactions).state(
            temperature_k, pressure_pa, "vapor",
        )
        component_enthalpy = math.fsum(
            flow * (
                thermo_by_id[compound_id].ideal_gas_formation_enthalpy_j_per_kmol
                + _dwsim_midpoint_integrals(
                    correlation_by_id[compound_id],
                    thermo_by_id[compound_id].formation_temperature_k,
                    temperature_k,
                )[0]
            )
            for compound_id, flow in zip(compound_ids, component_flows)
        )
        return component_enthalpy + total * state.departure_enthalpy_j_per_kmol

    duty = enthalpy_flow(outlet_vector) - enthalpy_flow(inlet_vector)
    return EquilibriumReactorResult(
        outlet,
        tuple((reaction.id, extent) for reaction, extent in zip(reactions, extents)),
        conversions,
        tuple((reaction.id, value) for reaction, value in zip(reactions, last_residual or ())),
        reference_heat,
        duty,
        iteration,
    )


def gibbs_reactor(
    inlet_component_flows_kmol_s: tuple[tuple[str, float], ...],
    thermochemistry: tuple[CompoundThermochemistry, ...],
    compounds: tuple[Compound, ...],
    correlations: tuple[IdealCorrelations, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    tolerance: float = 1.0e-10,
    max_iterations: int = 80,
) -> GibbsReactorResult:
    """Minimize vapor Gibbs energy subject to exact elemental balances."""
    inlet_items = tuple(inlet_component_flows_kmol_s)
    if not inlet_items or len({compound for compound, _ in inlet_items}) != len(inlet_items):
        raise ValidationError("Gibbs reactor inlet compound IDs must be non-empty and unique")
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value) or value <= 0.0
        for value in (temperature_k, pressure_pa, tolerance)
    ) or isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValidationError("Gibbs reactor state and solver controls must be finite and positive")

    inlet = {}
    for compound_id, flow in inlet_items:
        if not isinstance(compound_id, str) or not compound_id:
            raise ValidationError("Gibbs reactor compound IDs must be non-empty strings")
        if isinstance(flow, bool) or not isinstance(flow, (int, float)) or not math.isfinite(flow) or flow < 0.0:
            raise ValidationError("Gibbs reactor component flows must be finite and non-negative")
        inlet[compound_id] = float(flow)
    compound_ids = tuple(inlet)
    if tuple(compound.id for compound in compounds) != compound_ids:
        raise ValidationError("Gibbs reactor compound order must match the inlet flow order")
    thermo_by_id = {record.compound_id: record for record in thermochemistry}
    correlation_by_id = {record.compound_id: record for record in correlations}
    if any(compound_id not in thermo_by_id or compound_id not in correlation_by_id for compound_id in compound_ids):
        raise ValidationError("Gibbs reactor is missing thermochemistry or heat-capacity data")
    total_inlet = math.fsum(inlet.values())
    if total_inlet <= 0.0:
        raise ValidationError("Gibbs reactor requires positive inlet molar flow")

    elements = tuple(sorted({element for compound_id in compound_ids for element, _ in thermo_by_id[compound_id].elements}))
    element_matrix = tuple(
        tuple(dict(thermo_by_id[compound_id].elements).get(element, 0.0) for compound_id in compound_ids)
        for element in elements
    )
    inlet_vector = tuple(inlet.values())
    element_totals = tuple(
        math.fsum(coefficient * flow for coefficient, flow in zip(row, inlet_vector))
        for row in element_matrix
    )
    if any(total <= 0.0 for total in element_totals):
        raise ValidationError("Gibbs reactor elements must have positive inlet totals")
    standard_gibbs_rt = tuple(
        _standard_gibbs_rt(thermo_by_id[compound_id], correlation_by_id[compound_id], temperature_k)
        for compound_id in compound_ids
    )

    def chemical_potentials_rt(component_flows: tuple[float, ...]) -> tuple[float, ...]:
        total = math.fsum(component_flows)
        composition = tuple(flow / total for flow in component_flows)
        state = PengRobinsonMixture(compounds, composition, interactions).state(
            temperature_k, pressure_pa, "vapor",
        )
        return tuple(
            standard + math.log(fraction) + math.log(fugacity)
            + math.log(pressure_pa / DWSIM_REFERENCE_PRESSURE_PA)
            for standard, fraction, fugacity in zip(
                standard_gibbs_rt, composition, state.fugacity_coefficients,
            )
        )

    floor = max(total_inlet * 1.0e-10, 1.0e-18)
    variables = [math.log(max(flow, floor)) for flow in inlet_vector] + [0.0] * len(elements)

    def evaluate(values: list[float]) -> tuple[tuple[float, ...], tuple[float, ...]]:
        try:
            component_flows = tuple(math.exp(value) for value in values[:len(compound_ids)])
        except OverflowError as error:
            raise ValidationError("Gibbs reactor trial flows are outside the representable range") from error
        if not all(math.isfinite(flow) and flow > 0.0 for flow in component_flows):
            raise ValidationError("Gibbs reactor trial flows must be finite and positive")
        multipliers = values[len(compound_ids):]
        chemical_potentials = chemical_potentials_rt(component_flows)
        stationarity = tuple(
            chemical_potentials[index] - math.fsum(
                element_matrix[element_index][index] * multipliers[element_index]
                for element_index in range(len(elements))
            )
            for index in range(len(compound_ids))
        )
        balances = tuple(
            (math.fsum(coefficient * flow for coefficient, flow in zip(row, component_flows)) - total) / total
            for row, total in zip(element_matrix, element_totals)
        )
        return stationarity + balances, component_flows

    final_residual: tuple[float, ...] | None = None
    outlet_vector: tuple[float, ...] | None = None
    for iteration in range(1, max_iterations + 1):
        current, current_flows = evaluate(variables)
        norm = max(abs(value) for value in current)
        final_residual, outlet_vector = current, current_flows
        if norm <= tolerance:
            break
        size = len(variables)
        jacobian = [[0.0] * size for _ in range(size)]
        for column in range(size):
            step = 1.0e-6 * max(1.0, abs(variables[column]))
            trial = variables[:]
            trial[column] += step
            shifted, _ = evaluate(trial)
            for row in range(size):
                jacobian[row][column] = (shifted[row] - current[row]) / step
        direction = _solve_linear(jacobian, [-value for value in current])
        scale = 1.0
        accepted = False
        while scale >= 2.0**-24:
            trial = [value + scale * delta for value, delta in zip(variables, direction)]
            try:
                shifted, _ = evaluate(trial)
            except ValidationError:
                scale *= 0.5
                continue
            if max(abs(value) for value in shifted) < norm:
                variables = trial
                accepted = True
                break
            scale *= 0.5
        if not accepted:
            raise ConvergenceError("Gibbs reactor line search did not reduce the residual")
    else:
        raise ConvergenceError("Gibbs reactor did not converge")

    if outlet_vector is None or final_residual is None:
        raise ConvergenceError("Gibbs reactor produced no converged state")
    outlet = tuple(zip(compound_ids, outlet_vector))
    conversions = tuple(
        (compound_id, (inlet[compound_id] - outlet_vector[index]) / inlet[compound_id])
        for index, compound_id in enumerate(compound_ids)
        if inlet[compound_id] > 0.0 and outlet_vector[index] < inlet[compound_id]
    )
    element_residuals = tuple(
        (
            element,
            math.fsum(coefficient * flow for coefficient, flow in zip(row, outlet_vector)) - total,
        )
        for element, row, total in zip(elements, element_matrix, element_totals)
    )

    def gibbs_energy_flow(component_flows: tuple[float, ...]) -> float:
        potentials = chemical_potentials_rt(tuple(max(flow, floor) for flow in component_flows))
        return DWSIM_EQUILIBRIUM_R * temperature_k * math.fsum(
            flow * potential for flow, potential in zip(component_flows, potentials) if flow > 0.0
        )

    def enthalpy_flow(component_flows: tuple[float, ...]) -> float:
        total = math.fsum(component_flows)
        composition = tuple(flow / total for flow in component_flows)
        state = PengRobinsonMixture(compounds, composition, interactions).state(
            temperature_k, pressure_pa, "vapor",
        )
        ideal = math.fsum(
            flow * (
                thermo_by_id[compound_id].ideal_gas_formation_enthalpy_j_per_kmol
                + _dwsim_midpoint_integrals(
                    correlation_by_id[compound_id],
                    thermo_by_id[compound_id].formation_temperature_k,
                    temperature_k,
                )[0]
            )
            for compound_id, flow in zip(compound_ids, component_flows)
        )
        return ideal + total * state.departure_enthalpy_j_per_kmol

    return GibbsReactorResult(
        outlet,
        conversions,
        element_residuals,
        max(abs(value) for value in final_residual[:len(compound_ids)]),
        gibbs_energy_flow(inlet_vector),
        gibbs_energy_flow(outlet_vector),
        enthalpy_flow(outlet_vector) - enthalpy_flow(inlet_vector),
        iteration,
    )
