from dataclasses import dataclass
import math

from ..errors import ValidationError


def _finite_number(value: object) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
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
