import math
from dataclasses import dataclass

from ..compounds import Compound, PRInteractions
from ..errors import ConvergenceError, ValidationError
from ..streams import EnergyStream, PhaseState, StreamState, flash_stream
from ..thermo.flash import ph_flash, phase_enthalpy
from ..thermo.ideal import IdealCorrelations


def _positive_finite(value: float, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0:
        raise ValidationError(f"{name} must be finite and positive")


def mix_streams(
    inlets: tuple[PhaseState, ...],
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    temperature_bracket_k: tuple[float, float],
) -> PhaseState:
    """Adiabatically mix flashed streams at a specified pressure."""
    if not inlets:
        raise ValidationError("mixer requires at least one inlet")
    _positive_finite(outlet_pressure_pa, "mixer outlet pressure")
    compound_ids = tuple(compound.id for compound in compounds)
    if not compound_ids or any(inlet.stream.compound_ids != compound_ids for inlet in inlets):
        raise ValidationError("mixer inlet compound IDs must exactly match the supplied compound order")
    if outlet_pressure_pa > min(inlet.stream.pressure_pa for inlet in inlets):
        raise ValidationError("mixer outlet pressure must not exceed the lowest inlet pressure")

    total_flow = math.fsum(inlet.stream.molar_flow_kmol_s for inlet in inlets)
    if total_flow == 0.0:
        raise ValidationError("mixer requires a positive total inlet flow")
    composition = tuple(
        math.fsum(inlet.stream.molar_flow_kmol_s * inlet.stream.composition[index] for inlet in inlets) / total_flow
        for index in range(len(compound_ids))
    )
    target_enthalpy = math.fsum(
        inlet.stream.molar_flow_kmol_s * inlet.enthalpy_j_per_kmol for inlet in inlets
    ) / total_flow
    result = ph_flash(
        compounds, composition, interactions, correlations, outlet_pressure_pa, target_enthalpy, temperature_bracket_k,
    )
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "mixer PH flash did not converge")
    stream = StreamState(
        result.temperature_k, outlet_pressure_pa, total_flow, compound_ids, composition,
        result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
    )
    return PhaseState(stream, result.flash, result.enthalpy_j_per_kmol)


def split_stream(stream: StreamState, fractions: tuple[float, ...]) -> tuple[StreamState, ...]:
    """Split a stream without changing its intensive state or composition."""
    if not fractions or any(
        isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or fraction < 0.0
        for fraction in fractions
    ) or not math.isclose(math.fsum(fractions), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("splitter fractions must be finite, non-negative, and sum to one")
    return tuple(
        StreamState(
            stream.temperature_k, stream.pressure_pa, stream.molar_flow_kmol_s * fraction,
            stream.compound_ids, stream.composition, stream.enthalpy_j_per_kmol, stream.vapor_fraction,
        )
        for fraction in fractions
    )


@dataclass(frozen=True, slots=True)
class ThermalOperationResult:
    outlet: PhaseState
    energy: EnergyStream


@dataclass(frozen=True, slots=True)
class HeatExchangerResult:
    hot_outlet: PhaseState
    cold_outlet: PhaseState
    heat_duty_w: float


def heat_exchanger(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    heat_duty_w: float,
    temperature_bracket_k: tuple[float, float],
) -> HeatExchangerResult:
    """Transfer a specified positive duty from hot inlet to cold inlet at constant pressure."""
    _positive_finite(heat_duty_w, "heat exchanger duty")
    compound_ids = tuple(compound.id for compound in compounds)
    if hot_inlet.stream.compound_ids != compound_ids or cold_inlet.stream.compound_ids != compound_ids:
        raise ValidationError("heat exchanger inlet compound IDs must exactly match the supplied compound order")
    _positive_finite(hot_inlet.stream.molar_flow_kmol_s, "heat exchanger hot inlet flow")
    _positive_finite(cold_inlet.stream.molar_flow_kmol_s, "heat exchanger cold inlet flow")

    def outlet(inlet: PhaseState, target_enthalpy: float) -> PhaseState:
        result = ph_flash(
            compounds, inlet.stream.composition, interactions, correlations, inlet.stream.pressure_pa,
            target_enthalpy, temperature_bracket_k,
        )
        if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
            raise ConvergenceError(result.report.failure_reason or "heat exchanger PH flash did not converge")
        stream = StreamState(
            result.temperature_k, inlet.stream.pressure_pa, inlet.stream.molar_flow_kmol_s,
            inlet.stream.compound_ids, inlet.stream.composition,
            result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
        )
        return PhaseState(stream, result.flash, result.enthalpy_j_per_kmol)

    return HeatExchangerResult(
        outlet(hot_inlet, hot_inlet.enthalpy_j_per_kmol - heat_duty_w / hot_inlet.stream.molar_flow_kmol_s),
        outlet(cold_inlet, cold_inlet.enthalpy_j_per_kmol + heat_duty_w / cold_inlet.stream.molar_flow_kmol_s),
        heat_duty_w,
    )


def heat_exchanger_ua(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    overall_coefficient_w_m2_k: float,
    area_m2: float,
    temperature_bracket_k: tuple[float, float],
) -> HeatExchangerResult:
    """Calculate a counter-current, no-loss exchanger from a specified UA."""
    _positive_finite(overall_coefficient_w_m2_k, "heat exchanger overall coefficient")
    _positive_finite(area_m2, "heat exchanger area")
    if hot_inlet.stream.temperature_k <= cold_inlet.stream.temperature_k:
        raise ValidationError("heat exchanger hot inlet temperature must exceed cold inlet temperature")
    ua = overall_coefficient_w_m2_k * area_m2

    def lmtd(result: HeatExchangerResult) -> float | None:
        first = hot_inlet.stream.temperature_k - result.cold_outlet.stream.temperature_k
        second = result.hot_outlet.stream.temperature_k - cold_inlet.stream.temperature_k
        if first <= 0.0 or second <= 0.0:
            return None
        if math.isclose(first, second, rel_tol=1e-12, abs_tol=1e-12):
            return (first + second) / 2.0
        return (first - second) / math.log(first / second)

    lower = 0.0
    upper = ua * (hot_inlet.stream.temperature_k - cold_inlet.stream.temperature_k)
    for _ in range(80):
        trial = heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, upper, temperature_bracket_k)
        difference = lmtd(trial)
        if difference is None or upper >= ua * difference:
            break
        upper *= 2.0
    else:
        raise ConvergenceError("heat exchanger UA duty bracket did not converge")

    for _ in range(60):
        duty = (lower + upper) / 2.0
        trial = heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, duty, temperature_bracket_k)
        difference = lmtd(trial)
        if difference is None or duty >= ua * difference:
            upper = duty
        else:
            lower = duty
    return heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, (lower + upper) / 2.0, temperature_bracket_k)


def _thermal_operation(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_temperature_k: float,
    kind: str,
) -> ThermalOperationResult:
    _positive_finite(outlet_temperature_k, f"{kind} outlet temperature")
    outlet = flash_stream(
        StreamState(
            outlet_temperature_k, inlet.stream.pressure_pa, inlet.stream.molar_flow_kmol_s,
            inlet.stream.compound_ids, inlet.stream.composition,
        ),
        compounds, interactions, correlations,
    )
    duty = inlet.stream.molar_flow_kmol_s * (outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol)
    if (kind == "heater" and duty < 0.0) or (kind == "cooler" and duty > 0.0):
        raise ValidationError(f"{kind} outlet temperature gives an invalid duty sign")
    return ThermalOperationResult(outlet, EnergyStream(duty))


def heater(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_temperature_k: float,
) -> ThermalOperationResult:
    """Heat a flashed stream at constant pressure to a specified temperature."""
    return _thermal_operation(inlet, compounds, interactions, correlations, outlet_temperature_k, "heater")


def cooler(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_temperature_k: float,
) -> ThermalOperationResult:
    """Cool a flashed stream at constant pressure to a specified temperature."""
    return _thermal_operation(inlet, compounds, interactions, correlations, outlet_temperature_k, "cooler")


def valve(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    temperature_bracket_k: tuple[float, float],
) -> PhaseState:
    """Throttle a flashed stream isenthalpically to a lower pressure."""
    _positive_finite(outlet_pressure_pa, "valve outlet pressure")
    if outlet_pressure_pa >= inlet.stream.pressure_pa:
        raise ValidationError("valve outlet pressure must be below inlet pressure")
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("valve inlet compound IDs must exactly match the supplied compound order")
    result = ph_flash(
        compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa,
        inlet.enthalpy_j_per_kmol, temperature_bracket_k,
    )
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "valve PH flash did not converge")
    stream = StreamState(
        result.temperature_k, outlet_pressure_pa, inlet.stream.molar_flow_kmol_s, inlet.stream.compound_ids,
        inlet.stream.composition, result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
    )
    return PhaseState(stream, result.flash, result.enthalpy_j_per_kmol)


@dataclass(frozen=True, slots=True)
class EquilibriumSeparatorResult:
    liquid: StreamState | None
    vapor: StreamState | None


def equilibrium_separator(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    correlations: tuple[IdealCorrelations, ...],
) -> EquilibriumSeparatorResult:
    """Route phases from an existing converged flash without recalculation."""
    flash = inlet.flash
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("separator inlet compound IDs must exactly match the supplied compound order")
    if not flash.report.converged:
        raise ConvergenceError(flash.report.failure_reason or "separator requires a converged flash")

    def phase_stream(composition: tuple[float, ...], state, flow: float, vapor_fraction: float) -> StreamState:
        return StreamState(
            flash.temperature_k, flash.pressure_pa, flow, inlet.stream.compound_ids, composition,
            phase_enthalpy(compounds, composition, correlations, flash.temperature_k, state), vapor_fraction,
        )

    if flash.phase == "liquid":
        return EquilibriumSeparatorResult(
            phase_stream(flash.liquid_composition, flash.liquid_state, inlet.stream.molar_flow_kmol_s, 0.0), None,
        )
    if flash.phase == "vapor":
        return EquilibriumSeparatorResult(
            None, phase_stream(flash.vapor_composition, flash.vapor_state, inlet.stream.molar_flow_kmol_s, 1.0),
        )
    if flash.phase == "two-phase" and flash.vapor_fraction is not None:
        return EquilibriumSeparatorResult(
            phase_stream(
                flash.liquid_composition, flash.liquid_state,
                inlet.stream.molar_flow_kmol_s * (1.0 - flash.vapor_fraction), 0.0,
            ),
            phase_stream(
                flash.vapor_composition, flash.vapor_state,
                inlet.stream.molar_flow_kmol_s * flash.vapor_fraction, 1.0,
            ),
        )
    raise ValidationError(f"unsupported separator flash phase: {flash.phase}")
