"""Pressure-changing unit operations."""
import math
from dataclasses import dataclass

from ..compounds import Compound, PRInteractions
from ..errors import ConvergenceError, ValidationError
from ..streams import EnergyStream, PhaseState, StreamState
from ..thermo.flash import flash_enthalpy, flash_entropy, ph_flash, ps_flash
from ..thermo.ideal import IdealCorrelations
from ..thermo.peng_robinson import R


@dataclass(frozen=True, slots=True)
class PumpResult:
    outlet: PhaseState
    energy: EnergyStream


@dataclass(frozen=True, slots=True)
class CompressorResult:
    outlet: PhaseState
    energy: EnergyStream


@dataclass(frozen=True, slots=True)
class ExpanderResult:
    outlet: PhaseState
    energy: EnergyStream


def pump(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    efficiency: float,
    temperature_bracket_k: tuple[float, float],
) -> PumpResult:
    """Adiabatically raise a liquid stream's pressure at constant efficiency."""
    if isinstance(outlet_pressure_pa, bool) or not isinstance(outlet_pressure_pa, (int, float)) or not math.isfinite(outlet_pressure_pa) or outlet_pressure_pa <= inlet.stream.pressure_pa:
        raise ValidationError("pump outlet pressure must be finite and above inlet pressure")
    if isinstance(efficiency, bool) or not isinstance(efficiency, (int, float)) or not math.isfinite(efficiency) or not 0.0 < efficiency <= 1.0:
        raise ValidationError("pump efficiency must be finite and within (0, 1]")
    if inlet.flash.phase != "liquid" or inlet.flash.liquid_state is None:
        raise ValidationError("pump requires a liquid inlet")
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("pump inlet compound IDs must exactly match the supplied compound order")

    molar_volume = (
        inlet.flash.liquid_state.compressibility * R * inlet.stream.temperature_k / inlet.stream.pressure_pa
    )
    target_enthalpy = inlet.enthalpy_j_per_kmol + molar_volume * (outlet_pressure_pa - inlet.stream.pressure_pa) / efficiency
    result = ph_flash(
        compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa,
        target_enthalpy, temperature_bracket_k,
    )
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "pump PH flash did not converge")
    outlet = PhaseState(
        StreamState(
            result.temperature_k, outlet_pressure_pa, inlet.stream.molar_flow_kmol_s,
            inlet.stream.compound_ids, inlet.stream.composition,
            result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
        ),
        result.flash,
        result.enthalpy_j_per_kmol,
    )
    return PumpResult(outlet, EnergyStream(inlet.stream.molar_flow_kmol_s * (outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol)))


def compressor(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    efficiency: float,
    temperature_bracket_k: tuple[float, float],
) -> CompressorResult:
    """Adiabatically compress a vapor stream at constant isentropic efficiency."""
    if isinstance(outlet_pressure_pa, bool) or not isinstance(outlet_pressure_pa, (int, float)) or not math.isfinite(outlet_pressure_pa) or outlet_pressure_pa <= inlet.stream.pressure_pa:
        raise ValidationError("compressor outlet pressure must be finite and above inlet pressure")
    if isinstance(efficiency, bool) or not isinstance(efficiency, (int, float)) or not math.isfinite(efficiency) or not 0.0 < efficiency <= 1.0:
        raise ValidationError("compressor efficiency must be finite and within (0, 1]")
    if inlet.flash.phase not in {"vapor", "single"} or inlet.flash.vapor_state is None:
        raise ValidationError("compressor requires a vapor or single-root inlet")
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("compressor inlet compound IDs must exactly match the supplied compound order")
    entropy = flash_entropy(compounds, correlations, inlet.flash)
    isentropic = ps_flash(compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa, entropy, temperature_bracket_k)
    if not isentropic.report.converged or isentropic.flash is None:
        raise ConvergenceError(isentropic.report.failure_reason or "compressor PS flash did not converge")
    target_enthalpy = inlet.enthalpy_j_per_kmol + (flash_enthalpy(compounds, correlations, isentropic.flash) - inlet.enthalpy_j_per_kmol) / efficiency
    result = ph_flash(compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa, target_enthalpy, temperature_bracket_k)
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "compressor PH flash did not converge")
    outlet = PhaseState(
        StreamState(result.temperature_k, outlet_pressure_pa, inlet.stream.molar_flow_kmol_s, inlet.stream.compound_ids, inlet.stream.composition, result.enthalpy_j_per_kmol, result.flash.vapor_fraction),
        result.flash,
        result.enthalpy_j_per_kmol,
    )
    return CompressorResult(outlet, EnergyStream(inlet.stream.molar_flow_kmol_s * (outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol)))


def expander(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    efficiency: float,
    temperature_bracket_k: tuple[float, float],
) -> ExpanderResult:
    """Adiabatically expand a vapor stream at constant isentropic efficiency."""
    if isinstance(outlet_pressure_pa, bool) or not isinstance(outlet_pressure_pa, (int, float)) or not math.isfinite(outlet_pressure_pa) or outlet_pressure_pa >= inlet.stream.pressure_pa:
        raise ValidationError("expander outlet pressure must be finite and below inlet pressure")
    if isinstance(efficiency, bool) or not isinstance(efficiency, (int, float)) or not math.isfinite(efficiency) or not 0.0 < efficiency <= 1.0:
        raise ValidationError("expander efficiency must be finite and within (0, 1]")
    if inlet.flash.phase not in {"vapor", "single"} or inlet.flash.vapor_state is None:
        raise ValidationError("expander requires a vapor or single-root inlet")
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("expander inlet compound IDs must exactly match the supplied compound order")
    entropy = flash_entropy(compounds, correlations, inlet.flash)
    isentropic = ps_flash(compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa, entropy, temperature_bracket_k)
    if not isentropic.report.converged or isentropic.flash is None:
        raise ConvergenceError(isentropic.report.failure_reason or "expander PS flash did not converge")
    target_enthalpy = inlet.enthalpy_j_per_kmol - efficiency * (inlet.enthalpy_j_per_kmol - flash_enthalpy(compounds, correlations, isentropic.flash))
    result = ph_flash(compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa, target_enthalpy, temperature_bracket_k)
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "expander PH flash did not converge")
    outlet = PhaseState(
        StreamState(result.temperature_k, outlet_pressure_pa, inlet.stream.molar_flow_kmol_s, inlet.stream.compound_ids, inlet.stream.composition, result.enthalpy_j_per_kmol, result.flash.vapor_fraction),
        result.flash,
        result.enthalpy_j_per_kmol,
    )
    return ExpanderResult(outlet, EnergyStream(inlet.stream.molar_flow_kmol_s * (outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol)))
