import math

from ..compounds import Compound, PRInteractions
from ..errors import ConvergenceError, ValidationError
from ..streams import PhaseState, StreamState
from ..thermo.flash import ph_flash
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
