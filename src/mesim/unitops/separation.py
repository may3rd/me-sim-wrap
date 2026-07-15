"""Component-separation unit operations."""
import math
from dataclasses import dataclass

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError
from ..streams import EnergyStream, PhaseState, StreamState, flash_stream
from ..thermo.ideal import IdealCorrelations


@dataclass(frozen=True, slots=True)
class ComponentSeparatorResult:
    specified: PhaseState
    remainder: PhaseState
    energy: EnergyStream


def component_separator(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    specified_component_fractions: tuple[float, ...],
) -> ComponentSeparatorResult:
    """Split each component to the specified outlet at unchanged temperature and pressure."""
    if len(specified_component_fractions) != len(inlet.stream.composition) or any(
        isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0
        for fraction in specified_component_fractions
    ):
        raise ValidationError("component separator fractions must match components and be within [0, 1]")
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("component separator inlet compound IDs must exactly match the supplied compound order")

    def outlet(fractions: tuple[float, ...]) -> PhaseState:
        flows = tuple(inlet.stream.molar_flow_kmol_s * composition * fraction for composition, fraction in zip(inlet.stream.composition, fractions))
        total = math.fsum(flows)
        composition = tuple(flow / total for flow in flows) if total else inlet.stream.composition
        return flash_stream(StreamState(inlet.stream.temperature_k, inlet.stream.pressure_pa, total, inlet.stream.compound_ids, composition), compounds, interactions, correlations)

    specified = outlet(specified_component_fractions)
    remainder = outlet(tuple(1.0 - fraction for fraction in specified_component_fractions))
    duty = (
        specified.stream.molar_flow_kmol_s * specified.enthalpy_j_per_kmol
        + remainder.stream.molar_flow_kmol_s * remainder.enthalpy_j_per_kmol
        - inlet.stream.molar_flow_kmol_s * inlet.enthalpy_j_per_kmol
    )
    return ComponentSeparatorResult(specified, remainder, EnergyStream(duty))
