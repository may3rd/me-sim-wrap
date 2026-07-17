"""Reaction unit operations."""
import math
from dataclasses import dataclass

from ..errors import ValidationError
from ..reactions import ReactionDefinition


@dataclass(frozen=True, slots=True)
class ConversionReactorResult:
    outlet_component_flows_kmol_s: tuple[tuple[str, float], ...]
    extent_kmol_s: float
    conversion_fraction: float
    reaction_heat_w: float

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
