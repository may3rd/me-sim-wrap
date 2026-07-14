import math
from dataclasses import dataclass

from .compounds import Compound, PRInteractions
from .errors import ConvergenceError, ValidationError
from .thermo.flash import TPFlashResult, flash_enthalpy, tp_flash
from .thermo.ideal import IdealCorrelations


@dataclass(frozen=True, slots=True)
class StreamState:
    temperature_k: float
    pressure_pa: float
    molar_flow_kmol_s: float
    compound_ids: tuple[str, ...]
    composition: tuple[float, ...]
    enthalpy_j_per_kmol: float | None = None
    vapor_fraction: float | None = None

    def __post_init__(self) -> None:
        if any(
            isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            for value in (self.temperature_k, self.pressure_pa, self.molar_flow_kmol_s)
        ):
            raise ValidationError("stream temperature, pressure, and molar flow must be finite numeric values")
        if self.temperature_k <= 0.0 or self.pressure_pa <= 0.0 or self.molar_flow_kmol_s < 0.0:
            raise ValidationError("stream temperature and pressure must be positive and molar flow non-negative")
        if not self.compound_ids or len(self.compound_ids) != len(self.composition):
            raise ValidationError("stream compound IDs and composition must have the same nonzero length")
        if len(set(self.compound_ids)) != len(self.compound_ids) or any(not isinstance(value, str) or not value for value in self.compound_ids):
            raise ValidationError("stream compound IDs must be unique non-empty strings")
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0.0 for value in self.composition):
            raise ValidationError("stream composition must be finite and non-negative")
        if not math.isclose(math.fsum(self.composition), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValidationError("stream composition must sum to one")
        for name, value in (("enthalpy", self.enthalpy_j_per_kmol), ("vapor fraction", self.vapor_fraction)):
            if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)):
                raise ValidationError(f"stream {name} must be finite when supplied")
        if self.vapor_fraction is not None and not 0.0 <= self.vapor_fraction <= 1.0:
            raise ValidationError("stream vapor fraction must be between zero and one")


@dataclass(frozen=True, slots=True)
class PhaseState:
    stream: StreamState
    flash: TPFlashResult
    enthalpy_j_per_kmol: float


def flash_stream(
    stream: StreamState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
) -> PhaseState:
    if tuple(compound.id for compound in compounds) != stream.compound_ids:
        raise ValidationError("stream compound IDs must exactly match the supplied compound order")
    flash = tp_flash(compounds, stream.composition, interactions, stream.temperature_k, stream.pressure_pa)
    if not flash.report.converged:
        raise ConvergenceError(flash.report.failure_reason or "TP flash did not converge")
    return PhaseState(stream, flash, flash_enthalpy(compounds, correlations, flash))
