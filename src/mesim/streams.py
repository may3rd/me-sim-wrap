import math
from dataclasses import dataclass

from .compounds import Compound, PRInteractions
from .errors import ConvergenceError, ValidationError
from .thermo.flash import TPFlashResult
from .thermo.ideal import IdealCorrelations
from .thermo.systems import PengRobinsonSystem


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
class EnergyStream:
    """Energy duty in W; positive duty enters the connected material stream."""
    duty_w: float

    def __post_init__(self) -> None:
        if isinstance(self.duty_w, bool) or not isinstance(self.duty_w, (int, float)) or not math.isfinite(self.duty_w):
            raise ValidationError("energy stream duty must be a finite numeric value")


@dataclass(frozen=True, slots=True)
class PhaseState:
    stream: StreamState
    flash: TPFlashResult
    enthalpy_j_per_kmol: float


def flash_stream(
    stream: StreamState,
    system: PengRobinsonSystem | tuple[Compound, ...],
    interactions: PRInteractions | None = None,
    correlations: tuple[IdealCorrelations, ...] | None = None,
) -> PhaseState:
    """Flash a stream through an extracted PR system.

    The three-data-argument form remains temporarily accepted for internal
    callers while they migrate to :class:`PengRobinsonSystem`.
    """
    if isinstance(system, PengRobinsonSystem):
        if interactions is not None or correlations is not None:
            raise ValidationError("a PR system cannot be combined with separate thermodynamic data")
        thermo = system
    else:
        if interactions is None or correlations is None:
            raise ValidationError("legacy stream flashing requires PR interactions and ideal correlations")
        thermo = PengRobinsonSystem(tuple(system), interactions, tuple(correlations))
    if thermo.compound_ids != stream.compound_ids:
        raise ValidationError("stream compound IDs must exactly match the supplied compound order")
    flash = thermo.tp_flash(
        stream.composition, stream.temperature_k, stream.pressure_pa
    )
    if not flash.report.converged:
        raise ConvergenceError(flash.report.failure_reason or "TP flash did not converge")
    return PhaseState(stream, flash, thermo.enthalpy(flash))
