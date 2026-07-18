"""Explicit thermodynamic-system boundaries and stable model registry."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Protocol, runtime_checkable

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError
from .activity import (
    NRTLPhaseEnthalpies,
    NRTLVLEData,
    NRTLVLEResult,
    NRTLCaloricRecord,
    nrtl_bubble_pressure,
    nrtl_bubble_temperature,
    nrtl_dew_pressure,
    nrtl_equilibrium_ratios,
    nrtl_phase_enthalpies,
)
from .flash import (
    PHFlashResult,
    PSFlashResult,
    TPFlashResult,
    flash_enthalpy,
    ph_flash,
    ps_flash,
    tp_flash,
)
from .ideal import IdealCorrelations
from .pure import SaturatedLiquidCorrelations
from .raoult import (
    RaoultTPFlashResult,
    RaoultVLEResult,
    raoult_bubble_pressure,
    raoult_dew_pressure,
    raoult_equilibrium_ratios,
    raoult_fugacity_coefficients,
    raoult_tp_flash,
)
from .soave_redlich_kwong import (
    SRKMixtureState,
    SRKTPFlashResult,
    SoaveRedlichKwongMixture,
    srk_tp_flash,
)
from .transport import TransportRecord


PENG_ROBINSON_CLASSIC = "peng-robinson-classic"
NRTL_ACETONE_METHANOL = "nrtl-acetone-methanol"
IDEAL_RAOULT = "ideal-raoult"
SOAVE_REDLICH_KWONG = "soave-redlich-kwong"


@runtime_checkable
class ThermodynamicSystem(Protocol):
    """Identity shared by systems with deliberately different capabilities."""

    model_id: str
    compound_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PengRobinsonSystem:
    """Classic-PR flash and caloric system for one ordered compound domain."""

    compounds: tuple[Compound, ...]
    interactions: PRInteractions
    correlations: tuple[IdealCorrelations, ...]
    transport_correlations: tuple[TransportRecord, ...] = ()
    saturated_liquid_correlations: tuple[SaturatedLiquidCorrelations, ...] = ()
    model_id: str = field(default=PENG_ROBINSON_CLASSIC, init=False)
    compound_ids: tuple[str, ...] = field(init=False)

    def __post_init__(self) -> None:
        try:
            compounds = tuple(self.compounds)
            correlations = tuple(self.correlations)
            transport = tuple(self.transport_correlations)
            saturated = tuple(self.saturated_liquid_correlations)
        except TypeError as error:
            raise ValidationError("PR thermodynamic-system inputs must be sequences") from error
        if (
            not compounds
            or not isinstance(self.interactions, PRInteractions)
            or self.interactions.model != "Peng-Robinson"
            or any(not isinstance(record, Compound) for record in compounds)
            or any(not isinstance(record, IdealCorrelations) for record in correlations)
            or any(not isinstance(record, TransportRecord) for record in transport)
            or any(not isinstance(record, SaturatedLiquidCorrelations) for record in saturated)
        ):
            raise ValidationError("PR thermodynamic-system inputs are invalid")
        compound_ids = tuple(record.id for record in compounds)
        if len(set(compound_ids)) != len(compound_ids):
            raise ValidationError("PR thermodynamic-system compound IDs must be unique")
        correlation_ids = {record.compound_id for record in correlations}
        if len(correlation_ids) != len(correlations):
            raise ValidationError("PR thermodynamic-system correlation IDs must be unique")
        missing = tuple(name for name in compound_ids if name not in correlation_ids)
        if missing:
            raise ValidationError(
                f"PR thermodynamic system is missing ideal correlations: {', '.join(missing)}"
            )
        for label, records in (
            ("transport", transport),
            ("saturated-liquid", saturated),
        ):
            record_ids = tuple(record.compound_id for record in records)
            if len(set(record_ids)) != len(record_ids):
                raise ValidationError(
                    f"PR thermodynamic-system {label} correlation IDs must be unique"
                )
            missing = tuple(name for name in compound_ids if name not in record_ids)
            if records and missing:
                raise ValidationError(
                    f"PR thermodynamic system is missing {label} correlations: "
                    f"{', '.join(missing)}"
                )
        for first_index, first in enumerate(compound_ids):
            for second in compound_ids[first_index + 1:]:
                self.interactions.get(first, second)
        object.__setattr__(self, "compounds", compounds)
        object.__setattr__(self, "correlations", correlations)
        object.__setattr__(self, "transport_correlations", transport)
        object.__setattr__(self, "saturated_liquid_correlations", saturated)
        object.__setattr__(self, "compound_ids", compound_ids)

    def ideal(self, compound_id: str) -> IdealCorrelations:
        return self._record(self.correlations, compound_id, "ideal")

    def transport(self, compound_id: str) -> TransportRecord:
        return self._record(self.transport_correlations, compound_id, "transport")

    def saturated_liquid(self, compound_id: str) -> SaturatedLiquidCorrelations:
        return self._record(
            self.saturated_liquid_correlations, compound_id, "saturated-liquid"
        )

    def _record(self, records, compound_id: str, label: str):
        if not isinstance(compound_id, str) or compound_id not in self.compound_ids:
            raise ValidationError(
                f"compound is outside PR thermodynamic-system domain: {compound_id}"
            )
        try:
            return next(record for record in records if record.compound_id == compound_id)
        except StopIteration as error:
            raise ValidationError(
                f"PR thermodynamic system has no {label} correlation for {compound_id}"
            ) from error

    def tp_flash(
        self,
        composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
        *,
        max_iterations: int = 100,
    ) -> TPFlashResult:
        return tp_flash(
            self.compounds,
            composition,
            self.interactions,
            temperature_k,
            pressure_pa,
            max_iterations=max_iterations,
        )

    def ph_flash(
        self,
        composition: tuple[float, ...],
        pressure_pa: float,
        target_enthalpy_j_per_kmol: float,
        temperature_bracket_k: tuple[float, float],
        *,
        max_iterations: int = 100,
    ) -> PHFlashResult:
        return ph_flash(
            self.compounds,
            composition,
            self.interactions,
            self.correlations,
            pressure_pa,
            target_enthalpy_j_per_kmol,
            temperature_bracket_k,
            max_iterations=max_iterations,
        )

    def ps_flash(
        self,
        composition: tuple[float, ...],
        pressure_pa: float,
        target_entropy_j_per_kmol_k: float,
        temperature_bracket_k: tuple[float, float],
        *,
        max_iterations: int = 100,
    ) -> PSFlashResult:
        return ps_flash(
            self.compounds,
            composition,
            self.interactions,
            self.correlations,
            pressure_pa,
            target_entropy_j_per_kmol_k,
            temperature_bracket_k,
            max_iterations=max_iterations,
        )

    def enthalpy(self, flash: TPFlashResult) -> float:
        return flash_enthalpy(self.compounds, self.correlations, flash)


@dataclass(frozen=True, slots=True)
class NRTLSystem:
    """Accepted binary NRTL stage-equilibrium and caloric system."""

    data: NRTLVLEData
    compound_ids: tuple[str, ...]
    model_id: str = field(default=NRTL_ACETONE_METHANOL, init=False)

    def __post_init__(self) -> None:
        try:
            compound_ids = tuple(self.compound_ids)
        except TypeError as error:
            raise ValidationError("NRTL thermodynamic-system compound IDs must be a sequence") from error
        if (
            not isinstance(self.data, NRTLVLEData)
            or len(compound_ids) < 2
            or len(set(compound_ids)) != len(compound_ids)
            or any(not isinstance(value, str) or not value for value in compound_ids)
        ):
            raise ValidationError("NRTL thermodynamic-system inputs are invalid")
        for compound_id in compound_ids:
            self.data.vapor_pressure(compound_id)
            self.data.caloric(compound_id)
        for first in compound_ids:
            for second in compound_ids:
                if first != second:
                    self.data.interaction(first, second)
        object.__setattr__(self, "compound_ids", compound_ids)

    def caloric(self, compound_id: str) -> NRTLCaloricRecord:
        return self.data.caloric(compound_id)

    def equilibrium_ratios(
        self,
        liquid_composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
    ) -> tuple[float, ...]:
        return nrtl_equilibrium_ratios(
            self.data,
            self.compound_ids,
            liquid_composition,
            temperature_k,
            pressure_pa,
        )

    def bubble_pressure(
        self,
        liquid_composition: tuple[float, ...],
        temperature_k: float,
    ) -> NRTLVLEResult:
        return nrtl_bubble_pressure(
            self.data, self.compound_ids, liquid_composition, temperature_k
        )

    def dew_pressure(
        self,
        vapor_composition: tuple[float, ...],
        temperature_k: float,
        *,
        max_iterations: int = 100,
        tolerance: float = 1.0e-12,
    ) -> NRTLVLEResult:
        return nrtl_dew_pressure(
            self.data,
            self.compound_ids,
            vapor_composition,
            temperature_k,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )

    def bubble_temperature(
        self,
        liquid_composition: tuple[float, ...],
        pressure_pa: float,
        bracket_k: tuple[float, float],
        *,
        max_iterations: int = 100,
        tolerance: float = 1.0e-10,
    ) -> NRTLVLEResult:
        return nrtl_bubble_temperature(
            self.data,
            self.compound_ids,
            liquid_composition,
            pressure_pa,
            bracket_k,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )

    def phase_enthalpies(
        self,
        liquid_composition: tuple[float, ...],
        vapor_composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
    ) -> NRTLPhaseEnthalpies:
        return nrtl_phase_enthalpies(
            self.data,
            self.compound_ids,
            liquid_composition,
            vapor_composition,
            temperature_k,
            pressure_pa,
        )


@dataclass(frozen=True, slots=True)
class IdealRaoultSystem:
    """Ideal Raoult-law equilibrium over an ordered pure-correlation domain."""

    correlations: tuple[IdealCorrelations, ...]
    model_id: str = field(default=IDEAL_RAOULT, init=False)
    compound_ids: tuple[str, ...] = field(init=False)

    def __post_init__(self) -> None:
        try:
            correlations = tuple(self.correlations)
        except TypeError as error:
            raise ValidationError(
                "Raoult thermodynamic-system correlations must be a sequence"
            ) from error
        if (
            not correlations
            or any(not isinstance(record, IdealCorrelations) for record in correlations)
        ):
            raise ValidationError("Raoult thermodynamic-system correlations are invalid")
        compound_ids = tuple(record.compound_id for record in correlations)
        if (
            len(set(compound_ids)) != len(compound_ids)
            or any(not isinstance(value, str) or not value for value in compound_ids)
        ):
            raise ValidationError(
                "Raoult thermodynamic-system compound IDs must be non-empty and unique"
            )
        object.__setattr__(self, "correlations", correlations)
        object.__setattr__(self, "compound_ids", compound_ids)

    def equilibrium_ratios(
        self,
        temperature_k: float,
        pressure_pa: float,
        *,
        allow_extrapolation: bool = False,
    ) -> tuple[float, ...]:
        return raoult_equilibrium_ratios(
            self.correlations,
            temperature_k,
            pressure_pa,
            allow_extrapolation=allow_extrapolation,
        )

    def fugacity_coefficients(
        self,
        temperature_k: float,
        pressure_pa: float,
        phase: str,
        *,
        allow_extrapolation: bool = False,
    ) -> tuple[float, ...]:
        return raoult_fugacity_coefficients(
            self.correlations,
            temperature_k,
            pressure_pa,
            phase,
            allow_extrapolation=allow_extrapolation,
        )

    def bubble_pressure(
        self,
        liquid_composition: tuple[float, ...],
        temperature_k: float,
        *,
        allow_extrapolation: bool = False,
    ) -> RaoultVLEResult:
        return raoult_bubble_pressure(
            self.correlations,
            liquid_composition,
            temperature_k,
            allow_extrapolation=allow_extrapolation,
        )

    def dew_pressure(
        self,
        vapor_composition: tuple[float, ...],
        temperature_k: float,
        *,
        allow_extrapolation: bool = False,
    ) -> RaoultVLEResult:
        return raoult_dew_pressure(
            self.correlations,
            vapor_composition,
            temperature_k,
            allow_extrapolation=allow_extrapolation,
        )

    def tp_flash(
        self,
        composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
        *,
        max_iterations: int = 100,
        tolerance: float = 1.0e-12,
        allow_extrapolation: bool = False,
    ) -> RaoultTPFlashResult:
        return raoult_tp_flash(
            self.correlations,
            composition,
            temperature_k,
            pressure_pa,
            max_iterations=max_iterations,
            tolerance=tolerance,
            allow_extrapolation=allow_extrapolation,
        )


@dataclass(frozen=True, slots=True)
class SoaveRedlichKwongSystem:
    """Classic SRK phase-state boundary for one ordered compound domain."""

    compounds: tuple[Compound, ...]
    interactions: PRInteractions
    model_id: str = field(default=SOAVE_REDLICH_KWONG, init=False)
    compound_ids: tuple[str, ...] = field(init=False)

    def __post_init__(self) -> None:
        try:
            compounds = tuple(self.compounds)
        except TypeError as error:
            raise ValidationError("SRK compounds must be a sequence") from error
        if (
            len(compounds) < 2
            or any(not isinstance(compound, Compound) for compound in compounds)
            or not isinstance(self.interactions, PRInteractions)
            or self.interactions.model != "Soave-Redlich-Kwong"
        ):
            raise ValidationError("SRK thermodynamic-system inputs are invalid")
        compound_ids = tuple(compound.id for compound in compounds)
        if len(set(compound_ids)) != len(compound_ids):
            raise ValidationError("SRK thermodynamic-system compound IDs must be unique")
        for first_index, first in enumerate(compound_ids):
            for second in compound_ids[first_index + 1:]:
                self.interactions.get(first, second)
        object.__setattr__(self, "compounds", compounds)
        object.__setattr__(self, "compound_ids", compound_ids)

    def state(
        self,
        composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
        phase: str,
    ) -> SRKMixtureState:
        return SoaveRedlichKwongMixture(
            self.compounds, composition, self.interactions
        ).state(temperature_k, pressure_pa, phase)

    def stable_state(
        self,
        composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
    ) -> SRKMixtureState:
        return SoaveRedlichKwongMixture(
            self.compounds, composition, self.interactions
        ).stable_state(temperature_k, pressure_pa)

    def tp_flash(
        self,
        composition: tuple[float, ...],
        temperature_k: float,
        pressure_pa: float,
        *,
        max_iterations: int = 100,
        tolerance: float = 1.0e-10,
    ) -> SRKTPFlashResult:
        return srk_tp_flash(
            self.compounds,
            composition,
            self.interactions,
            temperature_k,
            pressure_pa,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )


ThermoSystemConstructor = Callable[..., ThermodynamicSystem]
_THERMO_SYSTEM_CONSTRUCTORS: dict[str, ThermoSystemConstructor] = {
    PENG_ROBINSON_CLASSIC: PengRobinsonSystem,
    NRTL_ACETONE_METHANOL: NRTLSystem,
    IDEAL_RAOULT: IdealRaoultSystem,
    SOAVE_REDLICH_KWONG: SoaveRedlichKwongSystem,
}
THERMO_SYSTEM_CONSTRUCTORS: Mapping[str, ThermoSystemConstructor] = MappingProxyType(
    _THERMO_SYSTEM_CONSTRUCTORS
)


def create_thermo_system(model_id: str, **configuration) -> ThermodynamicSystem:
    """Construct one registered system; runtime plugin registration is unsupported."""
    if not isinstance(model_id, str) or not model_id:
        raise ValidationError("thermodynamic-system model ID must be a non-empty string")
    try:
        constructor = THERMO_SYSTEM_CONSTRUCTORS[model_id]
    except KeyError as error:
        raise ValidationError(f"unsupported thermodynamic-system model ID: {model_id}") from error
    try:
        return constructor(**configuration)
    except TypeError as error:
        raise ValidationError(
            f"invalid configuration for thermodynamic system: {model_id}"
        ) from error
