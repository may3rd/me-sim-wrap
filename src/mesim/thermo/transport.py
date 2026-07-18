import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from ..compounds import Compound
from ..errors import OutOfRangeError, ValidationError
from .correlations import evaluate_temperature_equation


# DWSIM PengRobinson.vb AUX_Ci fallback values for the catalog compounds.
_PENELOUX = {"Methane": -0.1595, "Ethane": -0.1134, "Propane": -0.0863, "N-butane": -0.0675, "N-pentane": -0.039}
R = 8314.46261815324  # J/kmol/K


@dataclass(frozen=True)
class TransportCorrelation:
    equation: int
    coefficients: tuple[float, float, float, float, float]
    minimum_k: float
    maximum_k: float
    unit: str

    def value(self, temperature_k: float, allow_extrapolation: bool = False) -> float:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        if not self.minimum_k <= temperature_k <= self.maximum_k:
            if not allow_extrapolation:
                raise OutOfRangeError(f"transport correlation extrapolated outside {self.minimum_k:g}..{self.maximum_k:g} K")
        value = evaluate_temperature_equation(
            self.equation, self.coefficients, temperature_k
        )
        if not math.isfinite(value) or value <= 0:
            raise ValidationError("transport correlation produced a non-positive value")
        return value


@dataclass(frozen=True)
class LiquidTransportCorrelation:
    equation: int
    coefficients: tuple[float, float, float, float, float]
    minimum_k: float | None
    maximum_k: float | None
    unit: str

    def value(self, temperature_k: float, allow_extrapolation: bool = False) -> float:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        if (
            self.minimum_k is not None and self.maximum_k is not None
            and not self.minimum_k <= temperature_k <= self.maximum_k
            and not allow_extrapolation
        ):
            raise OutOfRangeError(
                f"liquid transport correlation extrapolated outside {self.minimum_k:g}..{self.maximum_k:g} K",
            )
        value = evaluate_temperature_equation(
            self.equation, self.coefficients, temperature_k
        )
        if not math.isfinite(value) or value <= 0:
            raise ValidationError("liquid transport correlation produced a non-positive value")
        return value


@dataclass(frozen=True)
class TransportRecord:
    compound_id: str
    critical_volume_m3_per_kmol: float
    vapor_viscosity: TransportCorrelation
    vapor_thermal_conductivity: TransportCorrelation
    liquid_viscosity: LiquidTransportCorrelation
    liquid_thermal_conductivity: LiquidTransportCorrelation


@dataclass(frozen=True)
class VaporTransport:
    dynamic_viscosity_pa_s: float
    thermal_conductivity_w_per_m_k: float


@dataclass(frozen=True)
class LiquidTransport:
    dynamic_viscosity_pa_s: float
    thermal_conductivity_w_per_m_k: float


def load_transport_correlations(path: str | Path) -> tuple[TransportRecord, ...]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if data["schema_version"] != "transport-correlations-3":
            raise ValidationError("unsupported transport correlation schema")
        provenance = data["provenance"]
        imported = datetime.fromisoformat(provenance["imported_utc"])
        if not all(isinstance(value, str) and value for value in provenance.values()) or not provenance["imported_utc"].endswith("Z") or imported.utcoffset() != timedelta(0):
            raise ValidationError("transport correlation provenance must be non-empty and use a UTC import timestamp")
        records = tuple(
            TransportRecord(
                record["compound_id"],
                _positive(record["critical_volume"], "m3/kmol", "critical volume"),
                _correlation(record["vapor_viscosity"], (2, 3, 16, 102), "Pa.s"),
                _correlation(record["vapor_thermal_conductivity"], (3, 16, 102), "W/m/K"),
                _liquid_correlation(record["liquid_viscosity"], (10, 16, 101), "Pa.s"),
                _liquid_correlation(record["liquid_thermal_conductivity"], (3, 16, 100), "W/m/K"),
            )
            for record in data["correlations"]
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid transport correlation data: {error}") from error
    if not records or len({record.compound_id for record in records}) != len(records) or any(not isinstance(record.compound_id, str) or not record.compound_id for record in records):
        raise ValidationError("transport correlation compound IDs must be non-empty and unique")
    return records


def _positive(data: dict, unit: str, label: str) -> float:
    value = data["value"]
    if data["unit"] != unit or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValidationError(f"{label} must be a positive finite {unit} value")
    return value


def _correlation(data: dict, equations: tuple[int, ...], unit: str) -> TransportCorrelation:
    values = tuple(data[key] for key in "ABCDE")
    minimum, maximum = data["minimum_k"], data["maximum_k"]
    if data["equation"] not in equations or data["unit"] != unit or any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in (*values, minimum, maximum)) or minimum <= 0 or maximum <= minimum:
        raise ValidationError(f"invalid transport correlation in {unit}")
    return TransportCorrelation(data["equation"], values, minimum, maximum, unit)


def _liquid_correlation(
    data: dict, equations: tuple[int, ...], unit: str,
) -> LiquidTransportCorrelation:
    values = tuple(data[key] for key in "ABCDE")
    minimum = data.get("minimum_k")
    maximum = data.get("maximum_k")
    if (
        data["equation"] not in equations or data["unit"] != unit
        or any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value in values)
        or (
            isinstance(minimum, bool) or isinstance(maximum, bool)
            or not isinstance(minimum, (int, float)) or not isinstance(maximum, (int, float))
            or not math.isfinite(minimum) or not math.isfinite(maximum)
            or minimum <= 0 or maximum <= minimum
        )
    ):
        raise ValidationError(f"invalid liquid transport correlation in {unit}")
    return LiquidTransportCorrelation(data["equation"], values, minimum, maximum, unit)


def vapor_transport(compounds: tuple[Compound, ...], mole_fractions: tuple[float, ...], records: tuple[TransportRecord, ...], temperature_k: float, density_kg_per_m3: float) -> VaporTransport:
    if len(compounds) != len(mole_fractions) or len(compounds) != len(records) or not compounds:
        raise ValidationError("transport components, mole fractions, and records must have equal non-zero lengths")
    if any(compound.id != record.compound_id for compound, record in zip(compounds, records)):
        raise ValidationError("transport records must match component order")
    if not all(math.isfinite(value) and value >= 0 for value in mole_fractions) or not math.isclose(sum(mole_fractions), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("transport mole fractions must be finite, non-negative, and sum to one")
    if not math.isfinite(density_kg_per_m3) or density_kg_per_m3 <= 0:
        raise ValidationError("vapor density must be finite and positive")
    viscosity = sum(fraction * record.vapor_viscosity.value(temperature_k) for fraction, record in zip(mole_fractions, records))
    conductivity = sum(fraction * record.vapor_thermal_conductivity.value(temperature_k) for fraction, record in zip(mole_fractions, records))
    molecular_weight = sum(fraction * compound.molecular_weight.value for fraction, compound in zip(mole_fractions, compounds))
    critical_temperature = sum(fraction * compound.critical_temperature.value for fraction, compound in zip(mole_fractions, compounds))
    critical_pressure = sum(fraction * compound.critical_pressure.value for fraction, compound in zip(mole_fractions, compounds))
    critical_volume = sum(fraction * record.critical_volume_m3_per_kmol for fraction, record in zip(mole_fractions, records)) / 1000.0
    reduced_density = critical_volume / (molecular_weight / density_kg_per_m3 / 1000.0)
    xi = (critical_temperature / molecular_weight**3 / (critical_pressure / 101325.0) ** 4) ** (1.0 / 6.0)
    correction = (1.023 + 0.23364 * reduced_density + 0.58533 * reduced_density**2 - 0.40758 * reduced_density**3 + 0.093324 * reduced_density**4) ** 4
    viscosity = ((correction - 1.0) / xi + viscosity * 10_000_000.0) / 10_000_000.0
    if not all(math.isfinite(value) and value > 0 for value in (viscosity, conductivity)):
        raise ValidationError("vapor transport is outside the representable range")
    return VaporTransport(viscosity, conductivity)


def liquid_transport(
    compounds: tuple[Compound, ...], mole_fractions: tuple[float, ...],
    records: tuple[TransportRecord, ...], temperature_k: float,
    allow_extrapolation: bool = False,
) -> LiquidTransport:
    """Calculate DWSIM's mole-average viscosity and Li liquid conductivity."""
    if len(compounds) != len(mole_fractions) or len(compounds) != len(records) or not compounds:
        raise ValidationError("liquid transport components, mole fractions, and records must have equal non-zero lengths")
    if any(compound.id != record.compound_id for compound, record in zip(compounds, records)):
        raise ValidationError("liquid transport records must match component order")
    if not isinstance(allow_extrapolation, bool):
        raise ValidationError("liquid transport extrapolation flag must be boolean")
    if (
        not all(math.isfinite(value) and value >= 0 for value in mole_fractions)
        or not math.isclose(sum(mole_fractions), 1.0, rel_tol=0.0, abs_tol=1e-12)
    ):
        raise ValidationError("liquid transport fractions must be finite, non-negative, and sum to one")
    viscosities = tuple(
        record.liquid_viscosity.value(temperature_k, allow_extrapolation)
        for record in records
    )
    conductivities = tuple(
        record.liquid_thermal_conductivity.value(temperature_k, allow_extrapolation)
        for record in records
    )
    viscosity = math.fsum(
        fraction * value for fraction, value in zip(mole_fractions, viscosities)
    )
    volume_weights = tuple(
        fraction * record.critical_volume_m3_per_kmol
        for fraction, record in zip(mole_fractions, records)
    )
    weight_total = math.fsum(volume_weights)
    if not math.isfinite(weight_total) or weight_total <= 0:
        raise ValidationError("liquid transport critical-volume weighting is non-positive")
    volume_fractions = tuple(weight / weight_total for weight in volume_weights)
    conductivity = math.fsum(
        left_fraction * right_fraction * 2.0 / (1.0 / left_value + 1.0 / right_value)
        for left_fraction, left_value in zip(volume_fractions, conductivities)
        for right_fraction, right_value in zip(volume_fractions, conductivities)
    )
    if not all(math.isfinite(value) and value > 0 for value in (viscosity, conductivity)):
        raise ValidationError("liquid transport is outside the representable range")
    return LiquidTransport(viscosity, conductivity)


def translated_vapor_density(compounds: tuple[Compound, ...], mole_fractions: tuple[float, ...], temperature_k: float, pressure_pa: float, compressibility: float) -> float:
    if len(compounds) != len(mole_fractions) or not compounds or not all(math.isfinite(value) and value >= 0 for value in mole_fractions) or not math.isclose(sum(mole_fractions), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("translated-density components and mole fractions must be non-empty, finite, and sum to one")
    if not all(math.isfinite(value) and value > 0 for value in (temperature_k, pressure_pa, compressibility)):
        raise ValidationError("temperature, pressure, and compressibility must be finite and positive")
    try:
        translation = sum(
            fraction * _PENELOUX[compound.id] * 0.07780 * R * compound.critical_temperature.value / compound.critical_pressure.value
            for compound, fraction in zip(compounds, mole_fractions)
        )
    except KeyError as error:
        raise ValidationError(f"missing Peneloux coefficient: {error.args[0]}") from error
    molecular_weight = sum(fraction * compound.molecular_weight.value for compound, fraction in zip(compounds, mole_fractions))
    volume = R * compressibility * temperature_k / pressure_pa - translation
    if not math.isfinite(volume) or volume <= 0:
        raise ValidationError("translated vapor volume must be finite and positive")
    return molecular_weight / volume
