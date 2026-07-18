"""Scoped DWSIM Seawater equilibrium behavior for Water/Salt solutions."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import ValidationError


@dataclass(frozen=True, slots=True)
class SeawaterData:
    source_revision: str
    case_sha256: str
    runtime_assembly_sha256: str
    property_package_source_sha256: str
    model_source_sha256: str
    compound_ids: tuple[str, ...]
    molecular_weights: tuple[float, ...]
    salinity_limit: float
    vapor_pressure_coefficients: tuple[float, ...]
    probe_temperature_k: float
    probe_composition: tuple[float, ...]
    probe_mass_fractions: tuple[float, ...]
    probe_salinity: float
    current_stream_salinity: float
    probe_vapor_pressure_pa: float


@dataclass(frozen=True, slots=True)
class SeawaterTPFlashResult:
    liquid_fraction: float
    vapor_fraction: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    equilibrium_ratios: tuple[float, ...]
    iterations: int


def load_seawater_data(path: str | Path) -> SeawaterData:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        source = document["source"]
        probe = document["scoped_probe"]
        if (
            document["schema_version"] != "dwsim-seawater-data-1"
            or document["model"] != "Seawater"
        ):
            raise ValidationError("unsupported Seawater data schema or model")
        data = SeawaterData(
            source["revision"],
            source["case_sha256"],
            source["runtime_assembly_sha256"],
            source["property_package_source_sha256"],
            source["model_source_sha256"],
            tuple(document["compound_ids"]),
            tuple(document["molecular_weights"]),
            document["salinity_limit"],
            tuple(document["vapor_pressure_coefficients"]),
            probe["temperature_k"],
            tuple(probe["composition"]),
            tuple(probe["converted_mass_fractions"]),
            probe["calculated_salinity"],
            probe["current_stream_salinity"],
            probe["vapor_pressure_pa"],
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid Seawater data: {error}") from error
    hashes = (
        data.case_sha256,
        data.runtime_assembly_sha256,
        data.property_package_source_sha256,
        data.model_source_sha256,
    )
    numeric = (
        *data.molecular_weights,
        data.salinity_limit,
        *data.vapor_pressure_coefficients,
        data.probe_temperature_k,
        *data.probe_composition,
        *data.probe_mass_fractions,
        data.probe_salinity,
        data.current_stream_salinity,
        data.probe_vapor_pressure_pa,
    )
    if (
        not data.source_revision
        or data.compound_ids != ("Water", "Salt")
        or len(data.molecular_weights) != 2
        or len(data.vapor_pressure_coefficients) != 6
        or len(data.probe_composition) != 2
        or len(data.probe_mass_fractions) != 2
        or any(len(value) != 64 for value in hashes)
        or any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in numeric)
        or any(value <= 0.0 for value in data.molecular_weights)
        or not 0.0 < data.salinity_limit < 1.0
        or not 0.0 <= data.probe_salinity <= data.salinity_limit
        or not 0.0 <= data.current_stream_salinity <= data.salinity_limit
    ):
        raise ValidationError("Seawater source identity or scoped data are invalid")
    return data


def _composition(data: SeawaterData, composition) -> tuple[float, float]:
    if not isinstance(data, SeawaterData):
        raise ValidationError("invalid Seawater data")
    try:
        values = tuple(composition)
    except TypeError as error:
        raise ValidationError("Seawater composition must be a sequence") from error
    if (
        len(values) != 2
        or any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0.0 for value in values)
        or not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        or values[0] <= 0.0
    ):
        raise ValidationError("Seawater requires normalized Water/Salt mole fractions")
    return float(values[0]), float(values[1])


def seawater_salinity(data: SeawaterData, composition) -> float:
    """Return DWSIM's kg Salt/kg Water value, including its 0.12 ceiling."""
    water, salt = _composition(data, composition)
    water_mass = water * data.molecular_weights[0]
    salt_mass = salt * data.molecular_weights[1]
    salinity = salt_mass / water_mass
    return min(salinity, data.salinity_limit)


def seawater_vapor_pressure_pa(
    data: SeawaterData,
    temperature_k: float,
    *,
    composition=None,
    salinity: float | None = None,
) -> float:
    if isinstance(temperature_k, bool) or not isinstance(temperature_k, (int, float)) or not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValidationError("Seawater temperature must be positive and finite")
    if (composition is None) == (salinity is None):
        raise ValidationError("provide exactly one Seawater composition or salinity")
    if composition is not None:
        salinity_value = seawater_salinity(data, composition)
    else:
        if isinstance(salinity, bool) or not isinstance(salinity, (int, float)) or not math.isfinite(salinity) or not 0.0 <= salinity <= data.salinity_limit:
            raise ValidationError("Seawater salinity is outside the frozen source range")
        salinity_value = float(salinity)
    a1, a2, a3, a4, a5, a6 = data.vapor_pressure_coefficients
    try:
        pure_water = math.exp(
            a1 / temperature_k
            + a2
            + a3 * temperature_k
            + a4 * temperature_k**2
            + a5 * temperature_k**3
            + a6 * math.log(temperature_k)
        )
        salinity_g_kg = salinity_value * 1000.0
        value = pure_water / (
            1.0 + 0.57357 * salinity_g_kg / (1000.0 - salinity_g_kg)
        )
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError("Seawater vapor pressure is unrepresentable") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValidationError("Seawater vapor pressure is unrepresentable")
    return value


def seawater_fugacity_coefficients(
    data: SeawaterData,
    composition,
    temperature_k: float,
    pressure_pa: float,
    phase: str,
) -> tuple[float, float]:
    _composition(data, composition)
    if phase not in {"liquid", "vapor"} or isinstance(pressure_pa, bool) or not isinstance(pressure_pa, (int, float)) or not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise ValidationError("invalid Seawater phase or pressure")
    if phase == "vapor":
        return (1.0, 1.0)
    # DWSIM's indexed AUX_PVAPi path ignores Vx and reads CurrentMaterialStream.
    water_vapor_pressure = seawater_vapor_pressure_pa(
        data, temperature_k, salinity=data.current_stream_salinity
    )
    return (water_vapor_pressure / pressure_pa, 0.0)


def seawater_tp_flash(
    data: SeawaterData,
    composition,
    temperature_k: float,
    pressure_pa: float,
) -> SeawaterTPFlashResult:
    values = _composition(data, composition)
    seawater_fugacity_coefficients(data, values, temperature_k, pressure_pa, "liquid")
    if pressure_pa < seawater_vapor_pressure_pa(data, temperature_k, composition=values):
        raise ValidationError("Seawater TP flash is scoped to the captured all-liquid domain")
    return SeawaterTPFlashResult(
        1.0,
        0.0,
        values,
        (1.0e-10, 1.0e-10),
        (1.0e-10, 1.0e-10),
        1,
    )


__all__ = (
    "SeawaterData",
    "SeawaterTPFlashResult",
    "load_seawater_data",
    "seawater_salinity",
    "seawater_vapor_pressure_pa",
    "seawater_fugacity_coefficients",
    "seawater_tp_flash",
)
