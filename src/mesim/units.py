from __future__ import annotations

import math
from dataclasses import dataclass
from types import MappingProxyType

from .errors import OutOfRangeError, ValidationError


@dataclass(frozen=True, slots=True)
class _Unit:
    dimension: str
    factor: float
    offset: float = 0.0


def _aliases(dimension: str, factor: float, *names: str, offset: float = 0.0) -> dict[str, _Unit]:
    unit = _Unit(dimension, factor, offset)
    return {name: unit for name in names}


UNITS = MappingProxyType(
    {
        **_aliases("temperature", 1.0, "K"),
        **_aliases("temperature", 1.0, "degC", "°C", "C", offset=273.15),
        **_aliases("temperature", 5.0 / 9.0, "degF", "°F", "F", offset=459.67),
        **_aliases("temperature", 5.0 / 9.0, "degR", "°R", "R"),
        **_aliases("temperature_difference", 1.0, "delta_K", "delta_degC", "K.", "C."),
        **_aliases("temperature_difference", 5.0 / 9.0, "delta_degF", "delta_degR", "F.", "R."),
        **_aliases("pressure", 1.0, "Pa"),
        **_aliases("pressure", 1_000.0, "kPa"),
        **_aliases("pressure", 1_000_000.0, "MPa"),
        **_aliases("pressure", 100_000.0, "bar"),
        **_aliases("pressure", 100.0, "mbar"),
        **_aliases("pressure", 101_325.0, "atm"),
        **_aliases("pressure", 6_894.757293168, "psi", "psia"),
        **_aliases("pressure", 100_000.0, "barg", offset=1.01325),
        **_aliases("pressure", 1_000.0, "kPag", offset=101.325),
        **_aliases("pressure", 1_000_000.0, "MPag", offset=0.101325),
        **_aliases("pressure", 6_894.757293168, "psig", offset=101_325.0 / 6_894.757293168),
        **_aliases("time", 1.0, "s"),
        **_aliases("time", 60.0, "min"),
        **_aliases("time", 3_600.0, "h", "hr"),
        **_aliases("mass", 1.0, "kg"),
        **_aliases("mass", 0.001, "g"),
        **_aliases("mass", 0.45359237, "lb", "lbm"),
        **_aliases("amount", 1.0, "mol"),
        **_aliases("amount", 1_000.0, "kmol"),
        **_aliases("mass_flow", 1.0, "kg/s"),
        **_aliases("mass_flow", 1.0 / 3_600.0, "kg/h", "kg/hr"),
        **_aliases("mass_flow", 0.001, "g/s"),
        **_aliases("mass_flow", 0.45359237, "lb/s"),
        **_aliases("mass_flow", 0.45359237 / 3_600.0, "lb/h", "lb/hr"),
        **_aliases("molar_flow", 1.0, "mol/s"),
        **_aliases("molar_flow", 1.0 / 3_600.0, "mol/h", "mol/hr"),
        **_aliases("molar_flow", 1_000.0, "kmol/s"),
        **_aliases("molar_flow", 1_000.0 / 3_600.0, "kmol/h", "kmol/hr"),
        **_aliases("volumetric_flow", 1.0, "m3/s", "m^3/s"),
        **_aliases("volumetric_flow", 1.0 / 3_600.0, "m3/h", "m3/hr", "m^3/h", "m^3/hr"),
        **_aliases("volumetric_flow", 0.001, "L/s"),
        **_aliases("volumetric_flow", 0.001 / 60.0, "L/min"),
        **_aliases("volumetric_flow", 0.028316846592, "ft3/s", "ft^3/s"),
        **_aliases("volumetric_flow", 0.028316846592 / 60.0, "ft3/min", "ft^3/min"),
        **_aliases("volumetric_flow", 0.003785411784 / 60.0, "gpm"),
        **_aliases("energy", 1.0, "J"),
        **_aliases("energy", 1_000.0, "kJ"),
        **_aliases("energy", 1_000_000.0, "MJ"),
        **_aliases("energy", 4.184, "cal"),
        **_aliases("energy", 4_184.0, "kcal"),
        **_aliases("energy", 1_055.05585262, "Btu"),
        **_aliases("power", 1.0, "W"),
        **_aliases("power", 1_000.0, "kW"),
        **_aliases("power", 1_000_000.0, "MW"),
        **_aliases("enthalpy", 1.0, "J/kg"),
        **_aliases("enthalpy", 1_000.0, "kJ/kg"),
        **_aliases("molar_enthalpy", 1.0, "J/kmol"),
        **_aliases("molar_enthalpy", 1_000.0, "kJ/kmol", "J/mol"),
        **_aliases("molar_enthalpy", 1_000_000.0, "kJ/mol"),
        **_aliases("heat_capacity", 1.0, "J/kg/K"),
        **_aliases("heat_capacity", 1_000.0, "kJ/kg/K"),
        **_aliases("molar_heat_capacity", 1.0, "J/kmol/K"),
        **_aliases("molar_heat_capacity", 1_000.0, "kJ/kmol/K", "J/mol/K"),
        **_aliases("molar_heat_capacity", 1_000_000.0, "kJ/mol/K"),
        **_aliases("density", 1.0, "kg/m3", "kg/m^3"),
        **_aliases("density", 1_000.0, "g/cm3", "g/cm^3"),
        **_aliases("density", 1_000.0, "kg/L"),
        **_aliases("density", 16.01846337, "lb/ft3", "lb/ft^3"),
        **_aliases("viscosity", 1.0, "Pa.s", "Pa·s"),
        **_aliases("viscosity", 0.001, "mPa.s", "mPa·s", "cP"),
        **_aliases("specific_volume", 1.0, "m3/kg", "m^3/kg"),
        **_aliases("specific_volume", 0.028316846592 / 0.45359237, "ft3/lb", "ft^3/lb"),
        **_aliases("kinematic_viscosity", 1.0, "m2/s", "m^2/s"),
        **_aliases("kinematic_viscosity", 0.000001, "cSt"),
        **_aliases("kinematic_viscosity", 0.0001, "cm2/s", "cm^2/s"),
        **_aliases("surface_tension", 1.0, "N/m"),
        **_aliases("surface_tension", 0.001, "dyn/cm"),
        **_aliases("thermal_conductivity", 1.0, "W/m/K", "W/m·K"),
        **_aliases("thermal_conductivity", 0.001, "mW/m/K", "mW/m·K"),
        **_aliases("thermal_conductivity", 1_000.0, "kW/m/K", "kW/m·K"),
        **_aliases("area", 1.0, "m2", "m^2"),
        **_aliases("area", 0.0001, "cm2", "cm^2"),
        **_aliases("area", 0.09290304, "ft2", "ft^2"),
        **_aliases("volume", 1.0, "m3", "m^3"),
        **_aliases("volume", 0.001, "L", "l"),
        **_aliases("volume", 0.000001, "mL", "ml", "cm3", "cm^3"),
        **_aliases("volume", 0.028316846592, "ft3", "ft^3"),
        **_aliases("length", 1.0, "m"),
        **_aliases("length", 0.001, "mm"),
        **_aliases("length", 0.01, "cm"),
        **_aliases("length", 1_000.0, "km"),
        **_aliases("length", 0.0254, "in"),
        **_aliases("length", 0.3048, "ft"),
        **_aliases("velocity", 1.0, "m/s"),
        **_aliases("velocity", 0.3048, "ft/s"),
        **_aliases("velocity", 1_000.0 / 3_600.0, "km/h", "km/hr"),
        **_aliases("pressure_gradient", 1.0, "Pa/m"),
        **_aliases("pressure_gradient", 1_000.0, "kPa/m"),
        **_aliases("pressure_gradient", 100.0, "bar/km"),
        **_aliases("mass_flux", 1.0, "kg/m2/s", "kg/m^2/s"),
        **_aliases("mass_flux", 1.0 / 3_600.0, "kg/m2/h", "kg/m^2/h"),
        **_aliases("mass_flux", 0.45359237 / 0.09290304, "lb/ft2/s", "lb/ft^2/s"),
        **_aliases("heat_flux", 1.0, "W/m2", "W/m^2"),
        **_aliases("heat_flux", 1_000.0, "kW/m2", "kW/m^2"),
        **_aliases("heat_flux", 1_055.05585262 / 3_600.0 / 0.09290304, "Btu/h/ft2", "Btu/h/ft^2"),
        **_aliases("acceleration", 1.0, "m/s2", "m/s^2"),
        **_aliases("acceleration", 0.3048, "ft/s2", "ft/s^2"),
        **_aliases("force", 1.0, "N"),
        **_aliases("force", 1_000.0, "kN"),
        **_aliases("force", 4.4482216152605, "lbf"),
        **_aliases("heat_transfer_coefficient", 1.0, "W/m2/K", "W/m^2/K"),
        **_aliases("heat_transfer_coefficient", 1_000.0, "kW/m2/K", "kW/m^2/K"),
        **_aliases("dimensionless", 1.0, "1", "fraction"),
        **_aliases("dimensionless", 0.01, "%", "percent"),
    }
)

_DIMENSION_ALIASES = {
    "diffusivity": "kinematic_viscosity",
    "entropy": "heat_capacity",
    "molar_entropy": "molar_heat_capacity",
}


def _unit(unit: str) -> _Unit:
    if not isinstance(unit, str) or not unit.strip():
        raise ValidationError("unit must be a non-empty string")
    try:
        return UNITS[unit.strip()]
    except KeyError as exc:
        raise ValidationError(f"unknown unit: {unit!r}") from exc


def _number(value: float) -> float:
    if isinstance(value, bool):
        raise ValidationError("boolean values are not quantities")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"quantity value must be numeric: {value!r}") from exc
    if not math.isfinite(number):
        raise ValidationError("quantity value must be finite")
    return number


def _check_dimension(spec: _Unit, unit: str, dimension: str | None) -> None:
    expected = _DIMENSION_ALIASES.get(dimension, dimension)
    if expected is not None and spec.dimension != expected:
        raise ValidationError(f"unit {unit!r} has dimension {spec.dimension!r}, expected {dimension!r}")


def _check_temperature(si_value: float, dimension: str) -> None:
    if dimension == "temperature" and si_value < 0.0:
        raise OutOfRangeError("temperature cannot be below 0 K")


def unit_dimension(unit: str) -> str:
    return _unit(unit).dimension


def to_si(value: float, unit: str, dimension: str | None = None) -> float:
    spec = _unit(unit)
    _check_dimension(spec, unit, dimension)
    number = _number(value)
    si_value = (number + spec.offset) * spec.factor
    _check_temperature(si_value, spec.dimension)
    return si_value


def from_si(si_value: float, unit: str, dimension: str | None = None) -> float:
    spec = _unit(unit)
    _check_dimension(spec, unit, dimension)
    number = _number(si_value)
    _check_temperature(number, spec.dimension)
    return number / spec.factor - spec.offset


def convert(value: float, from_unit: str, dimension: str, to_unit: str) -> float:
    source = _unit(from_unit)
    target = _unit(to_unit)
    expected = _DIMENSION_ALIASES.get(dimension, dimension)
    if source.dimension != expected or target.dimension != expected:
        raise ValidationError(f"cannot convert {from_unit!r} to {to_unit!r} as {dimension!r}")
    return from_si(to_si(value, from_unit, dimension), to_unit, dimension)


@dataclass(frozen=True, slots=True)
class Quantity:
    value: float
    unit: str
    si_value: float

    @classmethod
    def from_value(cls, value: float, unit: str, dimension: str | None = None) -> "Quantity":
        return cls(value, unit, to_si(value, unit, dimension))
