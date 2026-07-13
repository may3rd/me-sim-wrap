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
        **_aliases("pressure", 1.0, "Pa"),
        **_aliases("pressure", 1_000.0, "kPa"),
        **_aliases("pressure", 1_000_000.0, "MPa"),
        **_aliases("pressure", 100_000.0, "bar"),
        **_aliases("pressure", 100.0, "mbar"),
        **_aliases("pressure", 101_325.0, "atm"),
        **_aliases("pressure", 6_894.757293168, "psi", "psia"),
        **_aliases("mass_flow", 1.0, "kg/s"),
        **_aliases("mass_flow", 1.0 / 3_600.0, "kg/h", "kg/hr"),
        **_aliases("mass_flow", 0.001, "g/s"),
        **_aliases("mass_flow", 0.45359237, "lb/s"),
        **_aliases("mass_flow", 0.45359237 / 3_600.0, "lb/h", "lb/hr"),
        **_aliases("molar_flow", 1.0, "mol/s"),
        **_aliases("molar_flow", 1.0 / 3_600.0, "mol/h", "mol/hr"),
        **_aliases("molar_flow", 1_000.0, "kmol/s"),
        **_aliases("molar_flow", 1_000.0 / 3_600.0, "kmol/h", "kmol/hr"),
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
        **_aliases("molar_enthalpy", 1.0, "J/mol"),
        **_aliases("molar_enthalpy", 1_000.0, "kJ/mol"),
        **_aliases("molar_enthalpy", 1_000_000.0, "kJ/kmol"),
        **_aliases("density", 1.0, "kg/m3", "kg/m^3"),
        **_aliases("density", 1_000.0, "g/cm3", "g/cm^3"),
        **_aliases("density", 1_000.0, "kg/L"),
        **_aliases("density", 16.01846337, "lb/ft3", "lb/ft^3"),
        **_aliases("viscosity", 1.0, "Pa.s", "Pa·s"),
        **_aliases("viscosity", 0.001, "mPa.s", "mPa·s", "cP"),
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
        **_aliases("heat_transfer_coefficient", 1.0, "W/m2/K", "W/m^2/K"),
        **_aliases("heat_transfer_coefficient", 1_000.0, "kW/m2/K", "kW/m^2/K"),
        **_aliases("dimensionless", 1.0, "1", "fraction"),
        **_aliases("dimensionless", 0.01, "%", "percent"),
    }
)


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
    if dimension is not None and spec.dimension != dimension:
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
    if source.dimension != dimension or target.dimension != dimension:
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
