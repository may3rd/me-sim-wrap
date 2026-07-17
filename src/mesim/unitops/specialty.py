"""Source-backed specialty clean-energy unit operations.

These functions reproduce the compact DWSIM renewable-energy equations.  Air
density and liquid density remain explicit thermodynamic inputs so this module
does not hide a weather service, humid-air package, or manufacturer curve.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..errors import ValidationError


DWSIM_GRAVITY_M_S2 = 9.8
DWSIM_WIND_POWER_COEFFICIENT = 8.0 / 27.0


def _finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValidationError(f"{name} must be finite")
    return float(value)


def _positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise ValidationError(f"{name} must be positive")
    return result


def _efficiency(value: float) -> float:
    result = _finite(value, "efficiency")
    if not 0.0 <= result <= 100.0:
        raise ValidationError("efficiency must be between zero and 100 percent")
    return result / 100.0


def _unit_count(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValidationError(f"{name} must be a positive integer")
    return value


@dataclass(frozen=True, slots=True)
class SolarPanelResult:
    generated_power_kw: float


def solar_panel_power(
    solar_irradiation_kw_m2: float,
    panel_area_m2: float,
    efficiency_percent: float,
    number_of_panels: int,
) -> SolarPanelResult:
    """Calculate DWSIM's area-times-irradiation photovoltaic power."""
    irradiation = _finite(solar_irradiation_kw_m2, "solar irradiation")
    if irradiation < 0.0:
        raise ValidationError("solar irradiation must be non-negative")
    area = _positive(panel_area_m2, "panel area")
    efficiency = _efficiency(efficiency_percent)
    count = _unit_count(number_of_panels, "number of panels")
    return SolarPanelResult(irradiation * area * count * efficiency)


@dataclass(frozen=True, slots=True)
class WindTurbineResult:
    rotor_diameter_m: float
    maximum_theoretical_power_kw: float
    generated_power_kw: float


def wind_turbine_power(
    air_density_kg_m3: float,
    wind_speed_m_s: float,
    disk_area_m2: float,
    efficiency_percent: float,
    number_of_turbines: int,
) -> WindTurbineResult:
    """Calculate DWSIM's Betz-limit wind power from an explicit air density."""
    density = _positive(air_density_kg_m3, "air density")
    wind_speed = _finite(wind_speed_m_s, "wind speed")
    if wind_speed < 0.0:
        raise ValidationError("wind speed must be non-negative")
    area = _positive(disk_area_m2, "disk area")
    efficiency = _efficiency(efficiency_percent)
    count = _unit_count(number_of_turbines, "number of turbines")
    rotor_diameter = math.sqrt(4.0 * area / math.pi)
    maximum_power = (
        count
        * DWSIM_WIND_POWER_COEFFICIENT
        * density
        * wind_speed ** 3
        * area
        / 1000.0
    )
    return WindTurbineResult(rotor_diameter, maximum_power, maximum_power * efficiency)


@dataclass(frozen=True, slots=True)
class HydroelectricTurbineResult:
    volumetric_flow_m3_s: float
    velocity_head_m: float
    total_head_m: float
    generated_power_kw: float
    outlet_specific_enthalpy_change_kj_kg: float
    energy_balance_residual_kw: float


def hydroelectric_turbine_power(
    mass_flow_kg_s: float,
    liquid_density_kg_m3: float,
    static_head_m: float,
    inlet_velocity_m_s: float,
    outlet_velocity_m_s: float,
    efficiency_percent: float,
    *,
    gravity_m_s2: float = DWSIM_GRAVITY_M_S2,
) -> HydroelectricTurbineResult:
    """Calculate DWSIM's hydrostatic plus velocity-head turbine power."""
    mass_flow = _positive(mass_flow_kg_s, "mass flow")
    density = _positive(liquid_density_kg_m3, "liquid density")
    static_head = _finite(static_head_m, "static head")
    inlet_velocity = _finite(inlet_velocity_m_s, "inlet velocity")
    outlet_velocity = _finite(outlet_velocity_m_s, "outlet velocity")
    if inlet_velocity < 0.0 or outlet_velocity < 0.0:
        raise ValidationError("hydroelectric velocities must be non-negative")
    gravity = _positive(gravity_m_s2, "gravity")
    efficiency = _efficiency(efficiency_percent)
    volumetric_flow = mass_flow / density
    velocity_head = (inlet_velocity ** 2 - outlet_velocity ** 2) / (2.0 * gravity)
    total_head = static_head + velocity_head
    if total_head <= 0.0:
        raise ValidationError("hydroelectric turbine total head must be positive")
    generated_power = efficiency * density * gravity * total_head * volumetric_flow / 1000.0
    enthalpy_change = -generated_power / mass_flow
    energy_residual = generated_power + mass_flow * enthalpy_change
    return HydroelectricTurbineResult(
        volumetric_flow,
        velocity_head,
        total_head,
        generated_power,
        enthalpy_change,
        energy_residual,
    )
