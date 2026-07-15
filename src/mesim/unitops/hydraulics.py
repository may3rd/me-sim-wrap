"""Single-phase hydraulic calculations."""
import math
from dataclasses import dataclass

from ..errors import ValidationError


def _dwsim_friction_factor(reynolds: float, diameter_m: float, roughness_m: float) -> float:
    if reynolds > 4000.0:
        a1 = math.log10(((roughness_m / diameter_m) ** 1.1096) / 2.8257 + (5.8506 / reynolds) ** 0.8961)
        b1 = -2.0 * math.log10(roughness_m / diameter_m / 3.7065 - 5.0452 * a1 / reynolds)
        return (1.0 / b1) ** 2
    if reynolds < 2100.0:
        return 64.0 / reynolds
    a = (8.0 / reynolds) ** 12
    b = (2.457 * math.log(1.0 / ((7.0 / reynolds) ** 0.9 + 0.27 * roughness_m / diameter_m))) ** 16
    c = (37530.0 / reynolds) ** 16
    return 8.0 * (a + 1.0 / (b + c) ** 1.5) ** (1.0 / 12.0)


@dataclass(frozen=True, slots=True)
class PipePressureDrop:
    velocity_m_s: float
    reynolds: float
    friction_factor: float
    friction_drop_pa: float
    static_drop_pa: float
    total_drop_pa: float


@dataclass(frozen=True, slots=True)
class OrificePressureDrop:
    reynolds: float
    discharge_coefficient: float
    orifice_drop_pa: float
    overall_drop_pa: float


@dataclass(frozen=True, slots=True)
class LockhartMartinelliPressureDrop:
    vapor_reynolds: float
    liquid_reynolds: float
    martinelli_parameter: float
    liquid_holdup: float
    friction_drop_pa: float
    static_drop_pa: float
    total_drop_pa: float


@dataclass(frozen=True, slots=True)
class BeggsBrillPressureDrop:
    flow_regime: str
    liquid_holdup: float
    friction_drop_pa: float
    static_drop_pa: float
    total_drop_pa: float


def beggs_brill_pressure_drop(
    diameter_m: float, length_m: float, elevation_m: float, roughness_m: float,
    vapor_flow_m3_s: float, liquid_flow_m3_s: float, vapor_density_kg_m3: float,
    liquid_density_kg_m3: float, vapor_viscosity_pa_s: float, liquid_viscosity_pa_s: float,
    surface_tension_n_m: float,
) -> BeggsBrillPressureDrop:
    """DWSIM Beggs-Brill two-phase pressure drop with supplied phase properties."""
    values = ((diameter_m, "pipe diameter"), (length_m, "pipe length"), (roughness_m, "pipe roughness"), (vapor_flow_m3_s, "vapor flow"), (liquid_flow_m3_s, "liquid flow"), (vapor_density_kg_m3, "vapor density"), (liquid_density_kg_m3, "liquid density"), (vapor_viscosity_pa_s, "vapor viscosity"), (liquid_viscosity_pa_s, "liquid viscosity"), (surface_tension_n_m, "surface tension"))
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0 for value, _ in values):
        raise ValidationError("Beggs-Brill inputs must be finite and positive")
    if not isinstance(elevation_m, (int, float)) or isinstance(elevation_m, bool) or not math.isfinite(elevation_m) or abs(elevation_m) >= length_m:
        raise ValidationError("Beggs-Brill elevation magnitude must be below pipe length")
    diameter_ft, length_ft, elevation_ft = diameter_m * 3.28084, length_m * 3.28084, elevation_m * 3.28084
    vapor_flow_ft3_s, liquid_flow_ft3_s = vapor_flow_m3_s * 35.314668, liquid_flow_m3_s * 35.314668
    vapor_density, liquid_density = vapor_density_kg_m3 * 0.062428, liquid_density_kg_m3 * 0.062428
    area = math.pi * diameter_ft**2 / 4.0
    velocity = (vapor_flow_ft3_s + liquid_flow_ft3_s) / area
    liquid_fraction = liquid_flow_ft3_s / (vapor_flow_ft3_s + liquid_flow_ft3_s)
    froude = velocity**2 / (32.2 * diameter_ft)
    limit1, limit2, limit3, limit4 = 316.0 * liquid_fraction**0.302, 0.0009252 * liquid_fraction**-2.4684, 0.1 * liquid_fraction**-1.4516, 0.5 * liquid_fraction**-6.738
    if (liquid_fraction < 0.01 and froude < limit1) or (liquid_fraction >= 0.01 and froude < limit2):
        regime = "Segregated"
    elif (liquid_fraction >= 0.01 and liquid_fraction < 0.4 and limit3 < froude <= limit1) or (liquid_fraction >= 0.4 and limit3 < froude <= limit4):
        regime = "Intermittent"
    elif (liquid_fraction < 0.4 and froude >= limit1) or (liquid_fraction >= 0.4 and froude > limit4):
        regime = "Distributed"
    elif liquid_fraction >= 0.01 and limit2 < froude < limit3:
        regime = "Transition"
    else:
        raise ValidationError("Beggs-Brill state is outside the DWSIM flow-pattern map")
    horizontal = {"Segregated": lambda: 0.98 * liquid_fraction**0.4846 / froude**0.0868, "Intermittent": lambda: 0.845 * liquid_fraction**0.5351 / froude**0.0173, "Distributed": lambda: 1.065 * liquid_fraction**0.5824 / froude**0.0609}
    if regime == "Transition":
        fraction = (limit3 - froude) / (limit3 - limit2)
        holdup_0 = fraction * horizontal["Segregated"]() + (1.0 - fraction) * horizontal["Intermittent"]()
    else:
        holdup_0 = horizontal[regime]()
    liquid_velocity = liquid_flow_ft3_s / area
    liquid_velocity_number = 1.938 * liquid_velocity * (liquid_density / (32.2 * surface_tension_n_m * 1000.0)) ** 0.25
    if elevation_ft > 0.0 and regime == "Segregated":
        beta = (1.0 - liquid_fraction) * math.log(0.011 * liquid_velocity_number**3.539 / (liquid_fraction**3.768 * froude**1.614))
    elif elevation_ft > 0.0 and regime == "Intermittent":
        beta = (1.0 - liquid_fraction) * math.log(2.96 * liquid_fraction**0.305 * froude**0.0978 / liquid_velocity_number**0.4473)
    elif elevation_ft < 0.0:
        beta = (1.0 - liquid_fraction) * math.log(4.7 * liquid_velocity_number**0.1244 / (liquid_fraction**0.3692 * froude**0.5056))
    else:
        beta = 0.0
    angle = math.atan(elevation_ft / (length_ft**2 - elevation_ft**2) ** 0.5)
    holdup = holdup_0 * (1.0 + max(beta, 0.0) * (math.sin(1.8 * angle) - 0.3333 * math.sin(1.8 * angle) ** 3))
    mixture_density = holdup * liquid_density + (1.0 - holdup) * vapor_density
    source_reynolds = (liquid_fraction * liquid_density + (1.0 - liquid_fraction) * vapor_density) * velocity * diameter_ft / ((liquid_fraction * liquid_viscosity_pa_s * 1000.0 + (1.0 - liquid_fraction) * vapor_viscosity_pa_s * 1000.0) * 0.00067197)
    friction = _dwsim_friction_factor(source_reynolds, diameter_m, roughness_m)
    y = math.log(liquid_fraction / holdup**2)
    correction = math.log(2.2 * math.exp(y) - 1.2) if 1.0 < y < 1.2 else y / (-0.0523 + 3.182 * y - 0.8725 * y**2 + 0.01853 * y**4)
    friction_drop = friction * math.exp(correction) * velocity**2 / 2.0 * (liquid_fraction * liquid_density + (1.0 - liquid_fraction) * vapor_density) * length_ft / (32.2 * diameter_ft) * 47.88
    static_drop = mixture_density * elevation_ft * 47.88
    return BeggsBrillPressureDrop(regime, holdup, friction_drop, static_drop, friction_drop + static_drop)


def lockhart_martinelli_pressure_drop(
    diameter_m: float, length_m: float, elevation_m: float, roughness_m: float,
    vapor_flow_m3_s: float, liquid_flow_m3_s: float, vapor_density_kg_m3: float,
    liquid_density_kg_m3: float, vapor_viscosity_pa_s: float, liquid_viscosity_pa_s: float,
) -> LockhartMartinelliPressureDrop:
    """DWSIM homogeneous Lockhart-Martinelli two-phase pressure drop."""
    vapor = pipe_pressure_drop(diameter_m, length_m, 0.0, roughness_m, vapor_flow_m3_s, vapor_density_kg_m3, vapor_viscosity_pa_s)
    liquid = pipe_pressure_drop(diameter_m, length_m, 0.0, roughness_m, liquid_flow_m3_s, liquid_density_kg_m3, liquid_viscosity_pa_s)
    if vapor_flow_m3_s <= 0.0 or liquid_flow_m3_s <= 0.0:
        raise ValidationError("Lockhart-Martinelli requires positive vapor and liquid flows")
    x = (liquid.friction_drop_pa / vapor.friction_drop_pa) ** 0.5
    liquid_multiplier = 1.0 + 20.0 / x + 1.0 / x**2
    vapor_multiplier = 1.0 + 20.0 * x + x**2
    vapor_fraction = vapor_flow_m3_s / (vapor_flow_m3_s + liquid_flow_m3_s)
    static_drop = (vapor_fraction * vapor_density_kg_m3 + (1.0 - vapor_fraction) * liquid_density_kg_m3) * 9.8 * elevation_m
    friction_drop = max(liquid_multiplier * liquid.friction_drop_pa, vapor_multiplier * vapor.friction_drop_pa)
    return LockhartMartinelliPressureDrop(vapor.reynolds, liquid.reynolds, x, (1.0 / liquid_multiplier) ** 0.5, friction_drop, static_drop, friction_drop + static_drop)


def orifice_pressure_drop(
    pipe_diameter_m: float, orifice_diameter_m: float, mass_flow_kg_s: float,
    density_kg_m3: float, viscosity_pa_s: float, tap: str, correction_factor: float = 1.0,
) -> OrificePressureDrop:
    """DWSIM ISO-5167-style incompressible orifice pressure drop."""
    values = ((pipe_diameter_m, "pipe diameter"), (orifice_diameter_m, "orifice diameter"), (mass_flow_kg_s, "mass flow"), (density_kg_m3, "density"), (viscosity_pa_s, "viscosity"), (correction_factor, "correction factor"))
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0 for value, _ in values):
        raise ValidationError("orifice diameters, flow, properties, and correction factor must be finite and positive")
    if orifice_diameter_m >= pipe_diameter_m:
        raise ValidationError("orifice diameter must be below pipe diameter")
    if tap == "corner":
        separation, l1, l2 = 0.0, 0.0, 0.0
    elif tap == "flange":
        separation, l1, l2 = 0.0508, 0.0254 / orifice_diameter_m, 0.0254 / orifice_diameter_m
    elif tap == "radius":
        separation, l1, l2 = 1.5 * orifice_diameter_m, 1.0, 0.47
    else:
        raise ValidationError("orifice tap must be corner, flange, or radius")
    beta = orifice_diameter_m / pipe_diameter_m
    pipe_area = 3.1416 * pipe_diameter_m**2 / 4.0
    orifice_area = 3.1416 * orifice_diameter_m**2 / 4.0
    reynolds = mass_flow_kg_s * pipe_diameter_m / (pipe_area * viscosity_pa_s)
    a = (19000.0 * beta / reynolds) ** 0.8
    m2 = 2.0 * l2 / (1.0 - beta)
    coefficient = 0.5961 + 0.0261 * beta**2 - 0.216 * beta**8 + 0.000521 * (1e6 * beta / reynolds) ** 0.7 + (0.0188 + 0.0063 * a) * beta**3.5 * (1e6 / reynolds) ** 0.3
    coefficient += (0.043 + 0.08 * math.exp(-10.0 * l1) - 0.123 * math.exp(-7.0 * l1)) * (1.0 - 0.11 * a) * beta**4 / (1.0 - beta**4) - 0.031 * (m2 - 0.8 * m2**1.1) * beta**1.3
    if pipe_diameter_m < 0.07112:
        coefficient += 0.011 * (0.75 - beta) * (2.8 - pipe_diameter_m / 0.0254)
    orifice_drop = density_kg_m3 / 2.0 * (mass_flow_kg_s / density_kg_m3 / (correction_factor * coefficient / (1.0 - beta**4) ** 0.5 * orifice_area)) ** 2 + density_kg_m3 * 9.8 * separation
    recovery = (1.0 - beta**4 * (1.0 - coefficient**2)) ** 0.5
    overall_drop = orifice_drop * (recovery - coefficient * beta**2) / (recovery + coefficient * beta**2)
    return OrificePressureDrop(reynolds, coefficient, orifice_drop, overall_drop)


def minor_loss_pressure_drop(loss_coefficient: float, density_kg_m3: float, velocity_m_s: float) -> float:
    """DWSIM fixed-K fitting pressure drop for one incompressible phase."""
    values = ((loss_coefficient, "loss coefficient"), (density_kg_m3, "density"), (velocity_m_s, "velocity"))
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0.0 for value, _ in values):
        raise ValidationError("fitting loss coefficient, density, and velocity must be finite and non-negative")
    return loss_coefficient * density_kg_m3 * velocity_m_s**2 / 2.0


def pipe_pressure_drop(
    diameter_m: float, length_m: float, elevation_m: float, roughness_m: float,
    volumetric_flow_m3_s: float, density_kg_m3: float, viscosity_pa_s: float,
) -> PipePressureDrop:
    """DWSIM Darcy-Weisbach pipe drop for one incompressible phase."""
    values = ((diameter_m, "pipe diameter"), (length_m, "pipe length"), (roughness_m, "pipe roughness"), (volumetric_flow_m3_s, "volumetric flow"), (density_kg_m3, "density"), (viscosity_pa_s, "viscosity"))
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) for value, _ in values):
        raise ValidationError("pipe inputs must be finite numbers")
    if diameter_m <= 0.0 or length_m <= 0.0 or volumetric_flow_m3_s < 0.0 or density_kg_m3 <= 0.0 or viscosity_pa_s <= 0.0 or roughness_m < 0.0:
        raise ValidationError("pipe dimensions and properties must be positive; roughness and flow may be zero")
    if not isinstance(elevation_m, (int, float)) or isinstance(elevation_m, bool) or not math.isfinite(elevation_m) or abs(elevation_m) > length_m:
        raise ValidationError("pipe elevation must be finite and no greater than pipe length")
    velocity = volumetric_flow_m3_s / (math.pi * diameter_m**2 / 4.0)
    reynolds = density_kg_m3 * velocity * diameter_m / viscosity_pa_s
    if reynolds == 0.0:
        friction = friction_drop = 0.0
    else:
        friction = _dwsim_friction_factor(reynolds, diameter_m, roughness_m)
        friction_drop = friction * length_m / diameter_m * velocity**2 * density_kg_m3 / 2.0
    static_drop = density_kg_m3 * 9.8 * elevation_m
    return PipePressureDrop(velocity, reynolds, friction, friction_drop, static_drop, friction_drop + static_drop)
