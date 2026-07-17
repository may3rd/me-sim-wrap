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
class PipePressureProfileResult:
    segment_results: tuple[PipePressureDrop, ...]
    outlet_pressure_pa: float
    friction_drop_pa: float
    static_drop_pa: float
    total_drop_pa: float


@dataclass(frozen=True, slots=True)
class PipeThermalResult:
    area_m2: float
    outlet_temperature_k: float
    heat_transfer_w: float


@dataclass(frozen=True, slots=True)
class PipeThermalProfileResult:
    segment_results: tuple[PipeThermalResult, ...]
    total_area_m2: float
    outlet_temperature_k: float
    heat_transfer_w: float


@dataclass(frozen=True, slots=True)
class LiquidPipeProfileResult:
    pressure: PipePressureProfileResult
    thermal: PipeThermalProfileResult


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


@dataclass(frozen=True, slots=True)
class TwoPhasePipeSegment:
    length_m: float
    elevation_m: float
    vapor_flow_m3_s: float
    liquid_flow_m3_s: float
    vapor_density_kg_m3: float
    liquid_density_kg_m3: float
    vapor_viscosity_pa_s: float
    liquid_viscosity_pa_s: float
    surface_tension_n_m: float = 0.0


@dataclass(frozen=True, slots=True)
class TwoPhasePressureProfileResult:
    segment_results: tuple[BeggsBrillPressureDrop | LockhartMartinelliPressureDrop, ...]
    outlet_pressure_pa: float
    friction_drop_pa: float
    static_drop_pa: float
    total_drop_pa: float


@dataclass(frozen=True, slots=True)
class ApiRp520Area:
    relieving_pressure_pa: float
    critical_pressure_pa: float
    choked: bool
    required_area_in2: float
    standard_orifice: str
    standard_area_in2: float


_API_526_ORIFICES = ("D", 0.11), ("E", 0.196), ("F", 0.307), ("G", 0.503), ("H", 0.785), ("J", 1.287), ("K", 1.838), ("L", 2.853), ("M", 3.6), ("N", 4.34), ("P", 6.38), ("Q", 11.05), ("R", 16.0), ("T", 26.0)


def api_rp520_vapor_required_area(
    temperature_k: float, set_pressure_pa: float, back_pressure_pa: float, mass_flow_kg_s: float,
    compressibility: float, molecular_weight_kg_kmol: float, heat_capacity_ratio: float,
    overpressure_percent: float, discharge_coefficient: float, backpressure_coefficient: float,
    installation_coefficient: float,
) -> ApiRp520Area:
    """DWSIM's API RP 520 vapor sizing utility; source-equation parity only."""
    positive = (temperature_k, set_pressure_pa, back_pressure_pa, mass_flow_kg_s, compressibility, molecular_weight_kg_kmol, heat_capacity_ratio, discharge_coefficient, backpressure_coefficient, installation_coefficient)
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0 for value in positive):
        raise ValidationError("API RP 520 vapor inputs must be finite and positive")
    if heat_capacity_ratio <= 1.0:
        raise ValidationError("API RP 520 heat-capacity ratio must exceed one")
    if isinstance(overpressure_percent, bool) or not isinstance(overpressure_percent, (int, float)) or not math.isfinite(overpressure_percent) or overpressure_percent < 0.0:
        raise ValidationError("API RP 520 overpressure must be finite and non-negative")
    relieving_pressure_kpa = (set_pressure_pa * 1.033 / 101325.0 * (1.0 + overpressure_percent / 100.0) + 1.033) * 101.325 / 1.033
    back_pressure_kpa = (back_pressure_pa * 1.033 / 101325.0 + 1.033) * 101.325 / 1.033
    if back_pressure_kpa >= relieving_pressure_kpa:
        raise ValidationError("API RP 520 back pressure must be below relieving pressure")
    critical_pressure_kpa = relieving_pressure_kpa * (2.0 / (heat_capacity_ratio + 1.0)) ** (heat_capacity_ratio / (heat_capacity_ratio - 1.0))
    mass_flow_kg_h = mass_flow_kg_s * 3600.0
    if back_pressure_kpa <= critical_pressure_kpa:
        coefficient = 520.0 * (heat_capacity_ratio * (2.0 / (heat_capacity_ratio + 1.0)) ** ((heat_capacity_ratio + 1.0) / (heat_capacity_ratio - 1.0))) ** 0.5
        area = 13160.0 * mass_flow_kg_h * (temperature_k * compressibility / molecular_weight_kg_kmol) ** 0.5 / (coefficient * discharge_coefficient * relieving_pressure_kpa * backpressure_coefficient * installation_coefficient)
        choked = True
    else:
        ratio = back_pressure_kpa / relieving_pressure_kpa
        coefficient = (heat_capacity_ratio / (heat_capacity_ratio - 1.0) * ratio ** (2.0 / heat_capacity_ratio) * (1.0 - ratio ** ((heat_capacity_ratio - 1.0) / heat_capacity_ratio)) / (1.0 - ratio)) ** 0.5
        area = 17.9 * mass_flow_kg_h / (coefficient * discharge_coefficient * installation_coefficient) * (compressibility * temperature_k / (molecular_weight_kg_kmol * relieving_pressure_kpa * (relieving_pressure_kpa - back_pressure_kpa))) ** 0.5
        choked = False
    required_area = area * 0.00155
    for designation, standard_area in _API_526_ORIFICES:
        if required_area <= standard_area:
            return ApiRp520Area(relieving_pressure_kpa * 1000.0, critical_pressure_kpa * 1000.0, choked, required_area, designation, standard_area)
    raise ValidationError("API RP 520 required area exceeds DWSIM's API 526 T orifice")


def api_rp520_liquid_required_area(
    set_pressure_pa: float, back_pressure_pa: float, volumetric_flow_m3_s: float,
    density_kg_m3: float, viscosity_pa_s: float, overpressure_percent: float,
    discharge_coefficient: float, installation_coefficient: float,
) -> ApiRp520Area:
    """DWSIM's API RP 520 liquid sizing utility; source-equation parity only."""
    positive = (set_pressure_pa, back_pressure_pa, volumetric_flow_m3_s, density_kg_m3, viscosity_pa_s, discharge_coefficient, installation_coefficient)
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0 for value in positive):
        raise ValidationError("API RP 520 liquid inputs must be finite and positive")
    if isinstance(overpressure_percent, bool) or not isinstance(overpressure_percent, (int, float)) or not math.isfinite(overpressure_percent) or overpressure_percent < 0.0:
        raise ValidationError("API RP 520 overpressure must be finite and non-negative")
    relieving_pressure_kpa = (set_pressure_pa * 1.033 / 101325.0 * (1.0 + overpressure_percent / 100.0) + 1.033) * 101.325 / 1.033 - 101.325
    back_pressure_kpa = (back_pressure_pa * 1.033 / 101325.0 + 1.033) * 101.325 / 1.033 - 101.325
    if back_pressure_kpa >= relieving_pressure_kpa:
        raise ValidationError("API RP 520 back pressure must be below relieving pressure")
    flow = 16.6667 / 24.0 * volumetric_flow_m3_s * 86400.0
    density = density_kg_m3 / 1000.0
    correction = 1.0
    # ponytail: DWSIM passes Pa.s into its cP-labelled utility equation; retain source behavior for parity.
    for _ in range(100):
        required_area = 11.78 * flow / (discharge_coefficient * installation_coefficient * correction) * (density / (relieving_pressure_kpa - back_pressure_kpa)) ** 0.5 * 0.00155
        standard = next(((letter, area) for letter, area in _API_526_ORIFICES if required_area <= area), None)
        if standard is None:
            raise ValidationError("API RP 520 required area exceeds DWSIM's API 526 T orifice")
        reynolds = flow * 18800.0 * density / (viscosity_pa_s * (standard[1] / 0.00155) ** 0.5)
        correction = (0.9935 + 2.878 / reynolds**0.5 + 342.75 / reynolds**1.5) ** -1
        required_area = 11.78 * flow / (discharge_coefficient * installation_coefficient * correction) * (density / (relieving_pressure_kpa - back_pressure_kpa)) ** 0.5 * 0.00155
        if required_area <= standard[1]:
            return ApiRp520Area(relieving_pressure_kpa * 1000.0, 0.0, False, required_area, standard[0], standard[1])
    raise ValidationError("API RP 520 liquid viscosity correction did not converge")


def api_rp520_two_phase_required_area(
    set_pressure_pa: float, back_pressure_pa: float, volumetric_flow_m3_s: float,
    vapor_mass_fraction: float, vapor_density_kg_m3: float, mixture_density_kg_m3: float,
    mixture_density_at_90_percent_pressure_kg_m3: float, overpressure_percent: float,
    discharge_coefficient: float, backpressure_coefficient: float, installation_coefficient: float,
) -> ApiRp520Area:
    """DWSIM's API RP 520 two-phase sizing utility; source-equation parity only."""
    positive = (set_pressure_pa, back_pressure_pa, volumetric_flow_m3_s, vapor_density_kg_m3, mixture_density_kg_m3, mixture_density_at_90_percent_pressure_kg_m3, discharge_coefficient, backpressure_coefficient, installation_coefficient)
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0 for value in positive):
        raise ValidationError("API RP 520 two-phase inputs must be finite and positive")
    if not isinstance(vapor_mass_fraction, (int, float)) or isinstance(vapor_mass_fraction, bool) or not math.isfinite(vapor_mass_fraction) or not 0.0 < vapor_mass_fraction < 1.0:
        raise ValidationError("API RP 520 two-phase vapor mass fraction must be between zero and one")
    if isinstance(overpressure_percent, bool) or not isinstance(overpressure_percent, (int, float)) or not math.isfinite(overpressure_percent) or overpressure_percent < 0.0:
        raise ValidationError("API RP 520 overpressure must be finite and non-negative")
    relieving_pressure_pa = set_pressure_pa * (1.0 + overpressure_percent / 100.0)
    if back_pressure_pa >= relieving_pressure_pa:
        raise ValidationError("API RP 520 back pressure must be below relieving pressure")
    pressure = relieving_pressure_pa * 1.033 / 101325.0 * 14.22
    back_pressure = back_pressure_pa * 1.033 / 101325.0 * 14.22
    specific_vapor_volume = 16.0185 / vapor_density_kg_m3
    specific_volume = 16.0185 / mixture_density_kg_m3
    specific_volume_90 = 16.0185 / mixture_density_at_90_percent_pressure_kg_m3
    expansion = 9.0 * (specific_volume_90 / specific_volume - 1.0)
    if expansion <= 0.0:
        raise ValidationError("API RP 520 two-phase density at 90 percent pressure must be lower")

    def residual(eta: float) -> float:
        return eta**2 + (expansion**2 - 2.0 * expansion) * (1.0 - eta)**2 + 2.0 * expansion**2 * math.log(eta) + 2.0 * expansion**2 * (1.0 - eta)

    low, high = 1e-12, 1.0
    if residual(low) * residual(high) >= 0.0:
        raise ValidationError("API RP 520 two-phase critical-flow equation has no bracket")
    for _ in range(100):
        middle = (low + high) / 2.0
        if residual(low) * residual(middle) <= 0.0:
            high = middle
        else:
            low = middle
    critical_ratio = (low + high) / 2.0
    critical_pressure = critical_ratio * pressure
    if critical_pressure >= back_pressure:
        flow_coefficient = 68.09 * critical_ratio * (pressure / (specific_volume * expansion)) ** 0.5
        choked = True
    else:
        flow_coefficient = 68.09 * (-2.0 * (expansion / math.log(back_pressure / pressure) + (expansion - 1.0) * (1.0 - back_pressure / pressure))) ** 0.5 * (pressure / specific_volume) ** 0.5 / (expansion * (pressure / back_pressure - 1.0) + 1.0)
        choked = False
    required_area = 0.04 * volumetric_flow_m3_s * 86400.0 * 2.20462 / (discharge_coefficient * backpressure_coefficient * installation_coefficient * flow_coefficient)
    for designation, standard_area in _API_526_ORIFICES:
        if required_area <= standard_area:
            return ApiRp520Area(relieving_pressure_pa, critical_ratio * relieving_pressure_pa, choked, required_area, designation, standard_area)
    raise ValidationError("API RP 520 required area exceeds DWSIM's API 526 T orifice")


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


def _two_phase_pressure_profile(
    inlet_pressure_pa: float,
    results: tuple[BeggsBrillPressureDrop | LockhartMartinelliPressureDrop, ...],
) -> TwoPhasePressureProfileResult:
    if (
        isinstance(inlet_pressure_pa, bool)
        or not isinstance(inlet_pressure_pa, (int, float))
        or not math.isfinite(inlet_pressure_pa)
        or inlet_pressure_pa <= 0.0
    ):
        raise ValidationError("pipe inlet pressure must be finite and positive")
    friction_drop = math.fsum(item.friction_drop_pa for item in results)
    static_drop = math.fsum(item.static_drop_pa for item in results)
    total_drop = friction_drop + static_drop
    if total_drop >= inlet_pressure_pa:
        raise ValidationError("pipe pressure drop must remain below inlet absolute pressure")
    return TwoPhasePressureProfileResult(
        results, inlet_pressure_pa - total_drop, friction_drop, static_drop, total_drop,
    )


def beggs_brill_pressure_drop_profile(
    inlet_pressure_pa: float, diameter_m: float, roughness_m: float,
    segments: tuple[TwoPhasePipeSegment, ...],
) -> TwoPhasePressureProfileResult:
    """Aggregate DWSIM Beggs-Brill drops over supplied two-phase segment states."""
    try:
        segment_states = tuple(segments)
    except TypeError as exc:
        raise ValidationError("Beggs-Brill profile segments must be a finite sequence") from exc
    if not segment_states or any(not isinstance(item, TwoPhasePipeSegment) for item in segment_states):
        raise ValidationError("Beggs-Brill profile requires at least one two-phase segment")
    results = tuple(
        beggs_brill_pressure_drop(
            diameter_m, item.length_m, item.elevation_m, roughness_m,
            item.vapor_flow_m3_s, item.liquid_flow_m3_s, item.vapor_density_kg_m3,
            item.liquid_density_kg_m3, item.vapor_viscosity_pa_s,
            item.liquid_viscosity_pa_s, item.surface_tension_n_m,
        )
        for item in segment_states
    )
    return _two_phase_pressure_profile(inlet_pressure_pa, results)


def lockhart_martinelli_pressure_drop_profile(
    inlet_pressure_pa: float, diameter_m: float, roughness_m: float,
    segments: tuple[TwoPhasePipeSegment, ...],
) -> TwoPhasePressureProfileResult:
    """Aggregate DWSIM Lockhart-Martinelli drops over supplied two-phase segment states."""
    try:
        segment_states = tuple(segments)
    except TypeError as exc:
        raise ValidationError("Lockhart-Martinelli profile segments must be a finite sequence") from exc
    if not segment_states or any(not isinstance(item, TwoPhasePipeSegment) for item in segment_states):
        raise ValidationError("Lockhart-Martinelli profile requires at least one two-phase segment")
    results = tuple(
        lockhart_martinelli_pressure_drop(
            diameter_m, item.length_m, item.elevation_m, roughness_m,
            item.vapor_flow_m3_s, item.liquid_flow_m3_s, item.vapor_density_kg_m3,
            item.liquid_density_kg_m3, item.vapor_viscosity_pa_s, item.liquid_viscosity_pa_s,
        )
        for item in segment_states
    )
    return _two_phase_pressure_profile(inlet_pressure_pa, results)


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


def pipe_defined_htc_heat_transfer(
    inlet_temperature_k: float, external_temperature_k: float, overall_htc_w_m2_k: float,
    outer_diameter_m: float, length_m: float, mass_flow_kg_s: float, heat_capacity_j_kg_k: float,
) -> PipeThermalResult:
    """Constant-property pipe heat transfer for a defined external temperature and overall HTC."""
    positive = (
        (inlet_temperature_k, "pipe inlet temperature"),
        (external_temperature_k, "pipe external temperature"),
        (outer_diameter_m, "pipe outer diameter"),
        (length_m, "pipe length"),
        (mass_flow_kg_s, "pipe mass flow"),
        (heat_capacity_j_kg_k, "pipe heat capacity"),
    )
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0
        for value, _ in positive
    ):
        raise ValidationError("pipe thermal temperatures, geometry, flow, and heat capacity must be finite and positive")
    if (
        isinstance(overall_htc_w_m2_k, bool)
        or not isinstance(overall_htc_w_m2_k, (int, float))
        or not math.isfinite(overall_htc_w_m2_k)
        or overall_htc_w_m2_k < 0.0
    ):
        raise ValidationError("pipe overall heat-transfer coefficient must be finite and non-negative")

    area = math.pi * outer_diameter_m * length_m
    if overall_htc_w_m2_k == 0.0 or inlet_temperature_k == external_temperature_k:
        return PipeThermalResult(area, inlet_temperature_k, 0.0)
    capacity_rate = mass_flow_kg_s * heat_capacity_j_kg_k
    effectiveness = -math.expm1(-overall_htc_w_m2_k * area / capacity_rate)
    temperature_change = (external_temperature_k - inlet_temperature_k) * effectiveness
    return PipeThermalResult(
        area,
        inlet_temperature_k + temperature_change,
        capacity_rate * temperature_change,
    )


def pipe_defined_htc_profile(
    inlet_temperature_k: float, external_temperature_k: float, overall_htc_w_m2_k: float,
    outer_diameter_m: float, segment_lengths_m: tuple[float, ...], mass_flow_kg_s: float,
    heat_capacities_j_kg_k: tuple[float, ...],
) -> PipeThermalProfileResult:
    """Advance a defined-HTC thermal state across supplied pipe segments."""
    try:
        lengths = tuple(segment_lengths_m)
        heat_capacities = tuple(heat_capacities_j_kg_k)
    except TypeError as exc:
        raise ValidationError("pipe thermal segment lengths and heat capacities must be finite sequences") from exc
    if not lengths:
        raise ValidationError("pipe thermal profile must contain at least one segment")
    if len(lengths) != len(heat_capacities):
        raise ValidationError("pipe thermal segment lengths and heat capacities must have equal counts")

    temperature = inlet_temperature_k
    results = []
    for length_m, heat_capacity_j_kg_k in zip(lengths, heat_capacities):
        result = pipe_defined_htc_heat_transfer(
            temperature, external_temperature_k, overall_htc_w_m2_k, outer_diameter_m,
            length_m, mass_flow_kg_s, heat_capacity_j_kg_k,
        )
        results.append(result)
        temperature = result.outlet_temperature_k
    segment_results = tuple(results)
    return PipeThermalProfileResult(
        segment_results,
        math.fsum(item.area_m2 for item in segment_results),
        temperature,
        math.fsum(item.heat_transfer_w for item in segment_results),
    )


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


def pipe_pressure_drop_profile(
    inlet_pressure_pa: float, diameter_m: float, roughness_m: float,
    segment_lengths_m: tuple[float, ...], segment_elevations_m: tuple[float, ...],
    volumetric_flows_m3_s: tuple[float, ...], densities_kg_m3: tuple[float, ...],
    viscosities_pa_s: tuple[float, ...],
) -> PipePressureProfileResult:
    """Aggregate a single-phase pressure profile from supplied segment states."""
    if (
        isinstance(inlet_pressure_pa, bool)
        or not isinstance(inlet_pressure_pa, (int, float))
        or not math.isfinite(inlet_pressure_pa)
        or inlet_pressure_pa <= 0.0
    ):
        raise ValidationError("pipe inlet pressure must be finite and positive")
    try:
        sequences = tuple(
            tuple(values)
            for values in (
                segment_lengths_m, segment_elevations_m, volumetric_flows_m3_s,
                densities_kg_m3, viscosities_pa_s,
            )
        )
    except TypeError as exc:
        raise ValidationError("pipe pressure profile inputs must be finite sequences") from exc
    if not sequences[0]:
        raise ValidationError("pipe pressure profile must contain at least one segment")
    if any(len(values) != len(sequences[0]) for values in sequences[1:]):
        raise ValidationError("pipe pressure profile inputs must have equal segment counts")

    results = tuple(
        pipe_pressure_drop(diameter_m, length_m, elevation_m, roughness_m, flow_m3_s, density_kg_m3, viscosity_pa_s)
        for length_m, elevation_m, flow_m3_s, density_kg_m3, viscosity_pa_s in zip(*sequences)
    )
    friction_drop = math.fsum(item.friction_drop_pa for item in results)
    static_drop = math.fsum(item.static_drop_pa for item in results)
    total_drop = friction_drop + static_drop
    if total_drop >= inlet_pressure_pa:
        raise ValidationError("pipe pressure drop must remain below inlet absolute pressure")
    return PipePressureProfileResult(
        results,
        inlet_pressure_pa - total_drop,
        friction_drop,
        static_drop,
        total_drop,
    )


def liquid_pipe_supplied_state_profile(
    inlet_pressure_pa: float, inlet_temperature_k: float,
    inner_diameter_m: float, outer_diameter_m: float, roughness_m: float,
    hydraulic_segment_lengths_m: tuple[float, ...], hydraulic_segment_elevations_m: tuple[float, ...],
    liquid_flows_m3_s: tuple[float, ...], liquid_densities_kg_m3: tuple[float, ...],
    liquid_viscosities_pa_s: tuple[float, ...], external_temperature_k: float,
    overall_htc_w_m2_k: float, thermal_segment_lengths_m: tuple[float, ...],
    mass_flow_kg_s: float, heat_capacities_j_kg_k: tuple[float, ...],
) -> LiquidPipeProfileResult:
    """Couple supplied-state liquid pressure and defined-HTC thermal profiles."""
    try:
        flows = tuple(liquid_flows_m3_s)
        densities = tuple(liquid_densities_kg_m3)
    except TypeError as exc:
        raise ValidationError("liquid pipe flow and density inputs must be finite sequences") from exc
    pressure = pipe_pressure_drop_profile(
        inlet_pressure_pa, inner_diameter_m, roughness_m,
        hydraulic_segment_lengths_m, hydraulic_segment_elevations_m,
        flows, densities, liquid_viscosities_pa_s,
    )
    thermal = pipe_defined_htc_profile(
        inlet_temperature_k, external_temperature_k, overall_htc_w_m2_k,
        outer_diameter_m, thermal_segment_lengths_m, mass_flow_kg_s,
        heat_capacities_j_kg_k,
    )
    if outer_diameter_m <= inner_diameter_m:
        raise ValidationError("pipe outer diameter must exceed inner diameter")
    if any(
        not math.isclose(flow_m3_s * density_kg_m3, mass_flow_kg_s, rel_tol=1e-6, abs_tol=1e-12)
        for flow_m3_s, density_kg_m3 in zip(flows, densities)
    ):
        raise ValidationError("pipe hydraulic and thermal profiles must use one coherent mass-flow basis")
    return LiquidPipeProfileResult(pressure, thermal)
