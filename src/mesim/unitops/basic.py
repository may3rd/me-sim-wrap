import math
from dataclasses import dataclass

from ..compounds import Compound, PRInteractions
from ..errors import ConvergenceError, ValidationError
from ..streams import EnergyStream, PhaseState, StreamState, flash_stream
from ..thermo.flash import mixture_heat_capacity, ph_flash, phase_enthalpy, tp_flash
from ..thermo.ideal import IdealCorrelations
from ..thermo.transport import TransportRecord, translated_vapor_density, vapor_transport


def _positive_finite(value: float, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0:
        raise ValidationError(f"{name} must be finite and positive")


def mix_streams(
    inlets: tuple[PhaseState, ...],
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    temperature_bracket_k: tuple[float, float],
) -> PhaseState:
    """Adiabatically mix flashed streams at a specified pressure."""
    if not inlets:
        raise ValidationError("mixer requires at least one inlet")
    _positive_finite(outlet_pressure_pa, "mixer outlet pressure")
    compound_ids = tuple(compound.id for compound in compounds)
    if not compound_ids or any(inlet.stream.compound_ids != compound_ids for inlet in inlets):
        raise ValidationError("mixer inlet compound IDs must exactly match the supplied compound order")
    if outlet_pressure_pa > min(inlet.stream.pressure_pa for inlet in inlets):
        raise ValidationError("mixer outlet pressure must not exceed the lowest inlet pressure")

    total_flow = math.fsum(inlet.stream.molar_flow_kmol_s for inlet in inlets)
    if total_flow == 0.0:
        raise ValidationError("mixer requires a positive total inlet flow")
    composition = tuple(
        math.fsum(inlet.stream.molar_flow_kmol_s * inlet.stream.composition[index] for inlet in inlets) / total_flow
        for index in range(len(compound_ids))
    )
    target_enthalpy = math.fsum(
        inlet.stream.molar_flow_kmol_s * inlet.enthalpy_j_per_kmol for inlet in inlets
    ) / total_flow
    result = ph_flash(
        compounds, composition, interactions, correlations, outlet_pressure_pa, target_enthalpy, temperature_bracket_k,
    )
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "mixer PH flash did not converge")
    stream = StreamState(
        result.temperature_k, outlet_pressure_pa, total_flow, compound_ids, composition,
        result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
    )
    return PhaseState(stream, result.flash, result.enthalpy_j_per_kmol)


def split_stream(stream: StreamState, fractions: tuple[float, ...]) -> tuple[StreamState, ...]:
    """Split a stream without changing its intensive state or composition."""
    if not fractions or any(
        isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not math.isfinite(fraction) or fraction < 0.0
        for fraction in fractions
    ) or not math.isclose(math.fsum(fractions), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValidationError("splitter fractions must be finite, non-negative, and sum to one")
    return tuple(
        StreamState(
            stream.temperature_k, stream.pressure_pa, stream.molar_flow_kmol_s * fraction,
            stream.compound_ids, stream.composition, stream.enthalpy_j_per_kmol, stream.vapor_fraction,
        )
        for fraction in fractions
    )


@dataclass(frozen=True, slots=True)
class ThermalOperationResult:
    outlet: PhaseState
    energy: EnergyStream


@dataclass(frozen=True, slots=True)
class HeatExchangerResult:
    hot_outlet: PhaseState
    cold_outlet: PhaseState
    heat_duty_w: float


@dataclass(frozen=True, slots=True)
class ShellTubeGeometry:
    shell_count: int
    tube_passes: int
    tube_inner_diameter_mm: float
    tube_outer_diameter_mm: float
    tube_length_m: float
    tube_count: int
    tube_pitch_mm: float


@dataclass(frozen=True, slots=True)
class VaporProperties:
    density_kg_m3: float
    viscosity_pa_s: float
    conductivity_w_m_k: float
    heat_capacity_j_kg_k: float


def shell_tube_vapor_properties(
    compounds: tuple[Compound, ...], composition: tuple[float, ...], correlations: tuple[IdealCorrelations, ...],
    transport: tuple[TransportRecord, ...], interactions: PRInteractions, temperature_k: float, pressure_pa: float,
) -> VaporProperties:
    """All-vapor PR properties required by shell-and-tube rating."""
    flash = tp_flash(compounds, composition, interactions, temperature_k, pressure_pa)
    if not flash.report.converged or flash.phase not in {"vapor", "single"} or flash.vapor_state is None:
        raise ValidationError("shell-tube rating requires a converged vapor state")
    density = translated_vapor_density(compounds, composition, temperature_k, pressure_pa, flash.vapor_state.compressibility)
    transport_values = vapor_transport(compounds, composition, transport, temperature_k, density)
    return VaporProperties(density, transport_values.dynamic_viscosity_pa_s, transport_values.thermal_conductivity_w_per_m_k, mixture_heat_capacity(compounds, composition, correlations, interactions, temperature_k, pressure_pa))



def shell_tube_area(geometry: ShellTubeGeometry) -> float:
    """DWSIM shell-and-tube external area, excluding two outer tube diameters."""
    if geometry.shell_count <= 0 or geometry.tube_passes <= 0 or geometry.tube_count <= 0:
        raise ValidationError("shell-tube counts and passes must be positive")
    for value, name in (
        (geometry.tube_inner_diameter_mm, "tube inner diameter"),
        (geometry.tube_outer_diameter_mm, "tube outer diameter"),
        (geometry.tube_length_m, "tube length"),
        (geometry.tube_pitch_mm, "tube pitch"),
    ):
        _positive_finite(value, name)
    outer = geometry.tube_outer_diameter_mm / 1000.0
    if geometry.tube_inner_diameter_mm >= geometry.tube_outer_diameter_mm or geometry.tube_pitch_mm / 1000.0 < outer or geometry.tube_length_m <= 2.0 * outer:
        raise ValidationError("shell-tube geometry has invalid tube diameters, pitch, or length")
    return geometry.tube_count * math.pi * outer * (geometry.tube_length_m - 2.0 * outer)


def shell_tube_tube_side(
    geometry: ShellTubeGeometry, mass_flow_kg_s: float, density_kg_m3: float, viscosity_pa_s: float,
    conductivity_w_m_k: float, heat_capacity_j_kg_k: float, roughness_mm: float, friction_multiplier: float,
) -> tuple[float, float, float, float]:
    """DWSIM tube-side pressure drop and Gnielinski heat-transfer coefficient."""
    shell_tube_area(geometry)
    for value, name in ((mass_flow_kg_s, "tube mass flow"), (density_kg_m3, "tube density"), (viscosity_pa_s, "tube viscosity"), (conductivity_w_m_k, "tube conductivity"), (heat_capacity_j_kg_k, "tube heat capacity"), (roughness_mm, "tube roughness"), (friction_multiplier, "tube friction multiplier")):
        _positive_finite(value, name)
    diameter = geometry.tube_inner_diameter_mm / 1000.0
    paths = geometry.tube_count / geometry.tube_passes
    velocity = mass_flow_kg_s / (density_kg_m3 * paths * math.pi * diameter**2 / 4.0)
    reynolds = density_kg_m3 * velocity * diameter / viscosity_pa_s
    if reynolds > 3250.0:
        a1 = math.log10(((roughness_mm / 1000.0 / diameter) ** 1.1096) / 2.8257 + (7.149 / reynolds) ** 0.8961)
        b1 = -2.0 * math.log10((roughness_mm / 1000.0 / diameter) / 3.7065 - 5.0452 * a1 / reynolds)
        friction = (1.0 / b1) ** 2
    else:
        friction = 64.0 / reynolds
    friction *= friction_multiplier
    prandtl = viscosity_pa_s * heat_capacity_j_kg_k / conductivity_w_m_k
    pressure_drop = friction * geometry.tube_length_m * geometry.tube_passes / diameter * velocity**2 * density_kg_m3 / 2.0
    coefficient = conductivity_w_m_k / diameter * (friction / 8.0) * reynolds * prandtl / (1.07 + 12.7 * (friction / 8.0) ** 0.5 * (prandtl ** (2.0 / 3.0) - 1.0))
    return reynolds, friction, pressure_drop, coefficient


def shell_tube_shell_side(
    geometry: ShellTubeGeometry, mass_flow_kg_s: float, density_kg_m3: float, viscosity_pa_s: float,
    conductivity_w_m_k: float, heat_capacity_j_kg_k: float, shell_diameter_mm: float,
    baffle_spacing_mm: float, baffle_cut_percent: float,
) -> tuple[float, float, float, float]:
    """DWSIM simplified Tinker shell-side correlation for tube layout 0."""
    shell_tube_area(geometry)
    for value, name in ((mass_flow_kg_s, "shell mass flow"), (density_kg_m3, "shell density"), (viscosity_pa_s, "shell viscosity"), (conductivity_w_m_k, "shell conductivity"), (heat_capacity_j_kg_k, "shell heat capacity"), (shell_diameter_mm, "shell diameter"), (baffle_spacing_mm, "baffle spacing")):
        _positive_finite(value, name)
    if not 0.0 < baffle_cut_percent < 100.0:
        raise ValidationError("shell baffle cut must be within (0, 100)")
    outer, pitch, shell, spacing = geometry.tube_outer_diameter_mm / 1000.0, geometry.tube_pitch_mm / 1000.0, shell_diameter_mm / 1000.0, baffle_spacing_mm / 1000.0
    if pitch / outer > 1.5:
        raise ValidationError("shell tube pitch ratio must not exceed 1.5")
    bundle = (1.1 * geometry.tube_count**0.5 - 1.0) * pitch + outer
    xx, yy = shell / spacing, pitch / outer
    nh = 0.9078565328950694 * xx**0.6633110612656448 * yy**-4.432976463965648
    y = 5.371855907482061 * xx**-0.33416765138071414 * yy**0.7267144209289168
    np = 0.5380765047084108 * xx**0.3761125784751041 * yy**-3.8741224386187474
    fp = 1.0 / (0.8 + np * (shell / pitch) ** 0.5)
    clearance = 0.97 * (pitch - outer) / pitch
    section = clearance * spacing * bundle
    pressure_section = section / fp
    if pressure_section <= 0.0:
        raise ValidationError("shell-tube pressure-drop flow area must be positive")
    mass_velocity = mass_flow_kg_s / pressure_section
    reynolds = mass_velocity * outer / viscosity_pa_s
    if reynolds < 100.0:
        friction = 276.46 * reynolds**-0.979
    elif reynolds < 1000.0:
        friction = 30.26 * reynolds**-0.523
    else:
        friction = 2.93 * reynolds**-0.186
    baffles = geometry.tube_length_m / spacing + 1.0
    pressure_drop = 4.0 * friction * mass_velocity**2 / (2.0 * density_kg_m3) * 1.154 * (1.0 - baffle_cut_percent / 100.0) * shell / pitch * baffles * (1.0 + y * pitch / shell) * geometry.shell_count
    corrected_section = section * 0.96 / (1.0 / (1.0 + nh * (shell / pitch) ** 0.5))
    corrected_reynolds = mass_flow_kg_s / corrected_section * outer / viscosity_pa_s
    jh = 0.497 * corrected_reynolds**0.54 if corrected_reynolds < 100.0 else 0.378 * corrected_reynolds**0.61
    coefficient = jh * conductivity_w_m_k * (viscosity_pa_s * heat_capacity_j_kg_k / conductivity_w_m_k) ** 0.34 / outer
    return reynolds, friction, pressure_drop, coefficient


def shell_tube_overall_coefficient(
    geometry: ShellTubeGeometry, tube_coefficient_w_m2_k: float, shell_coefficient_w_m2_k: float,
    tube_fouling_m2_k_w: float, shell_fouling_m2_k_w: float, tube_conductivity_w_m_k: float,
) -> float:
    """DWSIM overall HTC based on external tube area."""
    shell_tube_area(geometry)
    for value, name in ((tube_coefficient_w_m2_k, "tube coefficient"), (shell_coefficient_w_m2_k, "shell coefficient"), (tube_conductivity_w_m_k, "tube conductivity")):
        _positive_finite(value, name)
    if any(not math.isfinite(value) or value < 0.0 for value in (tube_fouling_m2_k_w, shell_fouling_m2_k_w)):
        raise ValidationError("shell-tube fouling resistances must be finite and non-negative")
    inner, outer = geometry.tube_inner_diameter_mm / 1000.0, geometry.tube_outer_diameter_mm / 1000.0
    resistance = outer / (tube_coefficient_w_m2_k * inner) + tube_fouling_m2_k_w * outer / inner + outer * math.log(outer / inner) / (2.0 * tube_conductivity_w_m_k) + shell_fouling_m2_k_w + 1.0 / shell_coefficient_w_m2_k
    return 1.0 / resistance


def shell_tube_lmtd_correction(ratio_r: float, ratio_p: float, shell_passes: int) -> float:
    """DWSIM shell-and-tube correction factor for counter-current LMTD."""
    if shell_passes <= 0 or not all(math.isfinite(value) and value > 0.0 for value in (ratio_r, ratio_p)):
        raise ValidationError("shell-tube LMTD ratios and shell passes must be positive")
    if math.isclose(ratio_r, 1.0, rel_tol=0.0, abs_tol=1e-12):
        s = ratio_p / (shell_passes * (1.0 - ratio_p) + ratio_p)
        factor = s * 2.0**0.5 / ((1.0 - s) * math.log((2.0 * (1.0 - s) + s * 2.0**0.5) / (2.0 * (1.0 - s) - s * 2.0**0.5)))
    else:
        alpha = ((1.0 - ratio_r * ratio_p) / (1.0 - ratio_p)) ** (1.0 / shell_passes)
        s = (alpha - 1.0) / (alpha - ratio_r)
        root = (ratio_r**2 + 1.0) ** 0.5
        factor = root * math.log((1.0 - s) / (1.0 - ratio_r * s)) / ((ratio_r - 1.0) * math.log((2.0 - s * (ratio_r + 1.0 - root)) / (2.0 - s * (ratio_r + 1.0 + root))))
    if not math.isfinite(factor) or not 0.0 < factor <= 1.0:
        raise ValidationError("shell-tube LMTD correction is outside its physical range")
    return factor


@dataclass(frozen=True, slots=True)
class ShellTubeRatingResult:
    result: "HeatExchangerResult"
    overall_coefficient_w_m2_k: float
    area_m2: float
    shell_reynolds: float
    tube_reynolds: float
    cold_pressure_drop_pa: float
    hot_pressure_drop_pa: float


def shell_tube_rating(
    hot_inlet: PhaseState, cold_inlet: PhaseState, compounds: tuple[Compound, ...], interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...], transport: tuple[TransportRecord, ...], geometry: ShellTubeGeometry,
    shell_diameter_mm: float, baffle_spacing_mm: float, baffle_cut_percent: float, tube_conductivity_w_m_k: float,
    tube_roughness_mm: float, tube_friction_multiplier: float,
) -> ShellTubeRatingResult:
    """All-vapor, layout-0, counter-current DWSIM shell-and-tube rating."""
    if hot_inlet.stream.temperature_k <= cold_inlet.stream.temperature_k:
        raise ValidationError("shell-tube hot inlet temperature must exceed cold inlet temperature")
    ids = tuple(compound.id for compound in compounds)
    if hot_inlet.stream.compound_ids != ids or cold_inlet.stream.compound_ids != ids or hot_inlet.stream.composition != cold_inlet.stream.composition:
        raise ValidationError("shell-tube inlets must share the supplied compound order and composition")
    area = shell_tube_area(geometry)
    molecular_weight = math.fsum(x * c.molecular_weight.value for x, c in zip(cold_inlet.stream.composition, compounds))
    cold_mass = cold_inlet.stream.molar_flow_kmol_s * molecular_weight
    hot_mass = hot_inlet.stream.molar_flow_kmol_s * molecular_weight
    duty = 0.0
    maximum_duty = _maximum_heat_duty(hot_inlet, cold_inlet, compounds, interactions, correlations) * (1.0 - 1e-9)
    correction = 1.0
    for _ in range(80):
        trial = heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, duty, (cold_inlet.stream.temperature_k, hot_inlet.stream.temperature_k)) if duty else HeatExchangerResult(hot_inlet, cold_inlet, 0.0)
        cold = shell_tube_vapor_properties(compounds, cold_inlet.stream.composition, correlations, transport, interactions, (cold_inlet.stream.temperature_k + trial.cold_outlet.stream.temperature_k) / 2, cold_inlet.stream.pressure_pa)
        hot = shell_tube_vapor_properties(compounds, hot_inlet.stream.composition, correlations, transport, interactions, (hot_inlet.stream.temperature_k + trial.hot_outlet.stream.temperature_k) / 2, hot_inlet.stream.pressure_pa)
        tube = shell_tube_tube_side(geometry, cold_mass, cold.density_kg_m3, cold.viscosity_pa_s, cold.conductivity_w_m_k, cold.heat_capacity_j_kg_k, tube_roughness_mm, tube_friction_multiplier)
        shell = shell_tube_shell_side(geometry, hot_mass, hot.density_kg_m3, hot.viscosity_pa_s, hot.conductivity_w_m_k, hot.heat_capacity_j_kg_k, shell_diameter_mm, baffle_spacing_mm, baffle_cut_percent)
        u = shell_tube_overall_coefficient(geometry, tube[3], shell[3], 0.0, 0.0, tube_conductivity_w_m_k)
        dt1, dt2 = hot_inlet.stream.temperature_k - trial.cold_outlet.stream.temperature_k, trial.hot_outlet.stream.temperature_k - cold_inlet.stream.temperature_k
        if dt1 <= 0 or dt2 <= 0:
            break
        lmtd = (dt1 + dt2) / 2 if math.isclose(dt1, dt2, rel_tol=1e-12) else (dt1 - dt2) / math.log(dt1 / dt2)
        r = (hot_inlet.stream.temperature_k - trial.hot_outlet.stream.temperature_k) / (trial.cold_outlet.stream.temperature_k - cold_inlet.stream.temperature_k) if duty else hot_mass / cold_mass
        p = (trial.cold_outlet.stream.temperature_k - cold_inlet.stream.temperature_k) / (hot_inlet.stream.temperature_k - cold_inlet.stream.temperature_k) if duty else 0.01
        try:
            correction = shell_tube_lmtd_correction(r, p, geometry.tube_passes)
        except (ValidationError, ValueError, ZeroDivisionError):
            pass
        target = min(maximum_duty, u * area * correction * lmtd)
        if abs(target - duty) <= 1e-5 * max(1.0, target):
            return ShellTubeRatingResult(trial, u, area, shell[0], tube[0], tube[2], shell[2])
        duty = 0.1 * target + 0.9 * duty
    raise ConvergenceError("shell-tube rating did not converge")


def heat_exchanger(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    heat_duty_w: float,
    temperature_bracket_k: tuple[float, float],
) -> HeatExchangerResult:
    """Transfer a specified positive duty from hot inlet to cold inlet at constant pressure."""
    _positive_finite(heat_duty_w, "heat exchanger duty")
    compound_ids = tuple(compound.id for compound in compounds)
    if hot_inlet.stream.compound_ids != compound_ids or cold_inlet.stream.compound_ids != compound_ids:
        raise ValidationError("heat exchanger inlet compound IDs must exactly match the supplied compound order")
    _positive_finite(hot_inlet.stream.molar_flow_kmol_s, "heat exchanger hot inlet flow")
    _positive_finite(cold_inlet.stream.molar_flow_kmol_s, "heat exchanger cold inlet flow")

    def outlet(inlet: PhaseState, target_enthalpy: float) -> PhaseState:
        result = ph_flash(
            compounds, inlet.stream.composition, interactions, correlations, inlet.stream.pressure_pa,
            target_enthalpy, temperature_bracket_k,
        )
        if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
            raise ConvergenceError(result.report.failure_reason or "heat exchanger PH flash did not converge")
        stream = StreamState(
            result.temperature_k, inlet.stream.pressure_pa, inlet.stream.molar_flow_kmol_s,
            inlet.stream.compound_ids, inlet.stream.composition,
            result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
        )
        return PhaseState(stream, result.flash, result.enthalpy_j_per_kmol)

    return HeatExchangerResult(
        outlet(hot_inlet, hot_inlet.enthalpy_j_per_kmol - heat_duty_w / hot_inlet.stream.molar_flow_kmol_s),
        outlet(cold_inlet, cold_inlet.enthalpy_j_per_kmol + heat_duty_w / cold_inlet.stream.molar_flow_kmol_s),
        heat_duty_w,
    )


def heat_exchanger_ua(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    overall_coefficient_w_m2_k: float,
    area_m2: float,
    temperature_bracket_k: tuple[float, float],
) -> HeatExchangerResult:
    """Calculate a counter-current, no-loss exchanger from a specified UA."""
    _positive_finite(overall_coefficient_w_m2_k, "heat exchanger overall coefficient")
    _positive_finite(area_m2, "heat exchanger area")
    if hot_inlet.stream.temperature_k <= cold_inlet.stream.temperature_k:
        raise ValidationError("heat exchanger hot inlet temperature must exceed cold inlet temperature")
    ua = overall_coefficient_w_m2_k * area_m2

    def lmtd(result: HeatExchangerResult) -> float | None:
        first = hot_inlet.stream.temperature_k - result.cold_outlet.stream.temperature_k
        second = result.hot_outlet.stream.temperature_k - cold_inlet.stream.temperature_k
        if first <= 0.0 or second <= 0.0:
            return None
        if math.isclose(first, second, rel_tol=1e-12, abs_tol=1e-12):
            return (first + second) / 2.0
        return (first - second) / math.log(first / second)

    lower = 0.0
    upper = ua * (hot_inlet.stream.temperature_k - cold_inlet.stream.temperature_k)
    for _ in range(80):
        trial = heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, upper, temperature_bracket_k)
        difference = lmtd(trial)
        if difference is None or upper >= ua * difference:
            break
        upper *= 2.0
    else:
        raise ConvergenceError("heat exchanger UA duty bracket did not converge")

    for _ in range(60):
        duty = (lower + upper) / 2.0
        trial = heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, duty, temperature_bracket_k)
        difference = lmtd(trial)
        if difference is None or duty >= ua * difference:
            upper = duty
        else:
            lower = duty
    return heat_exchanger(hot_inlet, cold_inlet, compounds, interactions, correlations, (lower + upper) / 2.0, temperature_bracket_k)


def _maximum_heat_duty(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
) -> float:
    compound_ids = tuple(compound.id for compound in compounds)
    if hot_inlet.stream.compound_ids != compound_ids or cold_inlet.stream.compound_ids != compound_ids:
        raise ValidationError("heat exchanger inlet compound IDs must exactly match the supplied compound order")
    _positive_finite(hot_inlet.stream.molar_flow_kmol_s, "heat exchanger hot inlet flow")
    _positive_finite(cold_inlet.stream.molar_flow_kmol_s, "heat exchanger cold inlet flow")
    if hot_inlet.stream.temperature_k <= cold_inlet.stream.temperature_k:
        raise ValidationError("heat exchanger hot inlet temperature must exceed cold inlet temperature")
    hot_at_cold_inlet = flash_stream(
        StreamState(cold_inlet.stream.temperature_k, hot_inlet.stream.pressure_pa, hot_inlet.stream.molar_flow_kmol_s, hot_inlet.stream.compound_ids, hot_inlet.stream.composition),
        compounds, interactions, correlations,
    )
    cold_at_hot_inlet = flash_stream(
        StreamState(hot_inlet.stream.temperature_k, cold_inlet.stream.pressure_pa, cold_inlet.stream.molar_flow_kmol_s, cold_inlet.stream.compound_ids, cold_inlet.stream.composition),
        compounds, interactions, correlations,
    )
    return min(
        hot_inlet.stream.molar_flow_kmol_s * (hot_inlet.enthalpy_j_per_kmol - hot_at_cold_inlet.enthalpy_j_per_kmol),
        cold_inlet.stream.molar_flow_kmol_s * (cold_at_hot_inlet.enthalpy_j_per_kmol - cold_inlet.enthalpy_j_per_kmol),
    )


def heat_exchanger_efficiency(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    efficiency_percent: float,
) -> HeatExchangerResult:
    """Calculate a no-loss exchanger from DWSIM's maximum-heat-transfer efficiency."""
    if isinstance(efficiency_percent, bool) or not isinstance(efficiency_percent, (int, float)) or not math.isfinite(efficiency_percent) or not 0.0 <= efficiency_percent <= 100.0:
        raise ValidationError("heat exchanger efficiency must be finite and within [0, 100]")
    maximum_duty = _maximum_heat_duty(hot_inlet, cold_inlet, compounds, interactions, correlations)
    if efficiency_percent == 0.0:
        return HeatExchangerResult(hot_inlet, cold_inlet, 0.0)
    return heat_exchanger(
        hot_inlet, cold_inlet, compounds, interactions, correlations,
        maximum_duty * efficiency_percent / 100.0,
        (cold_inlet.stream.temperature_k, hot_inlet.stream.temperature_k),
    )


def heat_exchanger_pinch(
    hot_inlet: PhaseState,
    cold_inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    minimum_approach_k: float,
) -> HeatExchangerResult:
    """Solve a phase-free counter-current no-loss minimum temperature approach."""
    if isinstance(minimum_approach_k, bool) or not isinstance(minimum_approach_k, (int, float)) or not math.isfinite(minimum_approach_k) or minimum_approach_k <= 0.0:
        raise ValidationError("heat exchanger minimum approach must be finite and positive")
    inlet_approach = hot_inlet.stream.temperature_k - cold_inlet.stream.temperature_k
    if minimum_approach_k >= inlet_approach:
        raise ValidationError("heat exchanger minimum approach must be below the inlet temperature difference")
    maximum_duty = _maximum_heat_duty(hot_inlet, cold_inlet, compounds, interactions, correlations)
    lower, upper = 0.0, maximum_duty
    for _ in range(60):
        duty = (lower + upper) / 2.0
        trial = heat_exchanger(
            hot_inlet, cold_inlet, compounds, interactions, correlations, duty,
            (cold_inlet.stream.temperature_k, hot_inlet.stream.temperature_k),
        )
        approach = min(
            hot_inlet.stream.temperature_k - trial.cold_outlet.stream.temperature_k,
            trial.hot_outlet.stream.temperature_k - cold_inlet.stream.temperature_k,
        )
        if approach > minimum_approach_k:
            lower = duty
        else:
            upper = duty
    return heat_exchanger(
        hot_inlet, cold_inlet, compounds, interactions, correlations, (lower + upper) / 2.0,
        (cold_inlet.stream.temperature_k, hot_inlet.stream.temperature_k),
    )


def _thermal_operation(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_temperature_k: float,
    kind: str,
) -> ThermalOperationResult:
    _positive_finite(outlet_temperature_k, f"{kind} outlet temperature")
    outlet = flash_stream(
        StreamState(
            outlet_temperature_k, inlet.stream.pressure_pa, inlet.stream.molar_flow_kmol_s,
            inlet.stream.compound_ids, inlet.stream.composition,
        ),
        compounds, interactions, correlations,
    )
    duty = inlet.stream.molar_flow_kmol_s * (outlet.enthalpy_j_per_kmol - inlet.enthalpy_j_per_kmol)
    if (kind == "heater" and duty < 0.0) or (kind == "cooler" and duty > 0.0):
        raise ValidationError(f"{kind} outlet temperature gives an invalid duty sign")
    return ThermalOperationResult(outlet, EnergyStream(duty))


def heater(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_temperature_k: float,
) -> ThermalOperationResult:
    """Heat a flashed stream at constant pressure to a specified temperature."""
    return _thermal_operation(inlet, compounds, interactions, correlations, outlet_temperature_k, "heater")


def cooler(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_temperature_k: float,
) -> ThermalOperationResult:
    """Cool a flashed stream at constant pressure to a specified temperature."""
    return _thermal_operation(inlet, compounds, interactions, correlations, outlet_temperature_k, "cooler")


def valve(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    interactions: PRInteractions,
    correlations: tuple[IdealCorrelations, ...],
    outlet_pressure_pa: float,
    temperature_bracket_k: tuple[float, float],
) -> PhaseState:
    """Throttle a flashed stream isenthalpically to a lower pressure."""
    _positive_finite(outlet_pressure_pa, "valve outlet pressure")
    if outlet_pressure_pa >= inlet.stream.pressure_pa:
        raise ValidationError("valve outlet pressure must be below inlet pressure")
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("valve inlet compound IDs must exactly match the supplied compound order")
    result = ph_flash(
        compounds, inlet.stream.composition, interactions, correlations, outlet_pressure_pa,
        inlet.enthalpy_j_per_kmol, temperature_bracket_k,
    )
    if not result.report.converged or result.flash is None or result.enthalpy_j_per_kmol is None:
        raise ConvergenceError(result.report.failure_reason or "valve PH flash did not converge")
    stream = StreamState(
        result.temperature_k, outlet_pressure_pa, inlet.stream.molar_flow_kmol_s, inlet.stream.compound_ids,
        inlet.stream.composition, result.enthalpy_j_per_kmol, result.flash.vapor_fraction,
    )
    return PhaseState(stream, result.flash, result.enthalpy_j_per_kmol)


@dataclass(frozen=True, slots=True)
class EquilibriumSeparatorResult:
    liquid: StreamState | None
    vapor: StreamState | None


def equilibrium_separator(
    inlet: PhaseState,
    compounds: tuple[Compound, ...],
    correlations: tuple[IdealCorrelations, ...],
) -> EquilibriumSeparatorResult:
    """Route phases from an existing converged flash without recalculation."""
    flash = inlet.flash
    if tuple(compound.id for compound in compounds) != inlet.stream.compound_ids:
        raise ValidationError("separator inlet compound IDs must exactly match the supplied compound order")
    if not flash.report.converged:
        raise ConvergenceError(flash.report.failure_reason or "separator requires a converged flash")

    def phase_stream(composition: tuple[float, ...], state, flow: float, vapor_fraction: float) -> StreamState:
        return StreamState(
            flash.temperature_k, flash.pressure_pa, flow, inlet.stream.compound_ids, composition,
            phase_enthalpy(compounds, composition, correlations, flash.temperature_k, state), vapor_fraction,
        )

    if flash.phase == "liquid":
        return EquilibriumSeparatorResult(
            phase_stream(flash.liquid_composition, flash.liquid_state, inlet.stream.molar_flow_kmol_s, 0.0), None,
        )
    if flash.phase == "vapor":
        return EquilibriumSeparatorResult(
            None, phase_stream(flash.vapor_composition, flash.vapor_state, inlet.stream.molar_flow_kmol_s, 1.0),
        )
    if flash.phase == "two-phase" and flash.vapor_fraction is not None:
        return EquilibriumSeparatorResult(
            phase_stream(
                flash.liquid_composition, flash.liquid_state,
                inlet.stream.molar_flow_kmol_s * (1.0 - flash.vapor_fraction), 0.0,
            ),
            phase_stream(
                flash.vapor_composition, flash.vapor_state,
                inlet.stream.molar_flow_kmol_s * flash.vapor_fraction, 1.0,
            ),
        )
    raise ValidationError(f"unsupported separator flash phase: {flash.phase}")
