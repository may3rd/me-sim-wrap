"""Classic Soave-Redlich-Kwong pure and mixture states."""

import math
from dataclasses import dataclass, replace

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError
from .advanced_cubic import EvaluatedAdvancedCubicInteractions
from .peng_robinson import R, _cubic_real_roots, _positive_exp, _reduced
from .flash import SolverReport, rachford_rice


@dataclass(frozen=True, slots=True)
class SRKParameters:
    kappa: float
    alpha: float
    a_pa_m6_per_kmol2: float
    b_m3_per_kmol: float
    da_dtemperature: float


@dataclass(frozen=True, slots=True)
class SRKMixtureParameters:
    a_pa_m6_per_kmol2: float
    b_m3_per_kmol: float
    da_dtemperature: float


@dataclass(frozen=True, slots=True)
class SRKState:
    phase: str
    compressibility: float
    fugacity_coefficient: float
    density_kg_per_m3: float
    departure_enthalpy_j_per_kmol: float
    departure_entropy_j_per_kmol_k: float


@dataclass(frozen=True, slots=True)
class SRKMixtureState:
    phase: str
    compressibility: float
    fugacity_coefficients: tuple[float, ...]
    density_kg_per_m3: float
    departure_enthalpy_j_per_kmol: float
    departure_entropy_j_per_kmol_k: float


@dataclass(frozen=True, slots=True)
class SRKTPFlashResult:
    report: SolverReport
    temperature_k: float
    pressure_pa: float
    phase: str
    vapor_fraction: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    equilibrium_ratios: tuple[float, ...]
    liquid_state: SRKMixtureState | None
    vapor_state: SRKMixtureState | None


def _validate_state(temperature_k: float, pressure_pa: float) -> None:
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
        for value in (temperature_k, pressure_pa)
    ):
        raise ValidationError("absolute temperature and pressure must be finite and positive")


def _physical_roots(a_reduced: float, b_reduced: float, label: str) -> tuple[float, ...]:
    try:
        roots = _cubic_real_roots(
            -1.0,
            a_reduced - b_reduced - b_reduced**2,
            -a_reduced * b_reduced,
        )
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError(f"{label} cubic is outside the representable range") from error
    physical = tuple(root for root in roots if math.isfinite(root) and root > b_reduced)
    if not physical:
        raise ValidationError(f"{label} EOS has no physical compressibility root")
    return physical


class SoaveRedlichKwong:
    def __init__(self, compound: Compound):
        if not isinstance(compound, Compound):
            raise ValidationError("SRK compound data is invalid")
        self.compound = compound

    def parameters(self, temperature_k: float) -> SRKParameters:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        critical_temperature = self.compound.critical_temperature.value
        critical_pressure = self.compound.critical_pressure.value
        acentric_factor = self.compound.acentric_factor.value
        try:
            kappa = 0.48 + 1.574 * acentric_factor - 0.176 * acentric_factor**2
            alpha_base = 1.0 + kappa * (
                1.0 - math.sqrt(temperature_k / critical_temperature)
            )
            alpha = alpha_base**2
            a0 = 0.42748 * R**2 * critical_temperature**2 / critical_pressure
            parameters = SRKParameters(
                kappa,
                alpha,
                a0 * alpha,
                0.08664 * R * critical_temperature / critical_pressure,
                -a0
                * kappa
                * alpha_base
                / math.sqrt(temperature_k * critical_temperature),
            )
        except (OverflowError, ZeroDivisionError) as error:
            raise ValidationError("SRK parameters are outside the representable range") from error
        if (
            not all(
                math.isfinite(value) and value > 0
                for value in (
                    parameters.alpha,
                    parameters.a_pa_m6_per_kmol2,
                    parameters.b_m3_per_kmol,
                )
            )
            or not math.isfinite(parameters.da_dtemperature)
        ):
            raise ValidationError("SRK parameters are outside the representable range")
        return parameters

    def _reduced_parameters(
        self, temperature_k: float, pressure_pa: float
    ) -> tuple[SRKParameters, float, float]:
        _validate_state(temperature_k, pressure_pa)
        parameters = self.parameters(temperature_k)
        a_reduced, b_reduced = _reduced(
            parameters.a_pa_m6_per_kmol2,
            parameters.b_m3_per_kmol,
            temperature_k,
            pressure_pa,
        )
        return parameters, a_reduced, b_reduced

    def roots(self, temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
        _, a_reduced, b_reduced = self._reduced_parameters(temperature_k, pressure_pa)
        return _physical_roots(a_reduced, b_reduced, "SRK")

    def state(self, temperature_k: float, pressure_pa: float, phase: str) -> SRKState:
        if phase not in {"vapor", "liquid"}:
            raise ValidationError("phase must be vapor or liquid")
        parameters, a_reduced, b_reduced = self._reduced_parameters(
            temperature_k, pressure_pa
        )
        roots = self.roots(temperature_k, pressure_pa)
        compressibility = max(roots) if phase == "vapor" else min(roots)
        log_ratio = math.log((compressibility + b_reduced) / compressibility)
        log_fugacity = (
            compressibility
            - 1.0
            - math.log(compressibility - b_reduced)
            - a_reduced * log_ratio / b_reduced
        )
        departure_enthalpy = R * temperature_k * (compressibility - 1.0) + (
            temperature_k * parameters.da_dtemperature
            - parameters.a_pa_m6_per_kmol2
        ) * log_ratio / parameters.b_m3_per_kmol
        departure_entropy = R * math.log(compressibility - b_reduced) + (
            parameters.da_dtemperature * log_ratio / parameters.b_m3_per_kmol
        )
        density = _positive_exp(
            math.log(pressure_pa)
            + math.log(self.compound.molecular_weight.value)
            - math.log(compressibility)
            - math.log(R)
            - math.log(temperature_k),
            "SRK density",
        )
        if not all(
            math.isfinite(value)
            for value in (log_fugacity, departure_enthalpy, departure_entropy, density)
        ):
            raise ValidationError("SRK state is outside the representable range")
        return SRKState(
            phase,
            compressibility,
            _positive_exp(log_fugacity, "SRK fugacity coefficient"),
            density,
            departure_enthalpy,
            departure_entropy,
        )

    def stable_state(self, temperature_k: float, pressure_pa: float) -> SRKState:
        if len(self.roots(temperature_k, pressure_pa)) == 1:
            return replace(self.state(temperature_k, pressure_pa, "vapor"), phase="single")
        return min(
            (self.state(temperature_k, pressure_pa, phase) for phase in ("liquid", "vapor")),
            key=lambda state: state.fugacity_coefficient,
        )


class SoaveRedlichKwongMixture:
    def __init__(
        self,
        compounds: tuple[Compound, ...],
        mole_fractions: tuple[float, ...],
        interactions: PRInteractions,
    ):
        if len(compounds) < 2 or len(compounds) != len(mole_fractions):
            raise ValidationError(
                "SRK mixture compounds and mole fractions must have the same length of at least two"
            )
        if len({compound.id for compound in compounds}) != len(compounds):
            raise ValidationError("SRK mixture compounds must be unique")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
            for value in mole_fractions
        ):
            raise ValidationError("SRK mole fractions must be finite and non-negative")
        if not math.isclose(math.fsum(mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValidationError("SRK mole fractions must sum to one")
        if not isinstance(
            interactions, (PRInteractions, EvaluatedAdvancedCubicInteractions)
        ):
            raise ValidationError("SRK interaction data is invalid")
        self.components = tuple(SoaveRedlichKwong(compound) for compound in compounds)
        self.mole_fractions = tuple(mole_fractions)
        self.interactions = interactions

    def _mix(
        self, temperature_k: float
    ) -> tuple[
        SRKMixtureParameters,
        tuple[SRKParameters, ...],
        tuple[tuple[float, ...], ...],
    ]:
        pure = tuple(component.parameters(temperature_k) for component in self.components)
        cross = tuple(
            tuple(
                math.sqrt(first.a_pa_m6_per_kmol2)
                * math.sqrt(second.a_pa_m6_per_kmol2)
                * (
                    1.0
                    - self.interactions.get(
                        self.components[i].compound.id,
                        self.components[j].compound.id,
                    )
                )
                for j, second in enumerate(pure)
            )
            for i, first in enumerate(pure)
        )
        a_mixture = math.fsum(
            self.mole_fractions[i] * self.mole_fractions[j] * cross[i][j]
            for i in range(len(pure))
            for j in range(len(pure))
        )
        b_mixture = math.fsum(
            fraction * parameter.b_m3_per_kmol
            for fraction, parameter in zip(self.mole_fractions, pure)
        )
        da_dtemperature = math.fsum(
            self.mole_fractions[i]
            * self.mole_fractions[j]
            * cross[i][j]
            * 0.5
            * (
                pure[i].da_dtemperature / pure[i].a_pa_m6_per_kmol2
                + pure[j].da_dtemperature / pure[j].a_pa_m6_per_kmol2
            )
            for i in range(len(pure))
            for j in range(len(pure))
        )
        if (
            not all(math.isfinite(value) and value > 0 for value in (a_mixture, b_mixture))
            or not math.isfinite(da_dtemperature)
        ):
            raise ValidationError("SRK mixture parameters are outside the representable range")
        return SRKMixtureParameters(a_mixture, b_mixture, da_dtemperature), pure, cross

    def parameters(self, temperature_k: float) -> SRKMixtureParameters:
        return self._mix(temperature_k)[0]

    def roots(self, temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
        _validate_state(temperature_k, pressure_pa)
        parameters = self.parameters(temperature_k)
        a_reduced, b_reduced = _reduced(
            parameters.a_pa_m6_per_kmol2,
            parameters.b_m3_per_kmol,
            temperature_k,
            pressure_pa,
        )
        return _physical_roots(a_reduced, b_reduced, "SRK mixture")

    def state(
        self, temperature_k: float, pressure_pa: float, phase: str
    ) -> SRKMixtureState:
        if phase not in {"vapor", "liquid"}:
            raise ValidationError("phase must be vapor or liquid")
        _validate_state(temperature_k, pressure_pa)
        parameters, pure, cross = self._mix(temperature_k)
        a_reduced, b_reduced = _reduced(
            parameters.a_pa_m6_per_kmol2,
            parameters.b_m3_per_kmol,
            temperature_k,
            pressure_pa,
        )
        roots = self.roots(temperature_k, pressure_pa)
        compressibility = max(roots) if phase == "vapor" else min(roots)
        log_ratio = math.log((compressibility + b_reduced) / compressibility)
        log_fugacities = tuple(
            pure[i].b_m3_per_kmol
            / parameters.b_m3_per_kmol
            * (compressibility - 1.0)
            - math.log(compressibility - b_reduced)
            - a_reduced
            / b_reduced
            * (
                2.0
                * math.fsum(
                    self.mole_fractions[j] * cross[i][j] for j in range(len(pure))
                )
                / parameters.a_pa_m6_per_kmol2
                - pure[i].b_m3_per_kmol / parameters.b_m3_per_kmol
            )
            * log_ratio
            for i in range(len(pure))
        )
        departure_enthalpy = R * temperature_k * (compressibility - 1.0) + (
            temperature_k * parameters.da_dtemperature
            - parameters.a_pa_m6_per_kmol2
        ) * log_ratio / parameters.b_m3_per_kmol
        departure_entropy = R * math.log(compressibility - b_reduced) + (
            parameters.da_dtemperature * log_ratio / parameters.b_m3_per_kmol
        )
        molecular_weight = math.fsum(
            fraction * component.compound.molecular_weight.value
            for fraction, component in zip(self.mole_fractions, self.components)
        )
        density = _positive_exp(
            math.log(pressure_pa)
            + math.log(molecular_weight)
            - math.log(compressibility)
            - math.log(R)
            - math.log(temperature_k),
            "SRK mixture density",
        )
        if not all(
            math.isfinite(value)
            for value in (*log_fugacities, departure_enthalpy, departure_entropy, density)
        ):
            raise ValidationError("SRK mixture state is outside the representable range")
        return SRKMixtureState(
            phase,
            compressibility,
            tuple(_positive_exp(value, "SRK mixture fugacity coefficient") for value in log_fugacities),
            density,
            departure_enthalpy,
            departure_entropy,
        )

    def stable_state(self, temperature_k: float, pressure_pa: float) -> SRKMixtureState:
        if len(self.roots(temperature_k, pressure_pa)) == 1:
            return replace(self.state(temperature_k, pressure_pa, "vapor"), phase="single")
        return min(
            (self.state(temperature_k, pressure_pa, phase) for phase in ("liquid", "vapor")),
            key=lambda state: math.fsum(
                fraction * math.log(coefficient)
                for fraction, coefficient in zip(
                    self.mole_fractions, state.fugacity_coefficients
                )
            ),
        )


def srk_tp_flash(
    compounds: tuple[Compound, ...],
    composition: tuple[float, ...],
    interactions: PRInteractions,
    temperature_k: float,
    pressure_pa: float,
    *,
    max_iterations: int = 100,
    tolerance: float = 1.0e-10,
) -> SRKTPFlashResult:
    """Solve a classic-SRK VLE TP flash by fugacity-ratio iteration."""
    _validate_state(temperature_k, pressure_pa)
    if (
        isinstance(max_iterations, bool)
        or not isinstance(max_iterations, int)
        or max_iterations <= 0
        or isinstance(tolerance, bool)
        or not isinstance(tolerance, (int, float))
        or not math.isfinite(tolerance)
        or tolerance <= 0
    ):
        raise ValidationError("SRK flash limits must be finite and positive")
    # Constructor validates the ordered domain and composition before Wilson arithmetic.
    SoaveRedlichKwongMixture(compounds, composition, interactions)
    ratios = tuple(
        math.exp(
            math.log(compound.critical_pressure.value / pressure_pa)
            + 5.373
            * (1.0 + compound.acentric_factor.value)
            * (1.0 - compound.critical_temperature.value / temperature_k)
        )
        for compound in compounds
    )
    residual = math.inf
    liquid_state = None
    vapor_state = None
    flash = rachford_rice(composition, ratios, max_iterations=max_iterations)
    for iteration in range(1, max_iterations + 1):
        flash = rachford_rice(composition, ratios, max_iterations=max_iterations)
        if flash.phase != "two-phase":
            state = SoaveRedlichKwongMixture(
                compounds, composition, interactions
            ).stable_state(temperature_k, pressure_pa)
            report = SolverReport(
                True,
                iteration - 1,
                0.0,
                "SRK fugacity-ratio successive substitution",
                (),
                None,
            )
            return SRKTPFlashResult(
                report,
                temperature_k,
                pressure_pa,
                flash.phase,
                flash.vapor_fraction,
                flash.liquid_composition,
                flash.vapor_composition,
                ratios,
                state if flash.phase == "liquid" else None,
                state if flash.phase == "vapor" else None,
            )
        liquid_total = math.fsum(flash.liquid_composition)
        vapor_total = math.fsum(flash.vapor_composition)
        liquid_composition = tuple(value / liquid_total for value in flash.liquid_composition)
        vapor_composition = tuple(value / vapor_total for value in flash.vapor_composition)
        liquid_state = SoaveRedlichKwongMixture(
            compounds, liquid_composition, interactions
        ).state(temperature_k, pressure_pa, "liquid")
        vapor_state = SoaveRedlichKwongMixture(
            compounds, vapor_composition, interactions
        ).state(temperature_k, pressure_pa, "vapor")
        updated = tuple(
            liquid / vapor
            for liquid, vapor in zip(
                liquid_state.fugacity_coefficients,
                vapor_state.fugacity_coefficients,
            )
        )
        residual = max(abs(math.log(new / old)) for new, old in zip(updated, ratios))
        ratios = updated
        if residual <= tolerance:
            flash = rachford_rice(
                composition, ratios, max_iterations=max_iterations, tolerance=tolerance
            )
            liquid_total = math.fsum(flash.liquid_composition)
            vapor_total = math.fsum(flash.vapor_composition)
            liquid_composition = tuple(
                value / liquid_total for value in flash.liquid_composition
            )
            vapor_composition = tuple(
                value / vapor_total for value in flash.vapor_composition
            )
            liquid_state = SoaveRedlichKwongMixture(
                compounds, liquid_composition, interactions
            ).state(temperature_k, pressure_pa, "liquid")
            vapor_state = SoaveRedlichKwongMixture(
                compounds, vapor_composition, interactions
            ).state(temperature_k, pressure_pa, "vapor")
            report = SolverReport(
                True,
                iteration,
                residual,
                "SRK fugacity-ratio successive substitution",
                (),
                None,
            )
            return SRKTPFlashResult(
                report,
                temperature_k,
                pressure_pa,
                "two-phase",
                flash.vapor_fraction,
                liquid_composition,
                vapor_composition,
                ratios,
                liquid_state,
                vapor_state,
            )
    report = SolverReport(
        False,
        max_iterations,
        residual,
        "SRK fugacity-ratio successive substitution",
        (),
        "SRK TP flash iteration limit reached",
    )
    return SRKTPFlashResult(
        report,
        temperature_k,
        pressure_pa,
        flash.phase,
        flash.vapor_fraction,
        flash.liquid_composition,
        flash.vapor_composition,
        ratios,
        liquid_state,
        vapor_state,
    )
