import math
from dataclasses import dataclass, replace

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError


R = 8314.46261815324  # J/kmol/K
SQRT_2 = math.sqrt(2.0)


@dataclass(frozen=True)
class PRParameters:
    kappa: float
    alpha: float
    a_pa_m6_per_kmol2: float
    b_m3_per_kmol: float
    da_dtemperature: float


@dataclass(frozen=True)
class PRState:
    phase: str
    compressibility: float
    fugacity_coefficient: float
    density_kg_per_m3: float
    departure_enthalpy_j_per_kmol: float
    departure_entropy_j_per_kmol_k: float


@dataclass(frozen=True)
class PRMixtureParameters:
    a_pa_m6_per_kmol2: float
    b_m3_per_kmol: float
    da_dtemperature: float


@dataclass(frozen=True)
class PRMixtureState:
    phase: str
    compressibility: float
    fugacity_coefficients: tuple[float, ...]
    density_kg_per_m3: float
    departure_enthalpy_j_per_kmol: float
    departure_entropy_j_per_kmol_k: float


def _cubic_real_roots(c2: float, c1: float, c0: float) -> tuple[float, ...]:
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2**3 / 27.0 - c2 * c1 / 3.0 + c0
    discriminant = (q / 2.0) ** 2 + (p / 3.0) ** 3
    scale = max(1.0, abs((q / 2.0) ** 2), abs((p / 3.0) ** 3))
    if discriminant > 1e-14 * scale:
        root = math.copysign(abs(-q / 2.0 + math.sqrt(discriminant)) ** (1.0 / 3.0), -q / 2.0 + math.sqrt(discriminant))
        root += math.copysign(abs(-q / 2.0 - math.sqrt(discriminant)) ** (1.0 / 3.0), -q / 2.0 - math.sqrt(discriminant))
        return (root - c2 / 3.0,)
    if abs(p) < 1e-15:
        return (-c2 / 3.0,)
    argument = max(-1.0, min(1.0, -q / (2.0 * math.sqrt(-(p / 3.0) ** 3))))
    angle = math.acos(argument)
    roots = sorted(2.0 * math.sqrt(-p / 3.0) * math.cos((angle + 2.0 * index * math.pi) / 3.0) - c2 / 3.0 for index in range(3))
    return tuple(root for index, root in enumerate(roots) if index == 0 or not math.isclose(root, roots[index - 1], abs_tol=1e-12))


def _reduced(a: float, b: float, temperature_k: float, pressure_pa: float) -> tuple[float, float]:
    try:
        a_reduced = math.exp(math.log(a) + math.log(pressure_pa) - 2.0 * math.log(R) - 2.0 * math.log(temperature_k))
        b_reduced = math.exp(math.log(b) + math.log(pressure_pa) - math.log(R) - math.log(temperature_k))
    except (OverflowError, ValueError) as error:
        raise ValidationError("Peng-Robinson reduced parameters are outside the representable range") from error
    if not math.isfinite(a_reduced) or not math.isfinite(b_reduced) or a_reduced <= 0 or b_reduced <= 0:
        raise ValidationError("Peng-Robinson reduced parameters are outside the representable range")
    return a_reduced, b_reduced


def _positive_exp(value: float, label: str) -> float:
    try:
        result = math.exp(value)
    except OverflowError as error:
        raise ValidationError(f"{label} is outside the representable range") from error
    if not math.isfinite(result) or result <= 0:
        raise ValidationError(f"{label} is outside the representable range")
    return result


def _physical_roots(a_reduced: float, b_reduced: float, label: str) -> tuple[float, ...]:
    try:
        roots = _cubic_real_roots(
            -(1.0 - b_reduced),
            a_reduced - 3.0 * b_reduced**2 - 2.0 * b_reduced,
            -(a_reduced * b_reduced - b_reduced**2 - b_reduced**3),
        )
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError(f"{label} cubic is outside the representable range") from error
    physical = tuple(root for root in roots if math.isfinite(root) and root > b_reduced)
    if not physical:
        raise ValidationError(f"{label} EOS has no physical compressibility root")
    return physical


class PengRobinson:
    def __init__(self, compound: Compound):
        self.compound = compound

    @staticmethod
    def _validate_state(temperature_k: float, pressure_pa: float) -> None:
        if any(
            isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0
            for value in (temperature_k, pressure_pa)
        ):
            raise ValidationError("absolute temperature and pressure must be finite and positive")

    def parameters(self, temperature_k: float) -> PRParameters:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        critical_temperature = self.compound.critical_temperature.value
        critical_pressure = self.compound.critical_pressure.value
        acentric_factor = self.compound.acentric_factor.value
        try:
            kappa = 0.37464 + 1.54226 * acentric_factor - 0.26992 * acentric_factor**2
            alpha_base = 1.0 + kappa * (1.0 - math.sqrt(temperature_k / critical_temperature))
            alpha = alpha_base**2
            a0 = 0.45724 * R**2 * critical_temperature**2 / critical_pressure
            parameters = PRParameters(
                kappa,
                alpha,
                a0 * alpha,
                0.07780 * R * critical_temperature / critical_pressure,
                -a0 * kappa * alpha_base / math.sqrt(temperature_k * critical_temperature),
            )
        except (OverflowError, ZeroDivisionError) as error:
            raise ValidationError("Peng-Robinson parameters are outside the representable range") from error
        if not all(math.isfinite(value) and value > 0 for value in (parameters.alpha, parameters.a_pa_m6_per_kmol2, parameters.b_m3_per_kmol)) or not math.isfinite(parameters.da_dtemperature):
            raise ValidationError("Peng-Robinson parameters are outside the representable range")
        return parameters

    def _reduced_parameters(self, temperature_k: float, pressure_pa: float) -> tuple[PRParameters, float, float]:
        self._validate_state(temperature_k, pressure_pa)
        parameters = self.parameters(temperature_k)
        a_reduced, b_reduced = _reduced(parameters.a_pa_m6_per_kmol2, parameters.b_m3_per_kmol, temperature_k, pressure_pa)
        return parameters, a_reduced, b_reduced

    def roots(self, temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
        _, a_reduced, b_reduced = self._reduced_parameters(temperature_k, pressure_pa)
        return _physical_roots(a_reduced, b_reduced, "Peng-Robinson")

    def state(self, temperature_k: float, pressure_pa: float, phase: str) -> PRState:
        if phase not in {"vapor", "liquid"}:
            raise ValidationError("phase must be vapor or liquid")
        parameters, a_reduced, b_reduced = self._reduced_parameters(temperature_k, pressure_pa)
        roots = self.roots(temperature_k, pressure_pa)
        compressibility = max(roots) if phase == "vapor" else min(roots)
        log_ratio = math.log(
            (compressibility + (1.0 + SQRT_2) * b_reduced)
            / (compressibility + (1.0 - SQRT_2) * b_reduced)
        )
        log_fugacity = (
            compressibility
            - 1.0
            - math.log(compressibility - b_reduced)
            - a_reduced * log_ratio / (2.0 * SQRT_2 * b_reduced)
        )
        departure_enthalpy = R * temperature_k * (compressibility - 1.0) + (
            temperature_k * parameters.da_dtemperature - parameters.a_pa_m6_per_kmol2
        ) * log_ratio / (2.0 * SQRT_2 * parameters.b_m3_per_kmol)
        departure_entropy = R * math.log(compressibility - b_reduced) + (
            parameters.da_dtemperature * log_ratio / (2.0 * SQRT_2 * parameters.b_m3_per_kmol)
        )
        density = _positive_exp(
            math.log(pressure_pa)
            + math.log(self.compound.molecular_weight.value)
            - math.log(compressibility)
            - math.log(R)
            - math.log(temperature_k),
            "Peng-Robinson density",
        )
        values = (log_fugacity, departure_enthalpy, departure_entropy, density)
        if not all(math.isfinite(value) for value in values):
            raise ValidationError("Peng-Robinson state is outside the representable range")
        return PRState(
            phase,
            compressibility,
            _positive_exp(log_fugacity, "Peng-Robinson fugacity coefficient"),
            density,
            departure_enthalpy,
            departure_entropy,
        )

    def stable_state(self, temperature_k: float, pressure_pa: float) -> PRState:
        if len(self.roots(temperature_k, pressure_pa)) == 1:
            return replace(self.state(temperature_k, pressure_pa, "vapor"), phase="single")
        return min(
            (self.state(temperature_k, pressure_pa, phase) for phase in ("liquid", "vapor")),
            key=lambda state: state.fugacity_coefficient,
        )


class PengRobinsonMixture:
    def __init__(self, compounds: tuple[Compound, ...], mole_fractions: tuple[float, ...], interactions: PRInteractions):
        if len(compounds) < 2 or len(compounds) != len(mole_fractions):
            raise ValidationError("mixture compounds and mole fractions must have the same length of at least two")
        if len({compound.id for compound in compounds}) != len(compounds):
            raise ValidationError("mixture compounds must be unique")
        if not all(math.isfinite(value) and value >= 0 for value in mole_fractions):
            raise ValidationError("mole fractions must be finite and non-negative")
        if not math.isclose(sum(mole_fractions), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValidationError("mole fractions must sum to one")
        self.components = tuple(PengRobinson(compound) for compound in compounds)
        self.mole_fractions = mole_fractions
        self.interactions = interactions

    def _mix(self, temperature_k: float) -> tuple[PRMixtureParameters, tuple[PRParameters, ...], tuple[tuple[float, ...], ...]]:
        pure = tuple(component.parameters(temperature_k) for component in self.components)
        cross = tuple(
            tuple(
                math.sqrt(first.a_pa_m6_per_kmol2) * math.sqrt(second.a_pa_m6_per_kmol2)
                * (1.0 - self.interactions.get(self.components[i].compound.id, self.components[j].compound.id))
                for j, second in enumerate(pure)
            )
            for i, first in enumerate(pure)
        )
        a_mixture = sum(
            self.mole_fractions[i] * self.mole_fractions[j] * cross[i][j]
            for i in range(len(pure))
            for j in range(len(pure))
        )
        b_mixture = sum(fraction * parameter.b_m3_per_kmol for fraction, parameter in zip(self.mole_fractions, pure))
        da_dtemperature = sum(
            self.mole_fractions[i]
            * self.mole_fractions[j]
            * cross[i][j]
            * 0.5
            * (pure[i].da_dtemperature / pure[i].a_pa_m6_per_kmol2 + pure[j].da_dtemperature / pure[j].a_pa_m6_per_kmol2)
            for i in range(len(pure))
            for j in range(len(pure))
        )
        if not all(math.isfinite(value) and value > 0 for value in (a_mixture, b_mixture)) or not math.isfinite(da_dtemperature):
            raise ValidationError("Peng-Robinson mixture parameters are outside the representable range")
        return PRMixtureParameters(a_mixture, b_mixture, da_dtemperature), pure, cross

    def parameters(self, temperature_k: float) -> PRMixtureParameters:
        return self._mix(temperature_k)[0]

    def roots(self, temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
        PengRobinson._validate_state(temperature_k, pressure_pa)
        parameters = self.parameters(temperature_k)
        a_reduced, b_reduced = _reduced(parameters.a_pa_m6_per_kmol2, parameters.b_m3_per_kmol, temperature_k, pressure_pa)
        return _physical_roots(a_reduced, b_reduced, "Peng-Robinson mixture")

    def state(self, temperature_k: float, pressure_pa: float, phase: str) -> PRMixtureState:
        if phase not in {"vapor", "liquid"}:
            raise ValidationError("phase must be vapor or liquid")
        PengRobinson._validate_state(temperature_k, pressure_pa)
        parameters, pure, cross = self._mix(temperature_k)
        a_reduced, b_reduced = _reduced(parameters.a_pa_m6_per_kmol2, parameters.b_m3_per_kmol, temperature_k, pressure_pa)
        roots = self.roots(temperature_k, pressure_pa)
        compressibility = max(roots) if phase == "vapor" else min(roots)
        log_ratio = math.log(
            (compressibility + (1.0 + SQRT_2) * b_reduced)
            / (compressibility + (1.0 - SQRT_2) * b_reduced)
        )
        log_fugacities = tuple(
            pure[i].b_m3_per_kmol / parameters.b_m3_per_kmol * (compressibility - 1.0)
            - math.log(compressibility - b_reduced)
            - a_reduced
            / (2.0 * SQRT_2 * b_reduced)
            * (
                2.0 * sum(self.mole_fractions[j] * cross[i][j] for j in range(len(pure))) / parameters.a_pa_m6_per_kmol2
                - pure[i].b_m3_per_kmol / parameters.b_m3_per_kmol
            )
            * log_ratio
            for i in range(len(pure))
        )
        departure_enthalpy = R * temperature_k * (compressibility - 1.0) + (
            temperature_k * parameters.da_dtemperature - parameters.a_pa_m6_per_kmol2
        ) * log_ratio / (2.0 * SQRT_2 * parameters.b_m3_per_kmol)
        departure_entropy = R * math.log(compressibility - b_reduced) + (
            parameters.da_dtemperature * log_ratio / (2.0 * SQRT_2 * parameters.b_m3_per_kmol)
        )
        molecular_weight = sum(
            fraction * component.compound.molecular_weight.value
            for fraction, component in zip(self.mole_fractions, self.components)
        )
        density = _positive_exp(
            math.log(pressure_pa)
            + math.log(molecular_weight)
            - math.log(compressibility)
            - math.log(R)
            - math.log(temperature_k),
            "Peng-Robinson mixture density",
        )
        values = (*log_fugacities, departure_enthalpy, departure_entropy, density)
        if not all(math.isfinite(value) for value in values):
            raise ValidationError("Peng-Robinson mixture state is outside the representable range")
        fugacity_coefficients = tuple(_positive_exp(value, "Peng-Robinson mixture fugacity coefficient") for value in log_fugacities)
        return PRMixtureState(
            phase,
            compressibility,
            fugacity_coefficients,
            density,
            departure_enthalpy,
            departure_entropy,
        )

    def stable_state(self, temperature_k: float, pressure_pa: float) -> PRMixtureState:
        if len(self.roots(temperature_k, pressure_pa)) == 1:
            return replace(self.state(temperature_k, pressure_pa, "vapor"), phase="single")
        return min(
            (self.state(temperature_k, pressure_pa, phase) for phase in ("liquid", "vapor")),
            key=lambda state: sum(
                fraction * math.log(coefficient)
                for fraction, coefficient in zip(self.mole_fractions, state.fugacity_coefficients)
            ),
        )
