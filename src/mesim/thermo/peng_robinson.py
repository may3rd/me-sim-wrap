import math
from dataclasses import dataclass

from ..compounds import Compound
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


class PengRobinson:
    def __init__(self, compound: Compound):
        self.compound = compound

    @staticmethod
    def _validate_state(temperature_k: float, pressure_pa: float) -> None:
        if not math.isfinite(temperature_k) or not math.isfinite(pressure_pa) or temperature_k <= 0 or pressure_pa <= 0:
            raise ValidationError("absolute temperature and pressure must be finite and positive")

    def parameters(self, temperature_k: float) -> PRParameters:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        critical_temperature = self.compound.critical_temperature.value
        critical_pressure = self.compound.critical_pressure.value
        acentric_factor = self.compound.acentric_factor.value
        kappa = 0.37464 + 1.54226 * acentric_factor - 0.26992 * acentric_factor**2
        alpha_base = 1.0 + kappa * (1.0 - math.sqrt(temperature_k / critical_temperature))
        alpha = alpha_base**2
        a0 = 0.45724 * R**2 * critical_temperature**2 / critical_pressure
        return PRParameters(
            kappa,
            alpha,
            a0 * alpha,
            0.07780 * R * critical_temperature / critical_pressure,
            -a0 * kappa * alpha_base / math.sqrt(temperature_k * critical_temperature),
        )

    def _reduced_parameters(self, temperature_k: float, pressure_pa: float) -> tuple[PRParameters, float, float]:
        self._validate_state(temperature_k, pressure_pa)
        parameters = self.parameters(temperature_k)
        a_reduced = parameters.a_pa_m6_per_kmol2 * pressure_pa / (R**2 * temperature_k**2)
        b_reduced = parameters.b_m3_per_kmol * pressure_pa / (R * temperature_k)
        if not math.isfinite(a_reduced) or not math.isfinite(b_reduced):
            raise ValidationError("Peng-Robinson reduced parameters are outside the representable range")
        return parameters, a_reduced, b_reduced

    def roots(self, temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
        _, a_reduced, b_reduced = self._reduced_parameters(temperature_k, pressure_pa)
        roots = _cubic_real_roots(
            -(1.0 - b_reduced),
            a_reduced - 3.0 * b_reduced**2 - 2.0 * b_reduced,
            -(a_reduced * b_reduced - b_reduced**2 - b_reduced**3),
        )
        physical = tuple(root for root in roots if math.isfinite(root) and root > b_reduced)
        if not physical:
            raise ValidationError("Peng-Robinson EOS has no physical compressibility root")
        return physical

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
        density = math.exp(
            math.log(pressure_pa)
            + math.log(self.compound.molecular_weight.value)
            - math.log(compressibility)
            - math.log(R)
            - math.log(temperature_k)
        )
        values = (log_fugacity, departure_enthalpy, departure_entropy, density)
        if not all(math.isfinite(value) for value in values):
            raise ValidationError("Peng-Robinson state is outside the representable range")
        return PRState(
            phase,
            compressibility,
            math.exp(log_fugacity),
            density,
            departure_enthalpy,
            departure_entropy,
        )
