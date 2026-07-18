"""Peng-Robinson 1978 alpha-function variant used by DWSIM PR78."""

import math

from ..compounds import Compound, PRInteractions
from ..errors import ValidationError
from .peng_robinson import (
    R,
    PRMixtureState,
    PRParameters,
    PRState,
    PengRobinson,
    PengRobinsonMixture,
)


class PengRobinson1978(PengRobinson):
    """Classic PR cubic with the 1978 high-acentric-factor kappa branch."""

    def parameters(self, temperature_k: float) -> PRParameters:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        critical_temperature = self.compound.critical_temperature.value
        critical_pressure = self.compound.critical_pressure.value
        acentric_factor = self.compound.acentric_factor.value
        try:
            if acentric_factor <= 0.491:
                kappa = (
                    0.37464
                    + 1.5422 * acentric_factor
                    - 0.26992 * acentric_factor**2
                )
            else:
                kappa = (
                    0.379642
                    + 1.48503 * acentric_factor
                    - 0.164423 * acentric_factor**2
                    + 0.016666 * acentric_factor**3
                )
            alpha_base = 1.0 + kappa * (
                1.0 - math.sqrt(temperature_k / critical_temperature)
            )
            alpha = alpha_base**2
            a0 = 0.45724 * R**2 * critical_temperature**2 / critical_pressure
            parameters = PRParameters(
                kappa,
                alpha,
                a0 * alpha,
                0.07780 * R * critical_temperature / critical_pressure,
                -a0
                * kappa
                * alpha_base
                / math.sqrt(temperature_k * critical_temperature),
            )
        except (OverflowError, ZeroDivisionError) as error:
            raise ValidationError("PR78 parameters are outside the representable range") from error
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
            raise ValidationError("PR78 parameters are outside the representable range")
        return parameters

    def state(self, temperature_k: float, pressure_pa: float, phase: str) -> PRState:
        return super().state(temperature_k, pressure_pa, phase)


class PengRobinson1978Mixture(PengRobinsonMixture):
    def __init__(
        self,
        compounds: tuple[Compound, ...],
        mole_fractions: tuple[float, ...],
        interactions: PRInteractions,
    ):
        super().__init__(compounds, mole_fractions, interactions)
        self.components = tuple(PengRobinson1978(compound) for compound in compounds)

    def state(
        self, temperature_k: float, pressure_pa: float, phase: str
    ) -> PRMixtureState:
        return super().state(temperature_k, pressure_pa, phase)
