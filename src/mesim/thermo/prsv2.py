"""DWSIM-compatible Peng-Robinson-Stryjek-Vera 2 equilibrium models."""

from dataclasses import dataclass, replace
import json
import math
from pathlib import Path

from ..compounds import Compound
from ..errors import ValidationError
from .peng_robinson import SQRT_2, _physical_roots, _positive_exp


DWSIM_PRSV2_GAS_CONSTANT = 8314.0  # J/kmol/K; mirrors 8.314 in DWSIM's model
_MIXING_RULES = frozenset({"margules", "van-laar"})


@dataclass(frozen=True, slots=True)
class PRSV2AlphaParameters:
    compound: str
    kappa1: float
    kappa2: float
    kappa3: float


@dataclass(frozen=True, slots=True)
class PRSV2InteractionParameters:
    first: str
    second: str
    k12: float
    k21: float
    reference_temperature_k: float


@dataclass(frozen=True, slots=True)
class PRSV2Data:
    source_revision: str
    alpha_parameters: tuple[PRSV2AlphaParameters, ...]
    margules_interactions: tuple[PRSV2InteractionParameters, ...]
    van_laar_interactions: tuple[PRSV2InteractionParameters, ...]

    def alpha(self, compound: str) -> PRSV2AlphaParameters:
        """Return an exact-case source record or DWSIM's zero-record fallback."""
        if not isinstance(compound, str) or not compound:
            raise ValidationError("PRSV2 compound lookup key must be non-empty")
        return next(
            (record for record in self.alpha_parameters if record.compound == compound),
            PRSV2AlphaParameters(compound, 0.0, 0.0, 0.0),
        )

    def interaction(
        self, first: str, second: str, mixing_rule: str
    ) -> PRSV2InteractionParameters | None:
        if mixing_rule not in _MIXING_RULES:
            raise ValidationError(f"unsupported PRSV2 mixing rule: {mixing_rule}")
        records = (
            self.margules_interactions
            if mixing_rule == "margules"
            else self.van_laar_interactions
        )
        return next(
            (
                record
                for record in records
                if (record.first == first and record.second == second)
                or (record.first == second and record.second == first)
            ),
            None,
        )


@dataclass(frozen=True, slots=True)
class PRSV2ComponentParameters:
    correction: float
    alpha: float
    a_pa_m6_per_kmol2: float
    b_m3_per_kmol: float


@dataclass(frozen=True, slots=True)
class PRSV2MixtureState:
    phase: str
    compressibility: float
    fugacity_coefficients: tuple[float, ...]


def _finite_number(value, label: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise ValidationError(f"{label} must be finite")
    return float(value)


def load_prsv2_data(path: str | Path) -> PRSV2Data:
    """Load the frozen DWSIM PRSV2 alpha and asymmetric interaction tables."""
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if document["schema_version"] != "dwsim-prsv2-data-1":
            raise ValidationError("unsupported PRSV2 data schema")
        source = document["source"]
        alpha = tuple(
            PRSV2AlphaParameters(
                record["compound"],
                _finite_number(record["kappa1"], "PRSV2 kappa1"),
                _finite_number(record["kappa2"], "PRSV2 kappa2"),
                _finite_number(record["kappa3"], "PRSV2 kappa3"),
            )
            for record in document["alpha_parameters"]
        )

        def interactions(name: str) -> tuple[PRSV2InteractionParameters, ...]:
            return tuple(
                PRSV2InteractionParameters(
                    record["first"],
                    record["second"],
                    _finite_number(record["k12"], "PRSV2 k12"),
                    _finite_number(record["k21"], "PRSV2 k21"),
                    _finite_number(
                        record["reference_temperature_k"],
                        "PRSV2 reference temperature",
                    ),
                )
                for record in document["interactions"][name]
            )

        data = PRSV2Data(
            source["revision"],
            alpha,
            interactions("margules"),
            interactions("van_laar"),
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid PRSV2 data: {error}") from error

    if (
        not isinstance(data.source_revision, str)
        or not data.source_revision
        or len(data.alpha_parameters) != 90
        or len(data.margules_interactions) != 8
        or len(data.van_laar_interactions) != 8
    ):
        raise ValidationError("PRSV2 data has an invalid source identity or record count")
    alpha_names = tuple(record.compound for record in data.alpha_parameters)
    if (
        len(set(alpha_names)) != len(alpha_names)
        or any(not isinstance(name, str) or not name for name in alpha_names)
    ):
        raise ValidationError("PRSV2 alpha names must be non-empty and unique")
    for records in (data.margules_interactions, data.van_laar_interactions):
        pairs = tuple(frozenset((record.first, record.second)) for record in records)
        if (
            len(set(pairs)) != len(pairs)
            or any(
                not isinstance(record.first, str)
                or not record.first
                or not isinstance(record.second, str)
                or not record.second
                or record.first == record.second
                or record.reference_temperature_k <= 0.0
                for record in records
            )
        ):
            raise ValidationError("PRSV2 interaction records are invalid")
    return data


class PRSV2Mixture:
    """PRSV2 phase-fugacity model with DWSIM's composition-dependent mixing."""

    def __init__(
        self,
        compounds: tuple[Compound, ...],
        mole_fractions: tuple[float, ...],
        data: PRSV2Data,
        mixing_rule: str,
    ):
        try:
            compounds = tuple(compounds)
            mole_fractions = tuple(mole_fractions)
        except TypeError as error:
            raise ValidationError("PRSV2 mixture inputs must be sequences") from error
        if (
            len(compounds) < 2
            or len(compounds) != len(mole_fractions)
            or any(not isinstance(compound, Compound) for compound in compounds)
            or len({compound.id for compound in compounds}) != len(compounds)
            or not isinstance(data, PRSV2Data)
            or mixing_rule not in _MIXING_RULES
        ):
            raise ValidationError("PRSV2 mixture domain is invalid")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0.0
            for value in mole_fractions
        ) or not math.isclose(
            math.fsum(mole_fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12
        ):
            raise ValidationError("PRSV2 mole fractions must be non-negative and sum to one")
        self.compounds = compounds
        self.mole_fractions = mole_fractions
        self.data = data
        self.mixing_rule = mixing_rule

    @staticmethod
    def _validate_state(temperature_k: float, pressure_pa: float) -> None:
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0.0
            for value in (temperature_k, pressure_pa)
        ):
            raise ValidationError("PRSV2 temperature and pressure must be positive")

    def component_parameters(
        self, compound: Compound, temperature_k: float
    ) -> PRSV2ComponentParameters:
        self._validate_state(temperature_k, 1.0)
        critical_temperature = compound.critical_temperature.value
        critical_pressure = compound.critical_pressure.value
        acentric_factor = compound.acentric_factor.value
        reduced_temperature = temperature_k / critical_temperature
        square_root = math.sqrt(reduced_temperature)
        record = self.data.alpha(compound.id)
        try:
            if record.kappa1 * record.kappa2 * record.kappa3 != 0.0:
                correction = (
                    0.378893
                    + 1.4897153 * acentric_factor
                    - 0.17131848 * acentric_factor**2
                    + 0.0196544 * acentric_factor**3
                    + (
                        record.kappa1
                        + record.kappa2
                        * (record.kappa3 - reduced_temperature)
                        * (1.0 - square_root)
                    )
                    * (1.0 + square_root)
                    * (0.7 - reduced_temperature)
                )
            elif acentric_factor <= 0.491:
                correction = (
                    0.37464
                    + 1.5422 * acentric_factor
                    - 0.26992 * acentric_factor**2
                )
            else:
                correction = (
                    0.379642
                    + 1.48503 * acentric_factor
                    - 0.164423 * acentric_factor**2
                    + 0.016666 * acentric_factor**3
                )
            alpha = (1.0 + correction * (1.0 - square_root)) ** 2
            a = (
                0.45724
                * alpha
                * DWSIM_PRSV2_GAS_CONSTANT**2
                * critical_temperature**2
                / critical_pressure
            )
            b = (
                0.0778
                * DWSIM_PRSV2_GAS_CONSTANT
                * critical_temperature
                / critical_pressure
            )
        except (OverflowError, ZeroDivisionError) as error:
            raise ValidationError("PRSV2 component parameters are unrepresentable") from error
        if not all(math.isfinite(value) for value in (correction, alpha, a, b)) or any(
            value <= 0.0 for value in (alpha, a, b)
        ):
            raise ValidationError("PRSV2 component parameters are unrepresentable")
        return PRSV2ComponentParameters(correction, alpha, a, b)

    def _interaction_matrices(
        self,
    ) -> tuple[tuple[tuple[float, ...], ...], tuple[tuple[float, ...], ...]]:
        size = len(self.compounds)
        k12 = [[0.0 for _ in range(size)] for _ in range(size)]
        k21 = [[0.0 for _ in range(size)] for _ in range(size)]
        for i, first in enumerate(self.compounds):
            for j, second in enumerate(self.compounds):
                if i == j:
                    continue
                record = self.data.interaction(first.id, second.id, self.mixing_rule)
                if record is not None:
                    # DWSIM RET_KIJ/RET_KIJ2 preserve the stored orientation even
                    # when the compound lookup is reversed.
                    k12[i][j] = record.k12
                    k21[i][j] = record.k21
        return tuple(map(tuple, k12)), tuple(map(tuple, k21))

    def _mix(
        self, temperature_k: float
    ) -> tuple[
        tuple[PRSV2ComponentParameters, ...],
        tuple[tuple[float, ...], ...],
        tuple[float, ...],
        float,
        float,
    ]:
        pure = tuple(
            self.component_parameters(compound, temperature_k)
            for compound in self.compounds
        )
        k12, k21 = self._interaction_matrices()
        size = len(pure)
        composition = self.mole_fractions
        cross: list[list[float]] = [[0.0 for _ in range(size)] for _ in range(size)]
        for i in range(size):
            for j in range(size):
                geometric = math.sqrt(
                    pure[i].a_pa_m6_per_kmol2 * pure[j].a_pa_m6_per_kmol2
                )
                if self.mixing_rule == "margules":
                    factor = (
                        1.0
                        - composition[i] * k12[i][j]
                        - composition[j] * k21[j][i]
                    )
                else:
                    denominator = (
                        composition[i] * k12[i][j]
                        + composition[j] * k21[j][i]
                    )
                    factor = (
                        1.0 - k12[i][j] * k21[j][i] / denominator
                        if denominator != 0.0
                        else 1.0
                    )
                cross[i][j] = geometric * factor

        partial = []
        for i in range(size):
            sum1 = math.fsum(composition[j] * cross[i][j] for j in range(size))
            sum2 = 0.0
            sum3 = 0.0
            for j in range(size):
                if i != j:
                    geometric = math.sqrt(cross[i][i] * cross[j][j])
                    if self.mixing_rule == "margules":
                        sum2 += (
                            composition[i]
                            * composition[j]
                            * geometric
                            * (
                                composition[j] * k21[j][i]
                                - (1.0 - composition[i]) * k12[i][j]
                            )
                        )
                    else:
                        denominator = (
                            composition[i] * k12[i][j]
                            + composition[j] * k12[j][i]
                        )
                        if denominator != 0.0:
                            sum2 += (
                                composition[i]
                                * composition[j]
                                * geometric
                                * k12[i][j]
                                * k21[j][i]
                                * (
                                    (1.0 - composition[i]) * k12[i][j]
                                    - composition[j] * k21[j][i]
                                )
                                / denominator**2
                            )
                for k in range(size):
                    if i != j and k > j and k != i:
                        geometric = math.sqrt(cross[j][j] * cross[k][k])
                        if self.mixing_rule == "margules":
                            sum3 += (
                                composition[j]
                                * composition[k]
                                * geometric
                                * (
                                    composition[j] * k12[j][k]
                                    + composition[k] * k21[k][j]
                                )
                            )
                        else:
                            denominator = (
                                composition[j] * k12[j][k]
                                + composition[k] * k21[k][j]
                            )
                            if denominator != 0.0:
                                # The DWSIM source has sqrt(-a_jj*a_kk) here;
                                # for real positive EOS parameters that term is
                                # NaN and is intentionally excluded below.
                                sum3 = math.nan
            partial.append(sum1 + sum2 + (sum3 if math.isfinite(sum3) else 0.0))

        a_mixture = math.fsum(
            composition[i] * composition[j] * cross[i][j]
            for i in range(size)
            for j in range(size)
        )
        b_mixture = math.fsum(
            composition[i] * pure[i].b_m3_per_kmol for i in range(size)
        )
        values = (
            a_mixture,
            b_mixture,
            *partial,
            *(value for row in cross for value in row),
        )
        if not all(math.isfinite(value) for value in values) or any(
            value <= 0.0 for value in (a_mixture, b_mixture)
        ):
            raise ValidationError("PRSV2 mixture parameters are unrepresentable")
        return pure, tuple(map(tuple, cross)), tuple(partial), a_mixture, b_mixture

    def roots(self, temperature_k: float, pressure_pa: float) -> tuple[float, ...]:
        self._validate_state(temperature_k, pressure_pa)
        _, _, _, a_mixture, b_mixture = self._mix(temperature_k)
        a_reduced = (
            a_mixture
            * pressure_pa
            / (DWSIM_PRSV2_GAS_CONSTANT * temperature_k) ** 2
        )
        b_reduced = (
            b_mixture
            * pressure_pa
            / (DWSIM_PRSV2_GAS_CONSTANT * temperature_k)
        )
        return _physical_roots(a_reduced, b_reduced, "PRSV2 mixture")

    def state(
        self, temperature_k: float, pressure_pa: float, phase: str
    ) -> PRSV2MixtureState:
        if phase not in {"liquid", "vapor"}:
            raise ValidationError("PRSV2 phase must be liquid or vapor")
        self._validate_state(temperature_k, pressure_pa)
        pure, _, partial, a_mixture, b_mixture = self._mix(temperature_k)
        a_reduced = (
            a_mixture
            * pressure_pa
            / (DWSIM_PRSV2_GAS_CONSTANT * temperature_k) ** 2
        )
        b_reduced = (
            b_mixture
            * pressure_pa
            / (DWSIM_PRSV2_GAS_CONSTANT * temperature_k)
        )
        roots = _physical_roots(a_reduced, b_reduced, "PRSV2 mixture")
        compressibility = min(roots) if phase == "liquid" else max(roots)
        try:
            log_ratio = math.log(
                (compressibility + (1.0 + SQRT_2) * b_reduced)
                / (compressibility + (1.0 - SQRT_2) * b_reduced)
            )
            log_fugacities = tuple(
                pure[i].b_m3_per_kmol
                * (compressibility - 1.0)
                / b_mixture
                - math.log(compressibility - b_reduced)
                - a_reduced
                * (
                    2.0 * partial[i] / a_mixture
                    - pure[i].b_m3_per_kmol / b_mixture
                )
                * log_ratio
                / (2.0 * SQRT_2 * b_reduced)
                for i in range(len(pure))
            )
        except (OverflowError, ValueError, ZeroDivisionError) as error:
            raise ValidationError("PRSV2 state is unrepresentable") from error
        if not all(math.isfinite(value) for value in log_fugacities):
            raise ValidationError("PRSV2 state is unrepresentable")
        return PRSV2MixtureState(
            phase,
            compressibility,
            tuple(
                _positive_exp(value, "PRSV2 fugacity coefficient")
                for value in log_fugacities
            ),
        )

    def stable_state(
        self, temperature_k: float, pressure_pa: float
    ) -> PRSV2MixtureState:
        if len(self.roots(temperature_k, pressure_pa)) == 1:
            return replace(
                self.state(temperature_k, pressure_pa, "vapor"), phase="single"
            )
        return min(
            (
                self.state(temperature_k, pressure_pa, phase)
                for phase in ("liquid", "vapor")
            ),
            key=lambda state: math.fsum(
                fraction * math.log(coefficient)
                for fraction, coefficient in zip(
                    self.mole_fractions, state.fugacity_coefficients
                )
            ),
        )


__all__ = (
    "DWSIM_PRSV2_GAS_CONSTANT",
    "PRSV2AlphaParameters",
    "PRSV2ComponentParameters",
    "PRSV2Data",
    "PRSV2InteractionParameters",
    "PRSV2Mixture",
    "PRSV2MixtureState",
    "load_prsv2_data",
)
