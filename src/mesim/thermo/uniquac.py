"""Scoped DWSIM UNIQUAC activity-coefficient data and equations."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import MissingCompoundData, ValidationError


@dataclass(frozen=True, slots=True)
class UniquacCompoundBasis:
    compound_id: str
    chemsep_id: str
    q: float
    r: float


@dataclass(frozen=True, slots=True)
class UniquacResolvedInteraction:
    first: str
    second: str
    a12: float
    a21: float
    b12: float
    b21: float
    c12: float
    c21: float

    def log_tau(self, first: str, second: str, temperature_k: float, gas_constant: float) -> float:
        if first == second:
            return 0.0
        if first == self.first and second == self.second:
            a, b, c = self.a12, self.b12, self.c12
        elif first == self.second and second == self.first:
            a, b, c = self.a21, self.b21, self.c21
        else:
            raise ValidationError(f"missing UNIQUAC interaction: {first}/{second}")
        return (-a + b * temperature_k + c * temperature_k**2) / (
            gas_constant * temperature_k
        )


@dataclass(frozen=True, slots=True)
class UniquacSourceInteraction:
    first_chemsep_id: int
    second_chemsep_id: int
    a12: float
    a21: float
    comment: str


@dataclass(frozen=True, slots=True)
class UniquacData:
    source_revision: str
    resource_sha256: str
    gas_constant_cal_per_mol_k: float
    compound_basis: tuple[UniquacCompoundBasis, ...]
    interaction: UniquacResolvedInteraction
    source_interactions: tuple[UniquacSourceInteraction, ...]

    def compound(self, compound_id: str) -> UniquacCompoundBasis:
        for record in self.compound_basis:
            if record.compound_id == compound_id:
                return record
        raise MissingCompoundData(f"missing UNIQUAC compound basis: {compound_id}")


def load_uniquac_data(path: str | Path) -> UniquacData:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if (
            document["schema_version"] != "dwsim-uniquac-data-1"
            or document["model"] != "UNIQUAC"
            or document["source"]["interaction_basis"] != "cal/mol"
        ):
            raise ValidationError("unsupported UNIQUAC data schema, model, or basis")
        source = document["source"]
        basis = tuple(
            UniquacCompoundBasis(
                record["compound_id"], record["chemsep_id"], record["q"], record["r"]
            )
            for record in document["compound_basis"]
        )
        pair = document["resolved_interaction"]
        interaction = UniquacResolvedInteraction(
            pair["first"], pair["second"], pair["A12"], pair["A21"],
            pair["B12"], pair["B21"], pair["C12"], pair["C21"],
        )
        source_interactions = tuple(
            UniquacSourceInteraction(
                record["first_chemsep_id"], record["second_chemsep_id"],
                record["A12"], record["A21"], record["comment"],
            )
            for record in document["source_interactions"]
        )
        data = UniquacData(
            source["revision"], source["resource_sha256"],
            source["gas_constant_cal_per_mol_k"], basis, interaction,
            source_interactions,
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid UNIQUAC data: {error}") from error
    if (
        not data.source_revision
        or len(data.resource_sha256) != 64
        or data.gas_constant_cal_per_mol_k != 1.98721
        or len(data.compound_basis) != 2
        or len(data.source_interactions) != 376
    ):
        raise ValidationError("UNIQUAC source identity or record count is invalid")
    if len({record.compound_id for record in data.compound_basis}) != 2:
        raise ValidationError("UNIQUAC compound IDs must be unique")
    for record in data.compound_basis:
        if (
            not record.compound_id or not record.chemsep_id
            or any(
                isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or value <= 0.0
                for value in (record.q, record.r)
            )
        ):
            raise ValidationError("invalid UNIQUAC compound basis")
    if {data.interaction.first, data.interaction.second} != {
        record.compound_id for record in data.compound_basis
    } or any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value)
        for value in (
            data.interaction.a12, data.interaction.a21,
            data.interaction.b12, data.interaction.b21,
            data.interaction.c12, data.interaction.c21,
        )
    ):
        raise ValidationError("invalid UNIQUAC resolved interaction")
    # The source asset deliberately contains alternative regressions for some
    # ChemSep ID pairs. DWSIM retains the first row, so source order and all
    # alternatives are frozen instead of enforcing artificial uniqueness.
    if any(
        record.first_chemsep_id <= 0 or record.second_chemsep_id <= 0
        or not math.isfinite(record.a12) or not math.isfinite(record.a21)
        for record in data.source_interactions
    ):
        raise ValidationError("invalid UNIQUAC source interaction")
    return data


def _logsumexp(values: tuple[float, ...]) -> float:
    maximum = max(values)
    return maximum + math.log(math.fsum(math.exp(value - maximum) for value in values))


def uniquac_activity_coefficients(
    data: UniquacData, compound_ids, composition, temperature_k: float
) -> tuple[float, ...]:
    if not isinstance(data, UniquacData):
        raise ValidationError("UNIQUAC data is required")
    try:
        ids, fractions = tuple(compound_ids), tuple(composition)
    except TypeError as error:
        raise ValidationError("UNIQUAC compound IDs and composition must be sequences") from error
    if (
        len(ids) < 2 or len(ids) != len(fractions) or len(set(ids)) != len(ids)
        or any(not isinstance(value, str) or not value for value in ids)
        or isinstance(temperature_k, bool) or not isinstance(temperature_k, (int, float))
        or not math.isfinite(temperature_k) or temperature_k <= 0.0
        or any(
            isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or value < 0.0 for value in fractions
        )
        or not math.isclose(math.fsum(fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    ):
        raise ValidationError("invalid UNIQUAC compound domain, composition, or temperature")
    basis = tuple(data.compound(compound_id) for compound_id in ids)
    size = len(ids)
    r_average = math.fsum(fractions[i] * basis[i].r for i in range(size))
    q_average = math.fsum(fractions[i] * basis[i].q for i in range(size))
    theta = tuple(fractions[i] * basis[i].q / q_average for i in range(size))
    log_tau = tuple(
        tuple(
            data.interaction.log_tau(ids[i], ids[j], temperature_k,
                                     data.gas_constant_cal_per_mol_k)
            for j in range(size)
        ) for i in range(size)
    )
    try:
        log_s = tuple(
            _logsumexp(tuple(
                math.log(theta[j]) + log_tau[j][i]
                for j in range(size) if theta[j] > 0.0
            )) for i in range(size)
        )
        coefficients = []
        for i in range(size):
            sum1 = math.fsum(
                math.exp(math.log(theta[j]) + log_tau[i][j] - log_s[j])
                for j in range(size) if theta[j] > 0.0
            )
            r_ratio = basis[i].r / r_average
            fi_over_theta = r_ratio / (basis[i].q / q_average)
            log_combinatorial = (
                1.0 - r_ratio + math.log(r_ratio)
                - 5.0 * basis[i].q
                * (1.0 - fi_over_theta + math.log(fi_over_theta))
            )
            log_residual = basis[i].q * (1.0 - log_s[i] - sum1)
            gamma = math.exp(log_combinatorial + log_residual)
            if not math.isfinite(gamma) or gamma <= 0.0:
                raise ValidationError("UNIQUAC activity coefficient is unrepresentable")
            coefficients.append(gamma)
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError("UNIQUAC state is outside the representable range") from error
    return tuple(coefficients)


__all__ = (
    "UniquacCompoundBasis", "UniquacData", "UniquacResolvedInteraction",
    "UniquacSourceInteraction", "load_uniquac_data", "uniquac_activity_coefficients",
)
