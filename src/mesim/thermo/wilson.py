"""Scoped DWSIM Wilson activity-coefficient data and equations."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import MissingCompoundData, ValidationError


WILSON_GAS_CONSTANT_CAL_PER_MOL_K = 1.9872


@dataclass(frozen=True, slots=True)
class WilsonCompoundBasis:
    compound_id: str
    cas: str
    molar_volume_m3_per_kmol: float


@dataclass(frozen=True, slots=True)
class WilsonInteraction:
    first_cas: str
    second_cas: str
    a12_cal_per_mol: float
    a21_cal_per_mol: float


@dataclass(frozen=True, slots=True)
class WilsonData:
    source_revision: str
    resource_sha256: str
    compound_basis: tuple[WilsonCompoundBasis, ...]
    interactions: tuple[WilsonInteraction, ...]

    def compound(self, compound_id: str) -> WilsonCompoundBasis:
        for record in self.compound_basis:
            if record.compound_id == compound_id:
                return record
        raise MissingCompoundData(f"missing Wilson compound basis: {compound_id}")

    def energy(self, first_cas: str, second_cas: str) -> float:
        if first_cas == second_cas:
            return 0.0
        for record in self.interactions:
            if record.first_cas == first_cas and record.second_cas == second_cas:
                return record.a12_cal_per_mol
            if record.first_cas == second_cas and record.second_cas == first_cas:
                return record.a21_cal_per_mol
        raise ValidationError(f"missing Wilson interaction: {first_cas}/{second_cas}")


def load_wilson_data(path: str | Path) -> WilsonData:
    """Load the complete DWSIM 9.0.5 table and its scoped volume basis."""
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if (
            document["schema_version"] != "dwsim-wilson-data-1"
            or document["model"] != "Wilson"
        ):
            raise ValidationError("unsupported Wilson data schema or model")
        source = document["source"]
        if (
            source["interaction_basis"] != "cal/mol"
            or source["molar_volume_temperature_k"] != 298.15
        ):
            raise ValidationError("unsupported Wilson parameter basis")
        compound_basis = tuple(
            WilsonCompoundBasis(
                record["compound_id"],
                record["cas"],
                record["molar_volume_298_15_k"]["value"],
            )
            for record in document["compound_basis"]
            if record["molar_volume_298_15_k"]["unit"] == "m3/kmol"
        )
        interactions = tuple(
            WilsonInteraction(
                record["first_cas"],
                record["second_cas"],
                record["A12"]["value"],
                record["A21"]["value"],
            )
            for record in document["interactions"]
            if record["A12"]["unit"] == "cal/mol"
            and record["A21"]["unit"] == "cal/mol"
        )
        data = WilsonData(
            source["revision"],
            source["resource_sha256"],
            compound_basis,
            interactions,
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid Wilson data: {error}") from error
    if (
        not data.source_revision
        or len(data.resource_sha256) != 64
        or len(data.compound_basis) != len(document["compound_basis"])
        or len(data.interactions) != len(document["interactions"])
        or len(data.interactions) != 364
    ):
        raise ValidationError("Wilson source identity, unit, or record count is invalid")
    compound_ids = tuple(record.compound_id for record in data.compound_basis)
    cases = tuple(record.cas for record in data.compound_basis)
    pairs = tuple((record.first_cas, record.second_cas) for record in data.interactions)
    if (
        len(set(compound_ids)) != len(compound_ids)
        or len(set(cases)) != len(cases)
        or len(set(pairs)) != len(pairs)
    ):
        raise ValidationError("Wilson compound and interaction keys must be unique")
    for record in data.compound_basis:
        if (
            not record.compound_id
            or not record.cas
            or isinstance(record.molar_volume_m3_per_kmol, bool)
            or not isinstance(record.molar_volume_m3_per_kmol, (int, float))
            or not math.isfinite(record.molar_volume_m3_per_kmol)
            or record.molar_volume_m3_per_kmol <= 0.0
        ):
            raise ValidationError("invalid Wilson compound basis")
    for record in data.interactions:
        if (
            not record.first_cas
            or not record.second_cas
            or record.first_cas == record.second_cas
            or any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                for value in (record.a12_cal_per_mol, record.a21_cal_per_mol)
            )
        ):
            raise ValidationError("invalid Wilson interaction record")
    return data


def _logsumexp(values: tuple[float, ...]) -> float:
    maximum = max(values)
    return maximum + math.log(math.fsum(math.exp(value - maximum) for value in values))


def wilson_activity_coefficients(
    data: WilsonData,
    compound_ids,
    composition,
    temperature_k: float,
) -> tuple[float, ...]:
    """Reproduce DWSIM's Wilson equation on its exact cal/mol basis."""
    if not isinstance(data, WilsonData):
        raise ValidationError("Wilson data is required")
    try:
        ids = tuple(compound_ids)
        fractions = tuple(composition)
    except TypeError as error:
        raise ValidationError("Wilson compound IDs and composition must be sequences") from error
    if (
        len(ids) < 2
        or len(ids) != len(fractions)
        or len(set(ids)) != len(ids)
        or any(not isinstance(value, str) or not value for value in ids)
        or isinstance(temperature_k, bool)
        or not isinstance(temperature_k, (int, float))
        or not math.isfinite(temperature_k)
        or temperature_k <= 0.0
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0.0
            for value in fractions
        )
        or not math.isclose(math.fsum(fractions), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    ):
        raise ValidationError("invalid Wilson compound domain, composition, or temperature")

    basis = tuple(data.compound(compound_id) for compound_id in ids)
    size = len(ids)
    log_lambda = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(size):
            energy = data.energy(basis[i].cas, basis[j].cas)
            log_lambda[i][j] = (
                math.log(basis[j].molar_volume_m3_per_kmol)
                - math.log(basis[i].molar_volume_m3_per_kmol)
                - energy / (WILSON_GAS_CONSTANT_CAL_PER_MOL_K * temperature_k)
            )

    try:
        log_sums = tuple(
            _logsumexp(
                tuple(
                    math.log(fractions[j]) + log_lambda[i][j]
                    for j in range(size)
                    if fractions[j] > 0.0
                )
            )
            for i in range(size)
        )
        coefficients = []
        for i in range(size):
            sum2 = math.fsum(
                math.exp(math.log(fractions[k]) + log_lambda[k][i] - log_sums[k])
                for k in range(size)
                if fractions[k] > 0.0
            )
            gamma = math.exp(-log_sums[i] + 1.0 - sum2)
            if not math.isfinite(gamma) or gamma <= 0.0:
                raise ValidationError("Wilson activity coefficient is unrepresentable")
            coefficients.append(gamma)
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError("Wilson state is outside the representable range") from error
    return tuple(coefficients)


__all__ = (
    "WILSON_GAS_CONSTANT_CAL_PER_MOL_K",
    "WilsonCompoundBasis",
    "WilsonData",
    "WilsonInteraction",
    "load_wilson_data",
    "wilson_activity_coefficients",
)
