"""Original DWSIM UNIFAC group-contribution data and activity equation."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import MissingCompoundData, ValidationError


@dataclass(frozen=True, slots=True)
class UnifacGroup:
    primary_id: int
    secondary_id: int
    primary_name: str
    group_name: str
    r: float
    q: float


@dataclass(frozen=True, slots=True)
class UnifacInteraction:
    first_primary_id: int
    second_primary_id: int
    a_kelvin: float


@dataclass(frozen=True, slots=True)
class UnifacSurfaceFraction:
    secondary_id: int
    value: float


@dataclass(frozen=True, slots=True)
class UnifacCompoundBasis:
    compound_id: str
    q: float
    r: float
    group_surface_fractions: tuple[UnifacSurfaceFraction, ...]


@dataclass(frozen=True, slots=True)
class UnifacData:
    source_revision: str
    groups_sha256: str
    interactions_sha256: str
    compound_basis: tuple[UnifacCompoundBasis, ...]
    groups: tuple[UnifacGroup, ...]
    interactions: tuple[UnifacInteraction, ...]

    def compound(self, compound_id: str) -> UnifacCompoundBasis:
        for record in self.compound_basis:
            if record.compound_id == compound_id:
                return record
        raise MissingCompoundData(f"missing UNIFAC compound basis: {compound_id}")

    def group(self, secondary_id: int) -> UnifacGroup:
        for record in self.groups:
            if record.secondary_id == secondary_id:
                return record
        raise ValidationError(f"missing UNIFAC subgroup: {secondary_id}")

    def interaction(self, first_primary_id: int, second_primary_id: int) -> float:
        for record in self.interactions:
            if (
                record.first_primary_id == first_primary_id
                and record.second_primary_id == second_primary_id
            ):
                return record.a_kelvin
        # DWSIM's executable TAU method uses zero for an absent directed entry.
        return 0.0


def load_unifac_data(path: str | Path) -> UnifacData:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if document["schema_version"] != "dwsim-unifac-data-1" or document["model"] != "UNIFAC":
            raise ValidationError("unsupported UNIFAC data schema or model")
        source = document["source"]
        compounds = tuple(
            UnifacCompoundBasis(
                record["compound_id"], record["q"], record["r"],
                tuple(
                    UnifacSurfaceFraction(item["secondary_id"], item["value"])
                    for item in record["group_surface_fractions"]
                ),
            ) for record in document["compound_basis"]
        )
        groups = tuple(
            UnifacGroup(
                record["primary_id"], record["secondary_id"],
                record["primary_name"], record["group_name"], record["r"], record["q"],
            ) for record in document["groups"]
        )
        interactions = tuple(
            UnifacInteraction(
                record["first_primary_id"], record["second_primary_id"], record["a_kelvin"]
            ) for record in document["interactions"]
        )
        data = UnifacData(
            source["revision"], source["groups_sha256"], source["interactions_sha256"],
            compounds, groups, interactions,
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid UNIFAC data: {error}") from error
    if (
        not data.source_revision or len(data.groups_sha256) != 64
        or len(data.interactions_sha256) != 64 or len(data.compound_basis) != 2
        or len(data.groups) != 119 or len(data.interactions) != 1403
    ):
        raise ValidationError("UNIFAC source identity or record count is invalid")
    if (
        len({record.compound_id for record in data.compound_basis}) != 2
        or len({record.secondary_id for record in data.groups}) != len(data.groups)
        or len({(record.first_primary_id, record.second_primary_id) for record in data.interactions})
        != len(data.interactions)
    ):
        raise ValidationError("UNIFAC compound, subgroup, and interaction keys must be unique")
    for record in data.groups:
        if (
            record.primary_id <= 0 or record.secondary_id <= 0
            or not record.primary_name or not record.group_name
            or not math.isfinite(record.r) or record.r <= 0.0
            or not math.isfinite(record.q) or record.q < 0.0
        ):
            raise ValidationError("invalid UNIFAC group record")
    for record in data.interactions:
        if (
            record.first_primary_id <= 0 or record.second_primary_id <= 0
            or not math.isfinite(record.a_kelvin)
        ):
            raise ValidationError("invalid UNIFAC interaction record")
    for compound in data.compound_basis:
        if (
            not compound.compound_id or not math.isfinite(compound.q) or compound.q <= 0.0
            or not math.isfinite(compound.r) or compound.r <= 0.0
            or not compound.group_surface_fractions
            or len({item.secondary_id for item in compound.group_surface_fractions})
            != len(compound.group_surface_fractions)
            or any(
                not math.isfinite(item.value) or item.value <= 0.0
                for item in compound.group_surface_fractions
            )
            or not math.isclose(
                math.fsum(item.value for item in compound.group_surface_fractions),
                1.0, rel_tol=0.0, abs_tol=1.0e-12,
            )
        ):
            raise ValidationError("invalid UNIFAC compound group basis")
        for item in compound.group_surface_fractions:
            data.group(item.secondary_id)
    return data


def unifac_activity_coefficients(
    data: UnifacData, compound_ids, composition, temperature_k: float
) -> tuple[float, ...]:
    if not isinstance(data, UnifacData):
        raise ValidationError("UNIFAC data is required")
    try:
        ids, fractions = tuple(compound_ids), tuple(composition)
    except TypeError as error:
        raise ValidationError("UNIFAC compound IDs and composition must be sequences") from error
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
        raise ValidationError("invalid UNIFAC compound domain, composition, or temperature")
    basis = tuple(data.compound(compound_id) for compound_id in ids)
    subgroup_ids = tuple(dict.fromkeys(
        item.secondary_id for compound in basis for item in compound.group_surface_fractions
    ))
    eki = tuple(
        {item.secondary_id: item.value for item in compound.group_surface_fractions}
        for compound in basis
    )
    try:
        tau = {
            (first, second): math.exp(
                -data.interaction(
                    data.group(first).primary_id, data.group(second).primary_id
                ) / temperature_k
            )
            for first in subgroup_ids for second in subgroup_ids
        }
        beta = tuple({
            m: math.fsum(values.get(k, 0.0) * tau[(k, m)] for k in subgroup_ids)
            for m in subgroup_ids
        } for values in eki)
        sum_xq = math.fsum(fractions[i] * basis[i].q for i in range(len(basis)))
        theta = {
            k: math.fsum(
                fractions[i] * basis[i].q * eki[i].get(k, 0.0)
                for i in range(len(basis))
            ) / sum_xq for k in subgroup_ids
        }
        s = {
            k: math.fsum(theta[m] * tau[(m, k)] for m in subgroup_ids)
            for k in subgroup_ids
        }
        sum_xr = math.fsum(fractions[i] * basis[i].r for i in range(len(basis)))
        coefficients = []
        for i, compound in enumerate(basis):
            j_value = compound.r / sum_xr
            l_value = compound.q / sum_xq
            log_combinatorial = (
                1.0 - j_value + math.log(j_value)
                - 5.0 * compound.q
                * (1.0 - j_value / l_value + math.log(j_value / l_value))
            )
            residual_sum = math.fsum(
                theta[k] * beta[i][k] / s[k]
                - (eki[i][k] * math.log(beta[i][k] / s[k]) if k in eki[i] else 0.0)
                for k in subgroup_ids
            )
            gamma = math.exp(log_combinatorial + compound.q * (1.0 - residual_sum))
            if not math.isfinite(gamma) or gamma <= 0.0:
                raise ValidationError("UNIFAC activity coefficient is unrepresentable")
            coefficients.append(gamma)
    except (OverflowError, ValueError, ZeroDivisionError) as error:
        raise ValidationError("UNIFAC state is outside the representable range") from error
    return tuple(coefficients)


__all__ = (
    "UnifacCompoundBasis", "UnifacData", "UnifacGroup", "UnifacInteraction",
    "UnifacSurfaceFraction", "load_unifac_data", "unifac_activity_coefficients",
)
