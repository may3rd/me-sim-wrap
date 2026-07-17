"""Versioned reaction definitions and thermochemistry."""
import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from .errors import ValidationError


@dataclass(frozen=True, slots=True)
class ReactionProvenance:
    source: str
    source_revision: str
    selection: str
    imported_utc: str


@dataclass(frozen=True, slots=True)
class CompoundThermochemistry:
    compound_id: str
    elements: tuple[tuple[str, float], ...]
    formation_temperature_k: float
    ideal_gas_formation_enthalpy_j_per_kmol: float
    ideal_gas_formation_gibbs_energy_j_per_kmol: float
    ideal_gas_formation_entropy_j_per_kmol_k: float


@dataclass(frozen=True, slots=True)
class ArrheniusRateDefinition:
    model: str
    pre_exponential_factor: float
    activation_energy_j_per_mol: float
    activation_energy_unit: str
    orders: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class KineticRateDefinition:
    basis: str
    concentration_unit: str
    rate_unit: str
    forward: ArrheniusRateDefinition
    reverse: ArrheniusRateDefinition


@dataclass(frozen=True, slots=True)
class ReactionDefinition:
    id: str
    name: str
    reaction_type: str
    base_reactant: str
    phase: str
    stoichiometry: tuple[tuple[str, float], ...]
    reaction_heat_j_per_kmol: float
    conversion_fraction: float | None
    equilibrium_constant_model: str | None
    reaction_basis: str | None
    kinetics: KineticRateDefinition | None


@dataclass(frozen=True, slots=True)
class ReactionData:
    thermochemistry: tuple[CompoundThermochemistry, ...]
    reactions: tuple[ReactionDefinition, ...]
    provenance: ReactionProvenance


def _property(record: dict, name: str, unit: str) -> float:
    value = record[name]["value"]
    if record[name]["unit"] != unit:
        raise ValidationError(f"{name} must use {unit}")
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValidationError(f"{name} must be finite numeric data")
    return float(value)


def _arrhenius_rate(record: dict) -> ArrheniusRateDefinition:
    activation_energy = record["activation_energy"]
    if activation_energy["unit"] != "J/mol":
        raise ValidationError("kinetic activation energy must use J/mol")
    return ArrheniusRateDefinition(
        record["model"],
        float(record["pre_exponential_factor"]),
        float(activation_energy["value"]),
        activation_energy["unit"],
        tuple((compound, float(order)) for compound, order in record["orders"].items()),
    )


def _kinetics(record: dict) -> KineticRateDefinition | None:
    kinetics = record.get("kinetics")
    if kinetics is None:
        return None
    return KineticRateDefinition(
        kinetics["basis"], kinetics["concentration_unit"], kinetics["rate_unit"],
        _arrhenius_rate(kinetics["forward"]), _arrhenius_rate(kinetics["reverse"]),
    )


def load_reaction_data(path: str | Path) -> ReactionData:
    try:
        raw = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if raw["schema_version"] != "reaction-data-1":
            raise ValidationError("unsupported reaction data schema")
        provenance = ReactionProvenance(**raw["provenance"])
        thermochemistry = tuple(
            CompoundThermochemistry(
                record["compound_id"],
                tuple((element, float(count)) for element, count in record["elements"].items()),
                _property(record, "formation_temperature", "K"),
                _property(record, "ideal_gas_formation_enthalpy", "J/kmol"),
                _property(record, "ideal_gas_formation_gibbs_energy", "J/kmol"),
                _property(record, "ideal_gas_formation_entropy", "J/kmol/K"),
            )
            for record in raw["thermochemistry"]
        )
        reactions = tuple(
            ReactionDefinition(
                record["id"], record["name"], record["reaction_type"], record["base_reactant"], record["phase"],
                tuple((compound, float(coefficient)) for compound, coefficient in record["stoichiometry"].items()),
                _property(record, "reaction_heat", "J/kmol"),
                float(record["conversion_percent"]) / 100.0 if "conversion_percent" in record else None,
                record.get("equilibrium_constant_model"),
                record.get("reaction_basis"),
                _kinetics(record),
            )
            for record in raw["reactions"]
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid reaction data: {error}") from error

    provenance_values = (
        provenance.source, provenance.source_revision,
        provenance.selection, provenance.imported_utc,
    )
    if not all(isinstance(value, str) and value for value in provenance_values):
        raise ValidationError("reaction provenance fields must be non-empty strings")
    try:
        imported = datetime.fromisoformat(provenance.imported_utc)
    except ValueError as error:
        raise ValidationError("reaction provenance imported_utc must be an ISO 8601 timestamp") from error
    if not provenance.imported_utc.endswith("Z") or imported.utcoffset() != timedelta(0):
        raise ValidationError("reaction provenance imported_utc must be UTC")
    if len({record.compound_id for record in thermochemistry}) != len(thermochemistry):
        raise ValidationError("reaction thermochemistry compound IDs must be unique")

    thermo_by_id = {record.compound_id: record for record in thermochemistry}
    for record in thermochemistry:
        if not record.compound_id or not record.elements:
            raise ValidationError("reaction thermochemistry requires a compound ID and explicit elements")
        if record.formation_temperature_k <= 0.0:
            raise ValidationError("formation temperature must be positive")
        calculated_formation_entropy = (
            record.ideal_gas_formation_enthalpy_j_per_kmol
            - record.ideal_gas_formation_gibbs_energy_j_per_kmol
        ) / record.formation_temperature_k
        if not math.isclose(
            record.ideal_gas_formation_entropy_j_per_kmol_k,
            calculated_formation_entropy,
            rel_tol=1.0e-12,
            abs_tol=1.0e-6,
        ):
            raise ValidationError(f"{record.compound_id} formation H/G/S values are inconsistent")
        for element, count in record.elements:
            if not element or not math.isfinite(count) or count <= 0.0 or not count.is_integer():
                raise ValidationError("element counts must be positive integers")

    if len({reaction.id for reaction in reactions}) != len(reactions):
        raise ValidationError("reaction IDs must be unique")
    for reaction in reactions:
        stoichiometry = dict(reaction.stoichiometry)
        if (
            not reaction.id or not reaction.name
            or reaction.reaction_type not in {"conversion", "equilibrium", "kinetic"}
            or reaction.phase not in {"mixture", "vapor", "liquid", "solid"}
        ):
            raise ValidationError("reaction identity and phase are invalid")
        if len(stoichiometry) < 2 or any(not math.isfinite(value) or value == 0.0 for value in stoichiometry.values()):
            raise ValidationError("reaction stoichiometry requires finite non-zero coefficients")
        if reaction.base_reactant not in stoichiometry or stoichiometry[reaction.base_reactant] >= 0.0:
            raise ValidationError("reaction base reactant must have a negative coefficient")
        if reaction.reaction_type == "conversion":
            if reaction.conversion_fraction is None or not 0.0 <= reaction.conversion_fraction <= 1.0:
                raise ValidationError("conversion reaction conversion must be between zero and one")
            if reaction.equilibrium_constant_model is not None or reaction.reaction_basis is not None:
                raise ValidationError("conversion reaction cannot define equilibrium settings")
            if reaction.kinetics is not None:
                raise ValidationError("conversion reaction cannot define kinetics")
        elif reaction.reaction_type == "equilibrium":
            if reaction.conversion_fraction is not None:
                raise ValidationError("equilibrium reaction cannot define a fixed conversion")
            if reaction.equilibrium_constant_model != "gibbs" or reaction.reaction_basis != "fugacity":
                raise ValidationError("equilibrium reaction must use Gibbs fugacity equilibrium")
            if reaction.phase != "vapor":
                raise ValidationError("only vapor equilibrium reactions are supported")
            if reaction.kinetics is not None:
                raise ValidationError("equilibrium reaction cannot define kinetics")
        else:
            if reaction.conversion_fraction is not None:
                raise ValidationError("kinetic reaction cannot define a fixed conversion")
            if reaction.equilibrium_constant_model is not None or reaction.reaction_basis is not None:
                raise ValidationError("kinetic reaction cannot define equilibrium settings")
            kinetics = reaction.kinetics
            if kinetics is None:
                raise ValidationError("kinetic reaction requires explicit kinetics")
            if (
                kinetics.basis != "molar_concentration"
                or kinetics.concentration_unit != "kmol/m3"
                or kinetics.rate_unit not in {"kmol/[m3.s]", "kmol/[m3.h]"}
            ):
                raise ValidationError("unsupported kinetic basis or original expression units")
            for direction in (kinetics.forward, kinetics.reverse):
                if (
                    direction.model != "arrhenius"
                    or direction.activation_energy_unit != "J/mol"
                    or not math.isfinite(direction.pre_exponential_factor)
                    or direction.pre_exponential_factor < 0.0
                    or not math.isfinite(direction.activation_energy_j_per_mol)
                    or direction.activation_energy_j_per_mol < 0.0
                    or any(
                        compound not in stoichiometry
                        or not math.isfinite(order)
                        or order < 0.0
                        for compound, order in direction.orders
                    )
                ):
                    raise ValidationError("invalid Arrhenius kinetic definition")
            if kinetics.forward.pre_exponential_factor <= 0.0:
                raise ValidationError("kinetic reaction requires a positive forward factor")
        if any(compound not in thermo_by_id for compound in stoichiometry):
            raise ValidationError("reaction stoichiometry is missing explicit thermochemistry")
        elements = {element for compound in stoichiometry for element, _ in thermo_by_id[compound].elements}
        for element in elements:
            balance = math.fsum(
                coefficient * dict(thermo_by_id[compound].elements).get(element, 0.0)
                for compound, coefficient in reaction.stoichiometry
            )
            if not math.isclose(balance, 0.0, rel_tol=0.0, abs_tol=1.0e-12):
                raise ValidationError(f"reaction {reaction.id} is not balanced for {element}")
        calculated_heat = math.fsum(
            coefficient * thermo_by_id[compound].ideal_gas_formation_enthalpy_j_per_kmol
            for compound, coefficient in reaction.stoichiometry
        )
        if not math.isclose(calculated_heat, reaction.reaction_heat_j_per_kmol, rel_tol=1.0e-12, abs_tol=1.0e-6):
            raise ValidationError(f"reaction {reaction.id} heat does not match formation enthalpies")
    return ReactionData(thermochemistry, reactions, provenance)
