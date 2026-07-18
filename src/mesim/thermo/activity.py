import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from ..errors import MissingCompoundData, ValidationError


CAL_TO_J = 4.184
NRTL_GAS_CONSTANT_J_PER_MOL_K = 1.98721 * CAL_TO_J


@dataclass(frozen=True)
class NRTLProvenance:
    source: str
    source_revision: str
    selection: str
    imported_utc: str


@dataclass(frozen=True)
class NRTLVaporPressure:
    compound_id: str
    coefficients: tuple[float, float, float, float, float]
    minimum_k: float
    maximum_k: float

    def evaluate(self, temperature_k: float) -> float:
        if not math.isfinite(temperature_k) or temperature_k <= 0:
            raise ValidationError("absolute temperature must be finite and positive")
        if not self.minimum_k <= temperature_k <= self.maximum_k:
            raise ValidationError(
                f"{self.compound_id} vapor pressure is outside {self.minimum_k:g}..{self.maximum_k:g} K"
            )
        a, b, c, d, e = self.coefficients
        try:
            value = math.exp(a + b / temperature_k + c * math.log(temperature_k) + d * temperature_k**e)
        except OverflowError as error:
            raise ValidationError("vapor pressure is outside the representable range") from error
        if not math.isfinite(value) or value <= 0:
            raise ValidationError("vapor pressure must be finite and positive")
        return value


@dataclass(frozen=True)
class NRTLInteraction:
    compound_1: str
    compound_2: str
    a12_j_per_mol: float
    a21_j_per_mol: float
    b12_j_per_mol_k: float
    b21_j_per_mol_k: float
    c12_j_per_mol_k2: float
    c21_j_per_mol_k2: float
    alpha12: float


@dataclass(frozen=True)
class NRTLVLEData:
    vapor_pressures: tuple[NRTLVaporPressure, ...]
    interactions: tuple[NRTLInteraction, ...]
    provenance: NRTLProvenance

    def vapor_pressure(self, compound_id: str) -> NRTLVaporPressure:
        for correlation in self.vapor_pressures:
            if correlation.compound_id == compound_id:
                return correlation
        raise MissingCompoundData(f"missing NRTL VLE vapor-pressure data: {compound_id}")

    def interaction(self, first: str, second: str) -> NRTLInteraction:
        for interaction in self.interactions:
            if interaction.compound_1 == first and interaction.compound_2 == second:
                return interaction
        for interaction in self.interactions:
            if interaction.compound_1 == second and interaction.compound_2 == first:
                return NRTLInteraction(
                    first,
                    second,
                    interaction.a21_j_per_mol,
                    interaction.a12_j_per_mol,
                    interaction.b21_j_per_mol_k,
                    interaction.b12_j_per_mol_k,
                    interaction.c21_j_per_mol_k2,
                    interaction.c12_j_per_mol_k2,
                    interaction.alpha12,
                )
        raise MissingCompoundData(f"missing NRTL interaction: {first}/{second}")


@dataclass(frozen=True)
class NRTLVLEResult:
    converged: bool
    iterations: int
    residual: float
    algorithm: str
    warnings: tuple[str, ...]
    failure_reason: str | None
    kind: str
    temperature_k: float
    pressure_pa: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    activity_coefficients: tuple[float, ...]


def _number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValidationError(f"{name} must be finite numeric data")
    return float(value)


def _quantity(record: dict, name: str, unit: str) -> float:
    if record["unit"] != unit:
        raise ValidationError(f"{name} must use {unit}")
    return _number(record["value"], name)


def load_nrtl_vle_data(path: str | Path) -> NRTLVLEData:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if document["schema_version"] != "nrtl-vle-data-1" or document["model"] != "NRTL":
            raise ValidationError("unsupported NRTL VLE data schema or model")
        provenance = NRTLProvenance(**document["provenance"])
        vapor_pressures = tuple(
            NRTLVaporPressure(
                record["compound_id"],
                tuple(_number(record[key], f"{record['compound_id']}.{key}") for key in "ABCDE"),
                _number(record["minimum_k"], f"{record['compound_id']}.minimum_k"),
                _number(record["maximum_k"], f"{record['compound_id']}.maximum_k"),
            )
            for record in document["vapor_pressure_correlations"]
            if record["equation"] == 101 and record["unit"] == "Pa"
        )
        if len(vapor_pressures) != len(document["vapor_pressure_correlations"]):
            raise ValidationError("NRTL VLE vapor pressures must use equation 101 in Pa")
        interactions = tuple(
            NRTLInteraction(
                record["compound_1"],
                record["compound_2"],
                _quantity(record["A12"], "A12", "cal/mol") * CAL_TO_J,
                _quantity(record["A21"], "A21", "cal/mol") * CAL_TO_J,
                _quantity(record["B12"], "B12", "cal/mol/K") * CAL_TO_J,
                _quantity(record["B21"], "B21", "cal/mol/K") * CAL_TO_J,
                _quantity(record["C12"], "C12", "cal/mol/K2") * CAL_TO_J,
                _quantity(record["C21"], "C21", "cal/mol/K2") * CAL_TO_J,
                _quantity(record["alpha12"], "alpha12", "dimensionless"),
            )
            for record in document["interactions"]
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid NRTL VLE data: {error}") from error

    if not all(isinstance(value, str) and value for value in vars(provenance).values()):
        raise ValidationError("NRTL VLE provenance fields must be non-empty strings")
    try:
        imported = datetime.fromisoformat(provenance.imported_utc)
    except ValueError as error:
        raise ValidationError("NRTL VLE imported_utc must be an ISO 8601 timestamp") from error
    if not provenance.imported_utc.endswith("Z") or imported.utcoffset() != timedelta(0):
        raise ValidationError("NRTL VLE imported_utc must be UTC")
    if len({record.compound_id for record in vapor_pressures}) != len(vapor_pressures):
        raise ValidationError("NRTL VLE vapor-pressure compound IDs must be unique")
    if len({(record.compound_1, record.compound_2) for record in interactions}) != len(interactions):
        raise ValidationError("directed NRTL interaction pairs must be unique")
    if not vapor_pressures or not interactions:
        raise ValidationError("NRTL VLE data must include correlations and interactions")
    for record in vapor_pressures:
        if not record.compound_id or record.minimum_k <= 0 or record.maximum_k <= record.minimum_k:
            raise ValidationError("invalid NRTL VLE vapor-pressure range")
    for record in interactions:
        if not record.compound_1 or not record.compound_2 or record.compound_1 == record.compound_2:
            raise ValidationError("NRTL interaction compound IDs must be non-empty and distinct")
        if record.alpha12 < 0:
            raise ValidationError("NRTL alpha12 must be nonnegative")
    return NRTLVLEData(vapor_pressures, interactions, provenance)


def _validate_state(compound_ids, composition, temperature_k: float) -> tuple[tuple[str, ...], tuple[float, ...]]:
    ids = tuple(compound_ids)
    fractions = tuple(composition)
    if not ids or len(ids) != len(fractions) or len(set(ids)) != len(ids):
        raise ValidationError("compound IDs and composition must be non-empty, aligned, and unique")
    if any(not isinstance(value, str) or not value for value in ids):
        raise ValidationError("compound IDs must be non-empty strings")
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 for value in fractions):
        raise ValidationError("composition must contain finite nonnegative numbers")
    if not math.isclose(math.fsum(fractions), 1.0, abs_tol=1e-12):
        raise ValidationError("composition must sum to one")
    if not math.isfinite(temperature_k) or temperature_k <= 0:
        raise ValidationError("absolute temperature must be finite and positive")
    return ids, tuple(float(value) for value in fractions)


def nrtl_activity_coefficients(
    data: NRTLVLEData,
    compound_ids,
    composition,
    temperature_k: float,
) -> tuple[float, ...]:
    if not isinstance(data, NRTLVLEData):
        raise ValidationError("NRTL VLE data is required")
    ids, fractions = _validate_state(compound_ids, composition, temperature_k)
    count = len(ids)
    tau = [[0.0] * count for _ in range(count)]
    reverse_tau = [[0.0] * count for _ in range(count)]
    alpha = [[0.0] * count for _ in range(count)]
    try:
        for i in range(count):
            data.vapor_pressure(ids[i])
            for j in range(count):
                if i == j:
                    continue
                interaction = data.interaction(ids[i], ids[j])
                denominator = NRTL_GAS_CONSTANT_J_PER_MOL_K * temperature_k
                tau[i][j] = (
                    interaction.a12_j_per_mol
                    + interaction.b12_j_per_mol_k * temperature_k
                    + interaction.c12_j_per_mol_k2 * temperature_k**2
                ) / denominator
                reverse_tau[i][j] = (
                    interaction.a21_j_per_mol
                    + interaction.b21_j_per_mol_k * temperature_k
                    + interaction.c21_j_per_mol_k2 * temperature_k**2
                ) / denominator
                alpha[i][j] = interaction.alpha12

        g = [[math.exp(-alpha[i][j] * tau[i][j]) for j in range(count)] for i in range(count)]
        reverse_g = [[math.exp(-alpha[i][j] * reverse_tau[i][j]) for j in range(count)] for i in range(count)]
    except (OverflowError, ZeroDivisionError) as error:
        raise ValidationError("NRTL state is outside the representable range") from error
    sums = [math.fsum(fractions[j] * reverse_g[i][j] for j in range(count)) for i in range(count)]
    weighted = [
        math.fsum(fractions[j] * reverse_g[i][j] * reverse_tau[i][j] for j in range(count))
        for i in range(count)
    ]
    if any(value <= 0 or not math.isfinite(value) for value in sums):
        raise ValidationError("NRTL local-composition denominator is not positive and finite")

    coefficients = []
    ratios = [weighted[j] / sums[j] for j in range(count)]
    for i in range(count):
        log_gamma = weighted[i] / sums[i] + math.fsum(
            fractions[j] * g[i][j] * (tau[i][j] - ratios[j]) / sums[j]
            for j in range(count)
        )
        try:
            gamma = math.exp(log_gamma)
        except OverflowError as error:
            raise ValidationError("NRTL activity coefficient is outside the representable range") from error
        if not math.isfinite(gamma) or gamma <= 0:
            raise ValidationError("NRTL activity coefficient must be finite and positive")
        coefficients.append(gamma)
    return tuple(coefficients)


def nrtl_bubble_pressure(
    data: NRTLVLEData,
    compound_ids,
    liquid_composition,
    temperature_k: float,
) -> NRTLVLEResult:
    ids, liquid = _validate_state(compound_ids, liquid_composition, temperature_k)
    gamma = nrtl_activity_coefficients(data, ids, liquid, temperature_k)
    saturation = tuple(data.vapor_pressure(compound).evaluate(temperature_k) for compound in ids)
    terms = tuple(liquid[i] * gamma[i] * saturation[i] for i in range(len(ids)))
    pressure = math.fsum(terms)
    if not math.isfinite(pressure) or pressure <= 0:
        raise ValidationError("NRTL bubble pressure must be finite and positive")
    vapor = tuple(value / pressure for value in terms)
    residual = abs(math.fsum(vapor) - 1.0)
    return NRTLVLEResult(True, 1, residual, "NRTL modified-Raoult bubble pressure", (), None, "bubble", temperature_k, pressure, liquid, vapor, gamma)


def nrtl_equilibrium_ratios(
    data: NRTLVLEData,
    compound_ids,
    liquid_composition,
    temperature_k: float,
    pressure_pa: float,
) -> tuple[float, ...]:
    ids, liquid = _validate_state(compound_ids, liquid_composition, temperature_k)
    if not math.isfinite(pressure_pa) or pressure_pa <= 0:
        raise ValidationError("absolute pressure must be finite and positive")
    gamma = nrtl_activity_coefficients(data, ids, liquid, temperature_k)
    ratios = tuple(
        gamma[index] * data.vapor_pressure(compound).evaluate(temperature_k) / pressure_pa
        for index, compound in enumerate(ids)
    )
    if any(not math.isfinite(value) or value <= 0 for value in ratios):
        raise ValidationError("NRTL equilibrium ratios must be finite and positive")
    return ratios


def nrtl_bubble_temperature(
    data: NRTLVLEData,
    compound_ids,
    liquid_composition,
    pressure_pa: float,
    bracket_k: tuple[float, float],
    max_iterations: int = 100,
    tolerance: float = 1e-10,
) -> NRTLVLEResult:
    try:
        lower_k, upper_k = bracket_k
    except (TypeError, ValueError) as error:
        raise ValidationError("temperature bracket must contain two values") from error
    if not all(math.isfinite(value) and value > 0 for value in (lower_k, upper_k, pressure_pa)):
        raise ValidationError("temperature bracket and pressure must be finite and positive")
    if lower_k >= upper_k:
        raise ValidationError("temperature bracket must be increasing")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValidationError("max_iterations must be a positive integer")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValidationError("tolerance must be finite and positive")
    ids, liquid = _validate_state(compound_ids, liquid_composition, lower_k)
    lower = nrtl_bubble_pressure(data, ids, liquid, lower_k)
    upper = nrtl_bubble_pressure(data, ids, liquid, upper_k)
    lower_residual = lower.pressure_pa - pressure_pa
    upper_residual = upper.pressure_pa - pressure_pa
    if lower_residual == 0:
        return _replace_nrtl_result(lower, kind="bubble-temperature", algorithm="NRTL bubble-temperature bisection")
    if upper_residual == 0:
        return _replace_nrtl_result(upper, kind="bubble-temperature", algorithm="NRTL bubble-temperature bisection")
    if lower_residual * upper_residual > 0:
        raise ValidationError("temperature bracket does not enclose an NRTL bubble point")

    current = lower
    scaled_residual = abs(lower_residual) / pressure_pa
    for iteration in range(1, max_iterations + 1):
        middle_k = (lower_k + upper_k) / 2.0
        current = nrtl_bubble_pressure(data, ids, liquid, middle_k)
        signed_residual = current.pressure_pa - pressure_pa
        scaled_residual = abs(signed_residual) / pressure_pa
        if scaled_residual <= tolerance:
            return NRTLVLEResult(
                True,
                iteration,
                scaled_residual,
                "NRTL bubble-temperature bisection",
                (),
                None,
                "bubble-temperature",
                current.temperature_k,
                current.pressure_pa,
                current.liquid_composition,
                current.vapor_composition,
                current.activity_coefficients,
            )
        if lower_residual * signed_residual <= 0:
            upper_k = middle_k
            upper_residual = signed_residual
        else:
            lower_k = middle_k
            lower_residual = signed_residual
    return NRTLVLEResult(
        False,
        max_iterations,
        scaled_residual,
        "NRTL bubble-temperature bisection",
        (),
        "maximum iterations exceeded",
        "bubble-temperature",
        current.temperature_k,
        current.pressure_pa,
        current.liquid_composition,
        current.vapor_composition,
        current.activity_coefficients,
    )


def _replace_nrtl_result(result: NRTLVLEResult, *, kind: str, algorithm: str) -> NRTLVLEResult:
    return NRTLVLEResult(
        result.converged,
        result.iterations,
        result.residual,
        algorithm,
        result.warnings,
        result.failure_reason,
        kind,
        result.temperature_k,
        result.pressure_pa,
        result.liquid_composition,
        result.vapor_composition,
        result.activity_coefficients,
    )


def nrtl_dew_pressure(
    data: NRTLVLEData,
    compound_ids,
    vapor_composition,
    temperature_k: float,
    max_iterations: int = 100,
    tolerance: float = 1e-12,
) -> NRTLVLEResult:
    ids, vapor = _validate_state(compound_ids, vapor_composition, temperature_k)
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations <= 0:
        raise ValidationError("max_iterations must be a positive integer")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValidationError("tolerance must be finite and positive")
    saturation = tuple(data.vapor_pressure(compound).evaluate(temperature_k) for compound in ids)
    liquid = vapor
    pressure = math.nan
    gamma = tuple(1.0 for _ in ids)
    residual = math.inf
    for iteration in range(1, max_iterations + 1):
        gamma = nrtl_activity_coefficients(data, ids, liquid, temperature_k)
        denominator = math.fsum(vapor[i] / (gamma[i] * saturation[i]) for i in range(len(ids)))
        if not math.isfinite(denominator) or denominator <= 0:
            raise ValidationError("NRTL dew-pressure denominator must be finite and positive")
        pressure = 1.0 / denominator
        updated = tuple(vapor[i] * pressure / (gamma[i] * saturation[i]) for i in range(len(ids)))
        total = math.fsum(updated)
        if not math.isfinite(total) or total <= 0:
            raise ValidationError("NRTL dew liquid composition is not finite and positive")
        updated = tuple(value / total for value in updated)
        residual = max(abs(updated[i] - liquid[i]) for i in range(len(ids)))
        if residual <= tolerance:
            final_gamma = nrtl_activity_coefficients(data, ids, updated, temperature_k)
            return NRTLVLEResult(True, iteration, residual, "NRTL modified-Raoult dew fixed point", (), None, "dew", temperature_k, pressure, updated, vapor, final_gamma)
        liquid = updated
    return NRTLVLEResult(False, max_iterations, residual, "NRTL modified-Raoult dew fixed point", (), "maximum iterations exceeded", "dew", temperature_k, pressure, liquid, vapor, gamma)
