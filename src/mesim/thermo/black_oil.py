"""Scoped DWSIM Black Oil correlations and TP split."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import ValidationError
from .steam_tables import SteamTablesData, steam_saturation_pressure_pa


@dataclass(frozen=True, slots=True)
class BlackOilCompound:
    compound_id: str
    specific_gravity_oil: float
    specific_gravity_gas: float
    basic_sediment_water_percent: float
    gas_oil_ratio: float


@dataclass(frozen=True, slots=True)
class BlackOilData:
    source_revision: str
    case_sha256: str
    runtime_assembly_sha256: str
    property_package_source_sha256: str
    model_source_sha256: str
    flash_source_sha256: str
    petroleum_methods_source_sha256: str
    compounds: tuple[BlackOilCompound, ...]
    probe_temperature_k: float
    probe_pressure_pa: float
    probe_composition: tuple[float, ...]
    probe_vapor_pressures_pa: tuple[float, ...]
    probe_vaporized_fractions: tuple[float, ...]

    @property
    def compound_ids(self) -> tuple[str, ...]:
        return tuple(record.compound_id for record in self.compounds)

    def compound(self, compound_id: str) -> BlackOilCompound:
        try:
            return next(record for record in self.compounds if record.compound_id == compound_id)
        except StopIteration as error:
            raise ValidationError(f"compound is outside Black Oil domain: {compound_id}") from error


@dataclass(frozen=True, slots=True)
class BlackOilTPFlashResult:
    liquid_fraction: float
    vapor_fraction: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    equilibrium_ratios: tuple[float, ...]
    iterations: int


def load_black_oil_data(path: str | Path) -> BlackOilData:
    try:
        document = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        source = document["source"]
        probe = document["scoped_probe"]
        if document["schema_version"] != "dwsim-black-oil-data-1" or document["model"] != "Black Oil":
            raise ValidationError("unsupported Black Oil data schema or model")
        records = tuple(
            BlackOilCompound(
                value["compound_id"], value["specific_gravity_oil"],
                value["specific_gravity_gas"], value["basic_sediment_water_percent"],
                value["gas_oil_ratio"],
            )
            for value in document["compounds"]
        )
        data = BlackOilData(
            source["revision"], source["case_sha256"], source["runtime_assembly_sha256"],
            source["property_package_source_sha256"], source["model_source_sha256"],
            source["flash_source_sha256"], source["petroleum_methods_source_sha256"],
            records, probe["temperature_k"], probe["pressure_pa"],
            tuple(probe["composition"]), tuple(probe["component_vapor_pressures_pa"]),
            tuple(probe["component_vaporized_fractions"]),
        )
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid Black Oil data: {error}") from error
    hashes = (
        data.case_sha256, data.runtime_assembly_sha256,
        data.property_package_source_sha256, data.model_source_sha256,
        data.flash_source_sha256, data.petroleum_methods_source_sha256,
    )
    numeric = (
        data.probe_temperature_k, data.probe_pressure_pa, *data.probe_composition,
        *data.probe_vapor_pressures_pa, *data.probe_vaporized_fractions,
        *(value for record in data.compounds for value in (
            record.specific_gravity_oil, record.specific_gravity_gas,
            record.basic_sediment_water_percent, record.gas_oil_ratio,
        )),
    )
    if (
        not data.source_revision or len(data.compounds) != 2
        or data.compound_ids != ("n-Pentane", "n-Hexane")
        or len(set(data.compound_ids)) != len(data.compound_ids)
        or any(len(value) != 64 for value in hashes)
        or any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in numeric)
        or any(not 0.0 < record.specific_gravity_oil < 1.07 for record in data.compounds)
        or any(record.specific_gravity_gas <= 0.0 for record in data.compounds)
        or any(not 0.0 <= record.basic_sediment_water_percent < 100.0 for record in data.compounds)
        or any(record.gas_oil_ratio < 0.0 for record in data.compounds)
        or len(data.probe_composition) != 2 or len(data.probe_vapor_pressures_pa) != 2
        or len(data.probe_vaporized_fractions) != 2
    ):
        raise ValidationError("Black Oil source identity or scoped data are invalid")
    return data


def _composition(data: BlackOilData, compound_ids, composition) -> tuple[float, ...]:
    try:
        ids = tuple(compound_ids); values = tuple(composition)
    except TypeError as error:
        raise ValidationError("Black Oil IDs and composition must be sequences") from error
    if ids != data.compound_ids or len(values) != len(ids) or any(
        isinstance(value, bool) or not isinstance(value, (int, float))
        or not math.isfinite(value) or value < 0.0 for value in values
    ) or not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValidationError("Black Oil requires its exact normalized compound domain")
    return tuple(float(value) for value in values)


def _state(temperature_k: float, pressure_pa: float) -> None:
    if (
        isinstance(temperature_k, bool) or not isinstance(temperature_k, (int, float))
        or not math.isfinite(temperature_k) or temperature_k <= 0.0
        or isinstance(pressure_pa, bool) or not isinstance(pressure_pa, (int, float))
        or not math.isfinite(pressure_pa) or pressure_pa <= 0.0
    ):
        raise ValidationError("Black Oil temperature and pressure must be positive and finite")


def black_oil_vapor_pressure_pa(
    data: BlackOilData, steam_data: SteamTablesData, compound_id: str, temperature_k: float
) -> float:
    _state(temperature_k, 1.0)
    if not isinstance(steam_data, SteamTablesData):
        raise ValidationError("Black Oil requires frozen IAPWS region-4 data")
    record = data.compound(compound_id); sgo = record.specific_gravity_oil
    try:
        molecular_weight = ((math.log(1.07 - sgo) - 3.56073) / -2.93886) ** 10
        boiling_point = 1080.0 - math.exp(6.97996 - 0.01964 * molecular_weight ** (2.0 / 3.0))
        critical_temperature = 189.8 + 450.6 * sgo + (0.4244 + 0.1174 * sgo) * boiling_point + (0.1441 - 1.0069 * sgo) * 100000.0 / boiling_point
        critical_pressure = 100000.0 * math.exp(
            5.689 - 0.0566 / sgo
            - (0.43639 + 4.1216 / sgo + 0.21343 / sgo**2) * 0.001 * boiling_point
            + (0.47579 + 1.182 / sgo + 0.15302 / sgo**2) * 0.000001 * boiling_point**2
            - (2.4505 + 9.9099 / sgo**2) * 0.0000000001 * boiling_point**3
        )
        reduced_boiling_point = boiling_point / critical_temperature
        acentric_factor = (
            -math.log(critical_pressure / 101325.0) - 5.92714
            + 6.09648 / reduced_boiling_point + 1.28862 * math.log(reduced_boiling_point)
            - 0.169347 * reduced_boiling_point**6
        ) / (
            15.2518 - 15.6875 / reduced_boiling_point
            - 13.4721 * math.log(reduced_boiling_point)
            + 0.43577 * reduced_boiling_point**6
        )
        reduced_temperature = temperature_k / critical_temperature
        f0 = 5.92714 - 6.09648 / reduced_temperature - 1.28862 * math.log(reduced_temperature) + 0.169347 * reduced_temperature**6
        f1 = 15.2518 - 15.6875 / reduced_temperature - 13.4721 * math.log(reduced_temperature) + 0.43577 * reduced_temperature**6
        oil_pressure = critical_pressure * math.exp(f0 + acentric_factor * f1)
        water_numerical_bar = steam_saturation_pressure_pa(steam_data, temperature_k) / 100000.0
        bsw = record.basic_sediment_water_percent
        value = (100.0 - bsw) / 100.0 * oil_pressure + bsw / 100.0 * water_numerical_bar
    except (ValueError, OverflowError, ZeroDivisionError) as error:
        raise ValidationError("Black Oil vapor pressure is unrepresentable") from error
    if not math.isfinite(value) or value <= 0.0:
        raise ValidationError("Black Oil vapor pressure is unrepresentable")
    return value


def black_oil_component_vaporized_fraction(
    record: BlackOilCompound, temperature_k: float, pressure_pa: float
) -> float:
    _state(temperature_k, pressure_pa)
    sgo=record.specific_gravity_oil;sgg=record.specific_gravity_gas
    bsw=record.basic_sediment_water_percent;gor=record.gas_oil_ratio
    try:
        tf=(temperature_k-273.15)*9.0/5.0+32.0;trank=tf+459.67;ppsia=pressure_pa*0.000145038
        api=141.5/sgo-131.5;wor=bsw/(100.0-bsw);gorss=gor*5.6738
        rs=sgg*((ppsia/18.2+1.4)*10.0**(0.0125*api-0.00091*tf))**1.2048
        pb=18.2*((gorss/sgg)**(1.0/1.2048)*10.0**(0.00091*tf-0.0125*api)-1.4)
        bos=0.9759+0.00012*(rs*(sgg/sgo)**0.5+1.25*tf)**1.2
        sgfg100=sgg*(1.0+0.00005912*api*tf*math.log(ppsia/114.7)/math.log(10.0))
        c=0.0001*(2.81*gor+3.1*tf+171.0/sgfg100-118.0*sgfg100-1102.0)
        boss=bos*(pb/ppsia)**c;bo=bos if ppsia<pb else boss
        ppc=677.0+15.0*sgg-37.5*sgg**2;tpc=168.0+325.0*sgg-12.5*sgg**2
        ppr=ppsia/ppc;tpr=trank/tpc;z=1.0
        for _ in range(1002):
            reduced_density=0.27*ppr/(z*tpr)
            c1=0.3265-1.07/tpr-0.5339/tpr**3+0.01569/tpr**4-0.05165/tpr**5
            c2=0.5475-0.7361/tpr+0.1844/tpr**2;c3=-0.7361/tpr+0.1844/tpr**2
            previous=z
            z=1.0+c1*reduced_density+c2*reduced_density**2-0.1056*c3*reduced_density**5+0.6134*(1.0+0.721*reduced_density**2)*(reduced_density**2/tpr**3)*math.exp(-0.721*reduced_density**2)
            if abs(z-previous)<0.0001:break
        bg=0.02827*z*trank/ppsia;rhog0=sgg*1.22;rhoo0=sgo*997.0
        _rhoo=(rhoo0+rhog0*rs/5.6738)/bo;rhog=rhog0/bg
        denominator=rhog0*gor+rhoo0+997.0*wor
        value=(rhog*bg*(gorss-rs)/5.6738)/denominator
    except (ValueError, OverflowError, ZeroDivisionError) as error:
        raise ValidationError("Black Oil vaporized fraction is unrepresentable") from error
    if not math.isfinite(value):raise ValidationError("Black Oil vaporized fraction is unrepresentable")
    return value


def black_oil_fugacity_coefficients(data: BlackOilData, steam_data: SteamTablesData, compound_ids, composition, temperature_k: float, pressure_pa: float, phase: str) -> tuple[float, ...]:
    ids=tuple(compound_ids);_composition(data,ids,composition);_state(temperature_k,pressure_pa)
    if phase=="vapor":return tuple(1.0 for _ in ids)
    if phase!="liquid":raise ValidationError("Black Oil phase must be liquid or vapor")
    return tuple(black_oil_vapor_pressure_pa(data,steam_data,name,temperature_k)/pressure_pa for name in ids)


def black_oil_tp_flash(data: BlackOilData, compound_ids, composition, temperature_k: float, pressure_pa: float) -> BlackOilTPFlashResult:
    ids=tuple(compound_ids);z=_composition(data,ids,composition);_state(temperature_k,pressure_pa)
    fractions=tuple(black_oil_component_vaporized_fraction(data.compound(name),temperature_k,pressure_pa) for name in ids)
    vapor=sum(fraction*value for fraction,value in zip(fractions,z));liquid=1.0-vapor
    vapor_composition=tuple(fraction*value/vapor for fraction,value in zip(fractions,z)) if vapor>0.0 else tuple(0.0 for _ in z)
    liquid_composition=tuple(value*(1.0-fraction)/liquid for fraction,value in zip(fractions,z)) if liquid>0.0 else tuple(0.0 for _ in z)
    if vapor<=0.0:liquid=1.0;vapor=0.0;liquid_composition=z
    if vapor>=1.0:liquid=0.0;vapor=1.0;vapor_composition=z
    return BlackOilTPFlashResult(liquid,vapor,liquid_composition,vapor_composition,tuple(0.0 for _ in z),1)


__all__=("BlackOilCompound","BlackOilData","BlackOilTPFlashResult","load_black_oil_data","black_oil_vapor_pressure_pa","black_oil_component_vaporized_fraction","black_oil_fugacity_coefficients","black_oil_tp_flash")
