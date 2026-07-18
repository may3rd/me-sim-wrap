"""Scoped Chao-Seader liquid and vapor fugacity equations."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import MissingCompoundData, ValidationError
from .peng_robinson import _cubic_real_roots

R = 8.314

@dataclass(frozen=True, slots=True)
class ChaoSeaderCompound:
    compound_id: str
    molecular_weight: float
    critical_temperature_k: float
    critical_pressure_pa: float
    acentric_factor: float
    liquid_molar_volume: float
    chao_seader_acentricity: float
    solubility_parameter: float

@dataclass(frozen=True, slots=True)
class ChaoSeaderData:
    model: str
    source_revision: str
    case_sha256: str
    runtime_assembly_sha256: str
    property_package_source_sha256: str
    model_source_sha256: str
    compounds: tuple[ChaoSeaderCompound, ...]
    def compound(self, compound_id: str) -> ChaoSeaderCompound:
        for record in self.compounds:
            if record.compound_id == compound_id:
                return record
        raise MissingCompoundData(f"missing Chao-Seader compound: {compound_id}")

@dataclass(frozen=True, slots=True)
class ChaoSeaderTPFlashResult:
    liquid_fraction: float
    vapor_fraction: float
    liquid_composition: tuple[float, ...]
    vapor_composition: tuple[float, ...]
    equilibrium_ratios: tuple[float, ...]
    iterations: int

def load_chao_seader_data(path: str | Path) -> ChaoSeaderData:
    try:
        document=json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if document["schema_version"]!="dwsim-semi-empirical-data-1" or document["model"] not in {"Chao-Seader","Grayson-Streed"}:
            raise ValidationError("unsupported semi-empirical data schema or model")
        source=document["source"]
        compounds=tuple(ChaoSeaderCompound(**record) for record in document["compounds"])
        data=ChaoSeaderData(document["model"],source["revision"],source["case_sha256"],source["runtime_assembly_sha256"],source["property_package_source_sha256"],source["model_source_sha256"],compounds)
    except ValidationError:
        raise
    except (KeyError,TypeError,ValueError) as error:
        raise ValidationError(f"invalid Chao-Seader data: {error}") from error
    if not data.source_revision or len(data.compounds)!=2 or any(len(value)!=64 for value in (data.case_sha256,data.runtime_assembly_sha256,data.property_package_source_sha256,data.model_source_sha256)):
        raise ValidationError("Chao-Seader source identity or record count is invalid")
    if len({record.compound_id for record in data.compounds})!=len(data.compounds):
        raise ValidationError("Chao-Seader compound IDs must be unique")
    for record in data.compounds:
        positive=(record.molecular_weight,record.critical_temperature_k,record.critical_pressure_pa,record.liquid_molar_volume,record.solubility_parameter)
        finite=positive+(record.acentric_factor,record.chao_seader_acentricity)
        if not record.compound_id or any(not math.isfinite(value) for value in finite) or any(value<=0 for value in positive):
            raise ValidationError("invalid Chao-Seader compound record")
    return data

def _state(data: ChaoSeaderData, compound_ids, composition, temperature_k: float, pressure_pa: float):
    if not isinstance(data,ChaoSeaderData):
        raise ValidationError("Chao-Seader data is required")
    try:
        ids=tuple(compound_ids); fractions=tuple(composition)
    except TypeError as error:
        raise ValidationError("Chao-Seader inputs must be sequences") from error
    if len(ids)!=2 or len(fractions)!=2 or len(set(ids))!=2 or any(isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 for value in fractions) or not math.isclose(math.fsum(fractions),1.0,abs_tol=1e-12) or isinstance(temperature_k,bool) or not isinstance(temperature_k,(int,float)) or not math.isfinite(temperature_k) or temperature_k<=0 or isinstance(pressure_pa,bool) or not isinstance(pressure_pa,(int,float)) or not math.isfinite(pressure_pa) or pressure_pa<=0:
        raise ValidationError("invalid Chao-Seader state")
    return tuple(data.compound(compound_id) for compound_id in ids),fractions

def _positive_exp(log_value: float, label: str) -> float:
    try:
        value=math.exp(log_value)
    except OverflowError as error:
        raise ValidationError(f"{label} is outside the representable range") from error
    if not math.isfinite(value) or value<=0:
        raise ValidationError(f"{label} is outside the representable range")
    return value

def _semi_empirical_liquid_fugacity_coefficients(data: ChaoSeaderData,compound_ids,composition,temperature_k: float,pressure_pa: float)->tuple[float,...]:
    compounds,fractions=_state(data,compound_ids,composition,temperature_k,pressure_pa)
    try:
        volume_sum=math.fsum(x*c.liquid_molar_volume for x,c in zip(fractions,compounds))
        mixture_solubility=math.fsum(x*c.liquid_molar_volume*c.solubility_parameter for x,c in zip(fractions,compounds))/volume_sum
        output=[]
        for compound in compounds:
            reduced_temperature=temperature_k/compound.critical_temperature_k
            reduced_pressure=pressure_pa/compound.critical_pressure_pa
            molecular_weight_class=int(round(compound.molecular_weight))
            if data.model=="Chao-Seader" and molecular_weight_class==2:
                a=(1.96718,1.02972,-0.054009,0.0005288,0.0,0.008585,0.0,0.0,0.0,0.0)
            elif data.model=="Chao-Seader" and molecular_weight_class==16:
                a=(2.4384,-2.2455,-0.34084,0.00212,-0.00223,0.10486,-0.03691,0.0,0.0,0.0)
            elif data.model=="Chao-Seader":
                a=(5.75748,-3.01761,-4.985,2.02299,0.0,0.08427,0.26667,-0.31138,-0.02655,0.02883)
            elif molecular_weight_class==2:
                a=(1.50709,2.74283,-0.0211,0.00011,0.0,0.008585,0.0,0.0,0.0,0.0)
            elif molecular_weight_class==16:
                a=(1.36822,-1.54831,0.0,0.02889,-0.01076,0.10486,-0.02529,0.0,0.0,0.0)
            else:
                a=(2.05135,-2.10889,0.0,-0.19396,0.02282,0.08852,0.0,-0.00872,-0.00353,0.00203)
            log_nu=a[0]+a[1]/reduced_temperature+a[2]*reduced_temperature+a[3]*reduced_temperature**2+a[4]*reduced_temperature**3+(a[5]+a[6]*reduced_temperature+a[7]*reduced_temperature**2)*reduced_pressure+(a[8]+a[9]*reduced_temperature)*reduced_pressure**2-math.log10(reduced_pressure)
            correction=-4.23893+8.65808*reduced_temperature-1.2206/reduced_temperature-3.15224*reduced_temperature**3-0.025*(reduced_pressure-0.6)
            log_activity=compound.liquid_molar_volume*(compound.solubility_parameter-mixture_solubility)**2/(8314470.0*temperature_k)
            output.append(_positive_exp(math.log(10.0)*(log_nu+compound.chao_seader_acentricity*correction)+log_activity,"Chao-Seader liquid fugacity coefficient"))
        return tuple(output)
    except (ValueError,ZeroDivisionError,OverflowError) as error:
        raise ValidationError("Chao-Seader liquid state is outside the representable range") from error

def _semi_empirical_vapor_fugacity_coefficients(data: ChaoSeaderData,compound_ids,composition,temperature_k: float,pressure_pa: float)->tuple[float,...]:
    compounds,fractions=_state(data,compound_ids,composition,temperature_k,pressure_pa)
    try:
        pure_a=tuple(0.42748*R**2*c.critical_temperature_k**2.5/(c.critical_pressure_pa*math.sqrt(temperature_k)) for c in compounds)
        pure_b=tuple(0.08664*R*c.critical_temperature_k/c.critical_pressure_pa for c in compounds)
        mixture_a=math.fsum(fractions[i]*fractions[j]*math.sqrt(pure_a[i]*pure_a[j]) for i in range(2) for j in range(2))
        mixture_b=math.fsum(x*b for x,b in zip(fractions,pure_b))
        reduced_a=mixture_a*pressure_pa/(R*temperature_k)**2
        reduced_b=mixture_b*pressure_pa/(R*temperature_k)
        roots=_cubic_real_roots(-1.0,reduced_a-reduced_b-reduced_b**2,-reduced_a*reduced_b)
        compressibility=max(root for root in roots if root>reduced_b)
        log_ratio=math.log((compressibility+reduced_b)/compressibility)
        return tuple(_positive_exp(pure_b[i]*(compressibility-1.0)/mixture_b-math.log(compressibility-reduced_b)+(reduced_a/reduced_b)*(pure_b[i]/mixture_b-2.0*math.sqrt(pure_a[i]/mixture_a))*log_ratio,"Chao-Seader vapor fugacity coefficient") for i in range(2))
    except (ValueError,ZeroDivisionError,OverflowError) as error:
        raise ValidationError("Chao-Seader vapor state is outside the representable range") from error

def _rachford_rice_vapor_fraction(composition: tuple[float,...],ratios: tuple[float,...])->float:
    def residual(value: float)->float:
        return math.fsum(composition[i]*(ratios[i]-1.0)/(1.0+value*(ratios[i]-1.0)) for i in range(2))
    lower_residual=residual(0.0);upper_residual=residual(1.0)
    if lower_residual<=0.0 or upper_residual>=0.0:
        raise ValidationError("Chao-Seader scoped flash requires a two-phase state")
    lower=0.0;upper=1.0
    for _ in range(80):
        midpoint=(lower+upper)/2.0
        if residual(midpoint)>0.0: lower=midpoint
        else: upper=midpoint
    return (lower+upper)/2.0

def _semi_empirical_tp_flash(data: ChaoSeaderData,compound_ids,composition,temperature_k:float,pressure_pa:float,*,max_iterations:int=100,tolerance:float=1e-11)->ChaoSeaderTPFlashResult:
    compounds,feed=_state(data,compound_ids,composition,temperature_k,pressure_pa)
    if isinstance(max_iterations,bool) or not isinstance(max_iterations,int) or max_iterations<=0 or isinstance(tolerance,bool) or not isinstance(tolerance,(int,float)) or not math.isfinite(tolerance) or tolerance<=0:
        raise ValidationError("invalid Chao-Seader flash controls")
    try:
        ratios=tuple(c.critical_pressure_pa/pressure_pa*math.exp(5.373*(1.0+c.acentric_factor)*(1.0-c.critical_temperature_k/temperature_k)) for c in compounds)
        for iteration in range(1,max_iterations+1):
            vapor_fraction=_rachford_rice_vapor_fraction(feed,ratios)
            liquid=tuple(feed[i]/(1.0+vapor_fraction*(ratios[i]-1.0)) for i in range(2))
            vapor=tuple(ratios[i]*liquid[i] for i in range(2))
            liquid_total=math.fsum(liquid);vapor_total=math.fsum(vapor)
            liquid=tuple(value/liquid_total for value in liquid);vapor=tuple(value/vapor_total for value in vapor)
            liquid_phi=_semi_empirical_liquid_fugacity_coefficients(data,compound_ids,liquid,temperature_k,pressure_pa)
            vapor_phi=_semi_empirical_vapor_fugacity_coefficients(data,compound_ids,vapor,temperature_k,pressure_pa)
            updated=tuple(liquid_phi[i]/vapor_phi[i] for i in range(2))
            error=max(abs(math.log(updated[i]/ratios[i])) for i in range(2))
            ratios=updated
            if error<=tolerance:
                vapor_fraction=_rachford_rice_vapor_fraction(feed,ratios)
                liquid=tuple(feed[i]/(1.0+vapor_fraction*(ratios[i]-1.0)) for i in range(2))
                vapor=tuple(ratios[i]*liquid[i] for i in range(2))
                return ChaoSeaderTPFlashResult(1.0-vapor_fraction,vapor_fraction,liquid,vapor,ratios,iteration)
    except (ValueError,ZeroDivisionError,OverflowError) as error:
        raise ValidationError("Chao-Seader flash is outside the representable range") from error
    raise ValidationError("Chao-Seader flash did not converge")

def _require_model(data: ChaoSeaderData, expected: str) -> None:
    if not isinstance(data,ChaoSeaderData) or data.model!=expected:
        raise ValidationError(f"{expected} data is required")

def chao_seader_liquid_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa):
    _require_model(data,"Chao-Seader");return _semi_empirical_liquid_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa)
def chao_seader_vapor_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa):
    _require_model(data,"Chao-Seader");return _semi_empirical_vapor_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa)
def chao_seader_tp_flash(data,compound_ids,composition,temperature_k,pressure_pa,**controls):
    _require_model(data,"Chao-Seader");return _semi_empirical_tp_flash(data,compound_ids,composition,temperature_k,pressure_pa,**controls)
def grayson_streed_liquid_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa):
    _require_model(data,"Grayson-Streed");return _semi_empirical_liquid_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa)
def grayson_streed_vapor_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa):
    _require_model(data,"Grayson-Streed");return _semi_empirical_vapor_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa)
def grayson_streed_tp_flash(data,compound_ids,composition,temperature_k,pressure_pa,**controls):
    _require_model(data,"Grayson-Streed");return _semi_empirical_tp_flash(data,compound_ids,composition,temperature_k,pressure_pa,**controls)

__all__=("ChaoSeaderData","ChaoSeaderTPFlashResult","load_chao_seader_data","chao_seader_liquid_fugacity_coefficients","chao_seader_vapor_fugacity_coefficients","chao_seader_tp_flash","grayson_streed_liquid_fugacity_coefficients","grayson_streed_vapor_fugacity_coefficients","grayson_streed_tp_flash")
