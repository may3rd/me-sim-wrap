"""Scoped IAPWS-IF97 region-4 equilibrium boundary."""
from dataclasses import dataclass
import json,math
from pathlib import Path
from ..errors import ValidationError
from .ideal import IdealCorrelations
@dataclass(frozen=True,slots=True)
class SteamTablesData:
 source_revision:str;case_sha256:str;runtime_assembly_sha256:str;property_package_source_sha256:str;iapws_source_sha256:str;temperature_range_k:tuple[float,float];region4_coefficients:tuple[float,...]
@dataclass(frozen=True,slots=True)
class SteamTablesTPFlashResult:
 liquid_fraction:float;vapor_fraction:float;liquid_composition:tuple[float,...];vapor_composition:tuple[float,...];equilibrium_ratios:tuple[float,...];iterations:int
def load_steam_tables_data(path:str|Path)->SteamTablesData:
 try:
  d=json.loads(Path(path).read_text(encoding="utf-8-sig"));s=d["source"]
  if d["schema_version"]!="dwsim-steam-tables-data-1" or d["model"]!="Steam Tables":raise ValidationError("unsupported Steam Tables data schema or model")
  data=SteamTablesData(s["revision"],s["case_sha256"],s["runtime_assembly_sha256"],s["property_package_source_sha256"],s["iapws_source_sha256"],tuple(d["temperature_range_k"]),tuple(d["region4_coefficients"]))
 except ValidationError:raise
 except (KeyError,TypeError,ValueError) as error:raise ValidationError(f"invalid Steam Tables data: {error}") from error
 if not data.source_revision or len(data.temperature_range_k)!=2 or data.temperature_range_k!=(273.15,647.096) or len(data.region4_coefficients)!=10 or any(len(x)!=64 for x in (data.case_sha256,data.runtime_assembly_sha256,data.property_package_source_sha256,data.iapws_source_sha256)) or any(not math.isfinite(x) for x in data.region4_coefficients):raise ValidationError("Steam Tables source identity or coefficients are invalid")
 return data
def steam_saturation_pressure_pa(data:SteamTablesData,temperature_k:float)->float:
 if not isinstance(data,SteamTablesData) or isinstance(temperature_k,bool) or not isinstance(temperature_k,(int,float)) or not math.isfinite(temperature_k) or not data.temperature_range_k[0]<=temperature_k<=data.temperature_range_k[1]:raise ValidationError("temperature is outside the Steam Tables saturation range")
 n=(0.0,)+data.region4_coefficients;delta=temperature_k+n[9]/(temperature_k-n[10]);a=delta**2+n[1]*delta+n[2];b=n[3]*delta**2+n[4]*delta+n[5];c=n[6]*delta**2+n[7]*delta+n[8]
 try:value=(2*c/(-b+math.sqrt(b*b-4*a*c)))**4*1_000_000.0
 except (ValueError,OverflowError,ZeroDivisionError) as error:raise ValidationError("Steam Tables saturation pressure is unrepresentable") from error
 if not math.isfinite(value) or value<=0:raise ValidationError("Steam Tables saturation pressure is unrepresentable")
 return value
def steam_tables_fugacity_coefficients(data:SteamTablesData,water:IdealCorrelations,composition,temperature_k:float,pressure_pa:float,phase:str)->tuple[float,...]:
 try:z=tuple(composition)
 except TypeError as error:raise ValidationError("Steam Tables composition must be a sequence") from error
 if z!=(1,) and z!=(1.0,):raise ValidationError("Steam Tables supports pure Water only")
 if not isinstance(water,IdealCorrelations) or water.compound_id!="Water" or phase not in {"liquid","vapor"} or isinstance(pressure_pa,bool) or not isinstance(pressure_pa,(int,float)) or not math.isfinite(pressure_pa) or pressure_pa<=0:raise ValidationError("invalid Steam Tables state")
 return (water.vapor_pressure(temperature_k).value/pressure_pa,) if phase=="liquid" else (1.0,)
def steam_tables_tp_flash(data:SteamTablesData,water:IdealCorrelations,composition,temperature_k:float,pressure_pa:float)->SteamTablesTPFlashResult:
 steam_tables_fugacity_coefficients(data,water,composition,temperature_k,pressure_pa,"liquid");steam_tables_fugacity_coefficients(data,water,composition,temperature_k,pressure_pa,"vapor")
 liquid=1.0 if pressure_pa>=steam_saturation_pressure_pa(data,temperature_k) else 0.0
 return SteamTablesTPFlashResult(liquid,1-liquid,(1.0,),(1.0,),(0.0,),0)
__all__=("SteamTablesData","SteamTablesTPFlashResult","load_steam_tables_data","steam_saturation_pressure_pa","steam_tables_fugacity_coefficients","steam_tables_tp_flash")
