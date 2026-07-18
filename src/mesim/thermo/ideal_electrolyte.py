"""Scoped ideal aqueous-electrolyte fugacity and liquid TP classification."""
from dataclasses import dataclass
import json,math
from pathlib import Path
from ..errors import ValidationError
from .ideal import IdealCorrelations

@dataclass(frozen=True,slots=True)
class IdealElectrolyteCompound:
 compound_id:str;formula:str;molecular_weight:float;is_ion:bool;is_salt:bool;charge:float
@dataclass(frozen=True,slots=True)
class IdealElectrolyteData:
 source_revision:str;case_sha256:str;runtime_assembly_sha256:str;property_package_source_sha256:str;base_source_sha256:str;flash_source_sha256:str;database_source_sha256:str;compounds:tuple[IdealElectrolyteCompound,...];probe_temperature_k:float;probe_pressure_pa:float;probe_composition:tuple[float,...];probe_solvent_mass:float;probe_molalities:tuple[float,...]
 @property
 def compound_ids(self):return tuple(x.compound_id for x in self.compounds)
@dataclass(frozen=True,slots=True)
class IdealElectrolyteTPFlashResult:
 liquid_fraction:float;vapor_fraction:float;liquid_composition:tuple[float,...];vapor_composition:tuple[float,...];equilibrium_ratios:tuple[float,...];iterations:int
def load_ideal_electrolyte_data(path:str|Path)->IdealElectrolyteData:
 try:
  d=json.loads(Path(path).read_text(encoding="utf-8-sig"));s=d["source"];p=d["scoped_probe"]
  if d["schema_version"]!="dwsim-ideal-electrolyte-data-1" or d["model"]!="Ideal Solution (Aqueous Electrolytes)":raise ValidationError("unsupported ideal-electrolyte data schema or model")
  records=tuple(IdealElectrolyteCompound(x["compound_id"],x["formula"],x["molecular_weight"],x["is_ion"],x["is_salt"],x["charge"]) for x in d["compounds"])
  data=IdealElectrolyteData(s["revision"],s["case_sha256"],s["runtime_assembly_sha256"],s["property_package_source_sha256"],s["base_source_sha256"],s["flash_source_sha256"],s["database_source_sha256"],records,p["temperature_k"],p["pressure_pa"],tuple(p["composition"]),p["solvent_mass_kg_per_mol_mixture"],tuple(p["molalities"]))
 except ValidationError:raise
 except (KeyError,TypeError,ValueError) as error:raise ValidationError(f"invalid ideal-electrolyte data: {error}") from error
 hashes=(data.case_sha256,data.runtime_assembly_sha256,data.property_package_source_sha256,data.base_source_sha256,data.flash_source_sha256,data.database_source_sha256);numeric=(data.probe_temperature_k,data.probe_pressure_pa,*data.probe_composition,data.probe_solvent_mass,*data.probe_molalities,*(v for r in data.compounds for v in (r.molecular_weight,r.charge)))
 if not data.source_revision or data.compound_ids!=("Water","Sodium (ion)","Chloride (ion)") or len(data.compounds)!=3 or any(len(x)!=64 for x in hashes) or any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in numeric) or tuple((r.is_ion,r.is_salt) for r in data.compounds)!=((False,False),(True,False),(True,False)) or any(r.molecular_weight<=0 for r in data.compounds):raise ValidationError("ideal-electrolyte source identity or scoped data are invalid")
 return data
def _composition(data,compound_ids,composition):
 try:ids=tuple(compound_ids);z=tuple(composition)
 except TypeError as error:raise ValidationError("ideal-electrolyte inputs must be sequences") from error
 if not isinstance(data,IdealElectrolyteData) or ids!=data.compound_ids or len(z)!=3 or any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0 for x in z) or not math.isclose(sum(z),1,rel_tol=0,abs_tol=1e-12) or z[0]<=0 or not math.isclose(sum(x*r.charge for x,r in zip(z,data.compounds)),0,rel_tol=0,abs_tol=1e-12):raise ValidationError("ideal electrolyte requires normalized electroneutral Water/Na+/Cl- composition")
 return tuple(float(x) for x in z)
def ideal_electrolyte_molalities(data,compound_ids,composition):
 z=_composition(data,compound_ids,composition);solvent_mass=z[0]*data.compounds[0].molecular_weight/1000.0
 return tuple(x/solvent_mass for x in z)
def ideal_electrolyte_fugacity_coefficients(data,water,compound_ids,composition,temperature_k,pressure_pa,phase):
 z=_composition(data,compound_ids,composition)
 if not isinstance(water,IdealCorrelations) or water.compound_id!="Water" or phase not in {"liquid","vapor"} or isinstance(pressure_pa,bool) or not isinstance(pressure_pa,(int,float)) or not math.isfinite(pressure_pa) or pressure_pa<=0:raise ValidationError("invalid ideal-electrolyte state")
 if phase=="vapor":return (1.0,1.0e10,1.0e10)
 molality=ideal_electrolyte_molalities(data,compound_ids,z)
 return (water.vapor_pressure(temperature_k).value/pressure_pa,molality[1],molality[2])
def ideal_electrolyte_tp_flash(data,water,compound_ids,composition,temperature_k,pressure_pa):
 z=_composition(data,compound_ids,composition);ideal_electrolyte_fugacity_coefficients(data,water,compound_ids,z,temperature_k,pressure_pa,"liquid")
 if pressure_pa<water.vapor_pressure(temperature_k).value:raise ValidationError("ideal-electrolyte TP flash is scoped to the captured all-liquid domain")
 zero=(0.0,0.0,0.0);return IdealElectrolyteTPFlashResult(1.0,0.0,z,zero,zero,1)
__all__=("IdealElectrolyteCompound","IdealElectrolyteData","IdealElectrolyteTPFlashResult","load_ideal_electrolyte_data","ideal_electrolyte_molalities","ideal_electrolyte_fugacity_coefficients","ideal_electrolyte_tp_flash")
