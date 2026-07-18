"""Lee-Kesler-Plocker corresponding-states fugacity model."""

from dataclasses import dataclass
import json
import math
from pathlib import Path

from ..errors import MissingCompoundData, ValidationError

R=8.314
WH=0.3978
SIMPLE=(0.1181193,0.265728,0.15479,0.030323,0.0236744,0.0186984,0.0,0.042724,0.0000155488,0.0000623689,0.65392,0.060167)
HEAVY=(0.2026579,0.331511,0.027655,0.203488,0.0313385,0.0503618,0.016901,0.041577,0.000048736,0.00000740336,1.226,0.03754)

@dataclass(frozen=True,slots=True)
class LKPCompound:
    compound_id:str;molecular_weight:float;critical_temperature_k:float;critical_pressure_pa:float;acentric_factor:float;critical_volume:float
@dataclass(frozen=True,slots=True)
class LKPInteraction:
    first_compound_id:str;second_compound_id:str;kij:float
@dataclass(frozen=True,slots=True)
class LKPData:
    source_revision:str;case_sha256:str;runtime_assembly_sha256:str;interaction_resource_sha256:str;property_package_source_sha256:str;model_source_sha256:str
    compounds:tuple[LKPCompound,...];interaction_pairs:tuple[LKPInteraction,...]
    def compound(self,compound_id:str)->LKPCompound:
        for record in self.compounds:
            if record.compound_id==compound_id:return record
        raise MissingCompoundData(f"missing Lee-Kesler-Plocker compound: {compound_id}")
    def interaction(self,first:str,second:str)->float:
        if first==second:return 1.0
        for record in self.interaction_pairs:
            if (record.first_compound_id,record.second_compound_id) in {(first,second),(second,first)}:return record.kij if record.kij!=0 else 1.0
        return 1.0
@dataclass(frozen=True,slots=True)
class LKPTPFlashResult:
    liquid_fraction:float;vapor_fraction:float;liquid_composition:tuple[float,...];vapor_composition:tuple[float,...];equilibrium_ratios:tuple[float,...];iterations:int

def load_lkp_data(path:str|Path)->LKPData:
    try:
        d=json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if d["schema_version"]!="dwsim-lkp-data-1" or d["model"]!="Lee-Kesler-Plöcker":raise ValidationError("unsupported LKP data schema or model")
        s=d["source"];compounds=tuple(LKPCompound(**x) for x in d["compounds"]);pairs=tuple(LKPInteraction(**x) for x in d["interaction_pairs"])
        data=LKPData(s["revision"],s["case_sha256"],s["runtime_assembly_sha256"],s["interaction_resource_sha256"],s["property_package_source_sha256"],s["model_source_sha256"],compounds,pairs)
    except ValidationError:raise
    except (KeyError,TypeError,ValueError) as error:raise ValidationError(f"invalid LKP data: {error}") from error
    hashes=(data.case_sha256,data.runtime_assembly_sha256,data.interaction_resource_sha256,data.property_package_source_sha256,data.model_source_sha256)
    if not data.source_revision or len(data.compounds)!=2 or len(data.interaction_pairs)!=140 or any(len(x)!=64 for x in hashes):raise ValidationError("LKP source identity or record count is invalid")
    if len({x.compound_id for x in data.compounds})!=2 or len({(x.first_compound_id,x.second_compound_id) for x in data.interaction_pairs})!=140:raise ValidationError("LKP keys must be unique")
    for c in data.compounds:
        values=(c.molecular_weight,c.critical_temperature_k,c.critical_pressure_pa,c.critical_volume)
        if not c.compound_id or any(not math.isfinite(v) or v<=0 for v in values) or not math.isfinite(c.acentric_factor):raise ValidationError("invalid LKP compound")
    for p in data.interaction_pairs:
        if not p.first_compound_id or not p.second_compound_id or not math.isfinite(p.kij):raise ValidationError("invalid LKP interaction")
    return data

def _parameters(tr:float,constants):
    b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,beta,gamma=constants
    return b1-b2/tr-b3/tr**2-b4/tr**3,c1-c2/tr+c3/tr**3,d1+d2/tr,c4,beta,gamma
def _residual(vr:float,tr:float,pr:float,constants)->float:
    b,c,d,c4,beta,gamma=_parameters(tr,constants);v2=vr*vr
    return pr*vr/tr-(1+b/vr+c/v2+d/vr**5+c4/tr**3/v2*(beta+gamma/v2)*math.exp(-gamma/v2))
def _reduced_volume(phase:str,tr:float,pr:float,constants)->float:
    points=[10**(-5+i*(math.log10(2000)+5)/4000) for i in range(4001)];roots=[];left=points[0];fl=_residual(left,tr,pr,constants)
    for right in points[1:]:
        fr=_residual(right,tr,pr,constants)
        if fl*fr<0:
            a=left;b=right
            for _ in range(80):
                m=(a+b)/2;fm=_residual(m,tr,pr,constants)
                if fl*fm<=0:b=m
                else:a=m;fl=fm
            roots.append((a+b)/2)
        left=right;fl=fr
    if not roots:raise ValidationError("LKP reduced-volume root was not found")
    return min(roots) if phase=="liquid" else max(roots)
def _z_values(phase:str,tr:float,pr:float,w:float):
    zs=pr*_reduced_volume(phase,tr,pr,SIMPLE)/tr;zh=pr*_reduced_volume(phase,tr,pr,HEAVY)/tr
    return zs+w/WH*(zh-zs),zs,zh
def _e_term(vr:float,tr:float,constants)->float:
    *_,c4,_d1,_d2,beta,gamma=constants
    return c4/(2*tr**3*gamma)*(beta+1-(beta+1+gamma/vr**2)*math.exp(-gamma/vr**2))
def _h_lk(phase:str,tr:float,pr:float,w:float)->float:
    _,zs,zh=_z_values(phase,tr,pr,w);values=[]
    for z,constants in ((zs,SIMPLE),(zh,HEAVY)):
        _b1,b2,b3,b4,_c1,c2,c3,_c4,_d1,d2,_beta,_gamma=constants;vr=z*tr/pr;e=_e_term(vr,tr,constants)
        values.append(tr*(z-1-(b2+2*b3/tr+3*b4/tr**2)/(tr*vr)-(c2-3*c3/tr**2)/(2*tr*vr**2)+d2/(5*tr*vr**5)+3*e))
    return values[0]+w/WH*(values[1]-values[0])
def _ln_fugacity_mixture(phase:str,tr:float,pr:float,w:float):
    _,zs,zh=_z_values(phase,tr,pr,w);values=[]
    for z,constants in ((zs,SIMPLE),(zh,HEAVY)):
        b,c,d,*_= _parameters(tr,constants);vr=z*tr/pr
        values.append(z-1-math.log(z)+b/vr+c/(2*vr**2)+d/(5*vr**5)+_e_term(vr,tr,constants))
    return values[0]+w/WH*(values[1]-values[0]),values[0],values[1]
def _state(data,compound_ids,composition,temperature_k,pressure_pa):
    if not isinstance(data,LKPData):raise ValidationError("LKP data is required")
    try:ids=tuple(compound_ids);z=tuple(composition)
    except TypeError as error:raise ValidationError("LKP inputs must be sequences") from error
    if len(ids)!=2 or len(z)!=2 or len(set(ids))!=2 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in z) or not math.isclose(math.fsum(z),1,abs_tol=1e-12) or isinstance(temperature_k,bool) or not isinstance(temperature_k,(int,float)) or not math.isfinite(temperature_k) or temperature_k<=0 or isinstance(pressure_pa,bool) or not isinstance(pressure_pa,(int,float)) or not math.isfinite(pressure_pa) or pressure_pa<=0:raise ValidationError("invalid LKP state")
    return tuple(data.compound(x) for x in ids),z,ids
def _mixture(data,compounds,z,ids):
    vc=tuple(tuple((compounds[i].critical_volume**(1/3)+compounds[j].critical_volume**(1/3))**3/8/1000 for j in range(2)) for i in range(2));tc=tuple(tuple(math.sqrt(compounds[i].critical_temperature_k*compounds[j].critical_temperature_k)*data.interaction(ids[i],ids[j]) for j in range(2)) for i in range(2))
    vcm=math.fsum(z[i]*z[j]*vc[i][j] for i in range(2) for j in range(2));tcm=math.fsum(z[i]*z[j]*vc[i][j]**0.25*tc[i][j]/vcm**0.25 for i in range(2) for j in range(2));wm=math.fsum(z[i]*compounds[i].acentric_factor for i in range(2));pcm=(0.2905-0.085*wm)*R*tcm/vcm
    return tcm,pcm,vcm,wm,vc,tc

def lkp_fugacity_coefficients(data:LKPData,compound_ids,composition,temperature_k:float,pressure_pa:float,phase:str)->tuple[float,...]:
    if phase not in {"liquid","vapor"}:raise ValidationError("LKP phase must be liquid or vapor")
    compounds,z,ids=_state(data,compound_ids,composition,temperature_k,pressure_pa)
    try:
        tcm,pcm,vcm,wm,vc,tc=_mixture(data,compounds,z,ids);tr=temperature_k/tcm;pr=pressure_pa/pcm;zcm=pcm*vcm/(R*tcm);zm=_z_values(phase,tr,pr,wm)[0];lnm,lns,lnh=_ln_fugacity_mixture(phase,tr,pr,wm);dlnw=(lnh-lns)/WH;h=_h_lk(phase,tr,pr,wm);out=[]
        for i in range(2):
            suma=sumb=sumc=0.0
            for j in range(2):
                if j==i:continue
                sum1=math.fsum(z[l]*(vc[l][j]**0.25*tc[l][j]-vc[l][i]**0.25*tc[l][i]) for l in range(2));sum2=math.fsum(z[l]*(vc[l][j]-vc[l][i]) for l in range(2));dz=-0.085*(compounds[j].acentric_factor-compounds[i].acentric_factor);dv=2*sum2;dt=(2*sum1-0.25*vcm**(-0.75)*dv*tcm)/vcm**0.25;dp=pcm*(dz/zcm+dt/tcm-dv/vcm);suma+=z[j]*dt;sumb+=z[j]*dp;sumc+=z[j]*(compounds[j].acentric_factor-compounds[i].acentric_factor)
            logphi=lnm-h*suma/temperature_k+(zm-1)*sumb/pcm-dlnw*sumc;value=math.exp(logphi)
            if not math.isfinite(value) or value<=0:raise ValidationError("LKP fugacity coefficient is unrepresentable")
            out.append(value)
        return tuple(out)
    except (ValueError,OverflowError,ZeroDivisionError) as error:raise ValidationError("LKP state is outside the representable range") from error

def _vapor_fraction(feed,ratios):
    def f(v):return math.fsum(feed[i]*(ratios[i]-1)/(1+v*(ratios[i]-1)) for i in range(2))
    if f(0)<=0 or f(1)>=0:raise ValidationError("scoped LKP flash requires a two-phase state")
    lo=0.0;hi=1.0
    for _ in range(80):
        mid=(lo+hi)/2
        if f(mid)>0:lo=mid
        else:hi=mid
    return (lo+hi)/2
def lkp_tp_flash(data:LKPData,compound_ids,composition,temperature_k:float,pressure_pa:float,*,max_iterations:int=100,tolerance:float=1e-10)->LKPTPFlashResult:
    compounds,feed,ids=_state(data,compound_ids,composition,temperature_k,pressure_pa)
    if isinstance(max_iterations,bool) or not isinstance(max_iterations,int) or max_iterations<=0 or isinstance(tolerance,bool) or not isinstance(tolerance,(int,float)) or not math.isfinite(tolerance) or tolerance<=0:raise ValidationError("invalid LKP flash controls")
    ratios=tuple(c.critical_pressure_pa/pressure_pa*math.exp(5.373*(1+c.acentric_factor)*(1-c.critical_temperature_k/temperature_k)) for c in compounds)
    for iteration in range(1,max_iterations+1):
        vapor_fraction=_vapor_fraction(feed,ratios);liquid=tuple(feed[i]/(1+vapor_fraction*(ratios[i]-1)) for i in range(2));vapor=tuple(ratios[i]*liquid[i] for i in range(2));liquid=tuple(v/math.fsum(liquid) for v in liquid);vapor=tuple(v/math.fsum(vapor) for v in vapor)
        liquid_phi=lkp_fugacity_coefficients(data,ids,liquid,temperature_k,pressure_pa,"liquid");vapor_phi=lkp_fugacity_coefficients(data,ids,vapor,temperature_k,pressure_pa,"vapor");updated=tuple(liquid_phi[i]/vapor_phi[i] for i in range(2));error=max(abs(math.log(updated[i]/ratios[i])) for i in range(2));ratios=updated
        if error<=tolerance:
            vapor_fraction=_vapor_fraction(feed,ratios);liquid=tuple(feed[i]/(1+vapor_fraction*(ratios[i]-1)) for i in range(2));vapor=tuple(ratios[i]*liquid[i] for i in range(2));return LKPTPFlashResult(1-vapor_fraction,vapor_fraction,liquid,vapor,ratios,iteration)
    raise ValidationError("LKP flash did not converge")

__all__=("LKPData","LKPTPFlashResult","load_lkp_data","lkp_fugacity_coefficients","lkp_tp_flash")
