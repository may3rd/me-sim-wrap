"""Non-associating PC-SAFT methane/ethane equation slice."""
from dataclasses import dataclass
import json,math
from pathlib import Path
from ..errors import ValidationError

_A0=(.9105631445,.6361281449,2.6861347891,-26.547362491,97.759208784,-159.59154087,91.297774084)
_A1=(-.3084016918,.1860531159,-2.5030047259,21.419793629,-65.25588533,83.318680481,-33.74692293)
_A2=(-.0906148351,.4527842806,.5962700728,-1.7241829131,-4.1302112531,13.77663187,-8.6728470368)
_B0=(.7240946941,2.2382791861,-4.0025849485,-21.003576815,26.855641363,206.55133841,-355.60235612)
_B1=(-.5755498075,.6995095521,3.892567339,-17.215471648,192.67226447,-161.82646165,-165.20769346)
_B2=(.0976883116,-.2557574982,-9.155856153,20.642075974,-38.804430052,93.626774077,-29.666905585)
_KB=1.3806504e-23
@dataclass(frozen=True,slots=True)
class PCSAFTParameter:
 source_name:str;cas_number:str;molecular_weight:float;segment_count:float;segment_diameter_angstrom:float;dispersion_energy_k:float;association_volume:float;association_energy_k:float
@dataclass(frozen=True,slots=True)
class PCSAFTInteraction:
 cas_number_1:str;cas_number_2:str;kij:float
@dataclass(frozen=True,slots=True)
class PCSAFTData:
 source_revision:str;case_sha256:str;runtime_assembly_sha256:str;property_package_source_sha256:str;eos_source_sha256:str;pure_data_source_sha256:str;interaction_data_source_sha256:str;scoped_compounds:tuple[tuple[str,str],...];pure_parameters:tuple[PCSAFTParameter,...];interaction_parameters:tuple[PCSAFTInteraction,...]
 @property
 def compound_ids(self):return tuple(x[0] for x in self.scoped_compounds)
 def parameter(self,cas):
  try:return next(x for x in self.pure_parameters if x.cas_number==cas)
  except StopIteration as error:raise ValidationError(f"missing PC-SAFT pure parameter: {cas}") from error
 def kij(self,first,second):
  for x in self.interaction_parameters:
   if {x.cas_number_1,x.cas_number_2}=={first,second}:return x.kij
  return 0.0
@dataclass(frozen=True,slots=True)
class PCSAFTTPFlashResult:
 liquid_fraction:float;vapor_fraction:float;liquid_composition:tuple[float,...];vapor_composition:tuple[float,...];equilibrium_ratios:tuple[float,...];iterations:int;residual_norm:float
def load_pcsaft_data(path:str|Path)->PCSAFTData:
 try:
  d=json.loads(Path(path).read_text(encoding="utf-8-sig"));s=d["source"]
  if d["schema_version"]!="dwsim-pcsaft-data-1" or d["model"]!="PC-SAFT (with Association Support) (.NET Code)":raise ValidationError("unsupported PC-SAFT data schema or model")
  pure=tuple(PCSAFTParameter(x["source_name"],x["cas_number"],x["molecular_weight"],x["segment_count"],x["segment_diameter_angstrom"],x["dispersion_energy_k"],x["association_volume"],x["association_energy_k"]) for x in d["pure_parameters"])
  interactions=tuple(PCSAFTInteraction(x["cas_number_1"],x["cas_number_2"],x["kij"]) for x in d["interaction_parameters"])
  data=PCSAFTData(s["revision"],s["case_sha256"],s["runtime_assembly_sha256"],s["property_package_source_sha256"],s["eos_source_sha256"],s["pure_data_source_sha256"],s["interaction_data_source_sha256"],tuple((x["compound_id"],x["cas_number"]) for x in d["scoped_compounds"]),pure,interactions)
 except ValidationError:raise
 except (KeyError,TypeError,ValueError) as error:raise ValidationError(f"invalid PC-SAFT data: {error}") from error
 hashes=(data.case_sha256,data.runtime_assembly_sha256,data.property_package_source_sha256,data.eos_source_sha256,data.pure_data_source_sha256,data.interaction_data_source_sha256)
 if not data.source_revision or data.scoped_compounds!=(("Methane","74-82-8"),("Ethane","74-84-0")) or len(data.pure_parameters)!=94 or len(data.interaction_parameters)!=33 or len({x.cas_number for x in data.pure_parameters})!=93 or any(len(x)!=64 for x in hashes) or any(not x.source_name or not x.cas_number or x.molecular_weight<=0 or x.segment_count<=0 or x.segment_diameter_angstrom<=0 or x.dispersion_energy_k<=0 for x in data.pure_parameters) or any(not math.isfinite(x.kij) for x in data.interaction_parameters):raise ValidationError("PC-SAFT source identity or parameter domain is invalid")
 for _,cas in data.scoped_compounds:
  p=data.parameter(cas)
  if p.association_volume!=0 or p.association_energy_k!=0:raise ValidationError("scoped PC-SAFT compounds must be non-associating")
 return data
def _composition(data,compound_ids,composition):
 try:ids=tuple(compound_ids);z=tuple(composition)
 except TypeError as error:raise ValidationError("PC-SAFT inputs must be sequences") from error
 if not isinstance(data,PCSAFTData) or ids!=data.compound_ids or len(z)!=2 or any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0 for x in z) or not math.isclose(sum(z),1,rel_tol=0,abs_tol=1e-12):raise ValidationError("PC-SAFT requires its exact positive normalized compound domain")
 return tuple(float(x) for x in z)
def _model(data):
 p=tuple(data.parameter(cas) for _,cas in data.scoped_compounds);m=tuple(x.segment_count for x in p);sigma=tuple(x.segment_diameter_angstrom for x in p);epsilon=tuple(x.dispersion_energy_k for x in p);k=((0.0,data.kij(p[0].cas_number,p[1].cas_number)),(data.kij(p[0].cas_number,p[1].cas_number),0.0));return m,sigma,epsilon,k
def _state(temperature_k,pressure_pa):
 if isinstance(temperature_k,bool) or not isinstance(temperature_k,(int,float)) or not math.isfinite(temperature_k) or temperature_k<=0 or isinstance(pressure_pa,bool) or not isinstance(pressure_pa,(int,float)) or not math.isfinite(pressure_pa) or pressure_pa<=0:raise ValidationError("PC-SAFT state must be positive and finite")
def _properties(data,temperature,x,number_density):
 m,sigma,epsilon,kij=_model(data);diameter=tuple(sigma[i]*(1-.12*math.exp(-3*epsilon[i]/temperature)) for i in range(2));mean_m=sum(x[i]*m[i] for i in range(2));zeta=tuple(math.pi/6*number_density*sum(x[i]*m[i]*diameter[i]**j for i in range(2)) for j in range(4));eta=zeta[3]
 if not 0<eta<.99:raise ValidationError("PC-SAFT reduced density is outside the equation domain")
 g=[]
 for i in range(2):
  q=diameter[i]/2;g.append(1/(1-eta)+q*3*zeta[2]/(1-eta)**2+q*q*2*zeta[2]**2/(1-eta)**3)
 t1=3*zeta[1]*zeta[2]/(1-eta);t2=zeta[2]**3/(eta*(1-eta)**2);t3=(zeta[2]**3/eta**2-zeta[0])*math.log(1-eta);a_hs=(t1+t2+t3)/zeta[0];a_hc=mean_m*a_hs-sum(x[i]*(m[i]-1)*math.log(g[i]) for i in range(2))
 z_hs=eta/(1-eta)+3*zeta[1]*zeta[2]/(zeta[0]*(1-eta)**2)+(3*zeta[2]**3-eta*zeta[2]**3)/(zeta[0]*(1-eta)**3)
 dg=[]
 for i in range(2):
  q=diameter[i]/2;dg.append(eta/(1-eta)**2+q*(3*zeta[2]/(1-eta)**2+6*zeta[2]*eta/(1-eta)**3)+q*q*(4*zeta[2]**2/(1-eta)**3+6*zeta[2]**2*eta/(1-eta)**4))
 z_hc=mean_m*z_hs-sum(x[i]*(m[i]-1)/g[i]*dg[i] for i in range(2))
 a=tuple(_A0[j]+(mean_m-1)/mean_m*_A1[j]+(mean_m-1)/mean_m*(mean_m-2)/mean_m*_A2[j] for j in range(7));b=tuple(_B0[j]+(mean_m-1)/mean_m*_B1[j]+(mean_m-1)/mean_m*(mean_m-2)/mean_m*_B2[j] for j in range(7));i1=sum(a[j]*eta**j for j in range(7));i2=sum(b[j]*eta**j for j in range(7));term1=mean_m*(8*eta-2*eta**2)/(1-eta)**4;term2=(1-mean_m)*(20*eta-27*eta**2+12*eta**3-2*eta**4)/((1-eta)*(2-eta))**2;c1=1/(1+term1+term2)
 prom1=prom2=0.0
 for i in range(2):
  for j in range(2):
   sigma_ij=.5*(sigma[i]+sigma[j]);epsilon_ij=math.sqrt(epsilon[i]*epsilon[j])*(1-kij[i][j]);weight=x[i]*x[j]*m[i]*m[j]*sigma_ij**3;prom1+=weight*epsilon_ij/temperature;prom2+=weight*(epsilon_ij/temperature)**2
 a_disp=-2*math.pi*number_density*i1*prom1-math.pi*number_density*mean_m*c1*i2*prom2
 derivative_i1=sum(a[j]*(j+1)*eta**j for j in range(7));derivative_i2=sum(b[j]*(j+1)*eta**j for j in range(7));term1=mean_m*(-4*eta**2+20*eta+8)/(1-eta)**5;term2=(1-mean_m)*(2*eta**3+12*eta**2-48*eta+40)/((1-eta)*(2-eta))**3;c2=-c1**2*(term1+term2);z_disp=-2*math.pi*number_density*derivative_i1*prom1-math.pi*number_density*mean_m*(c1*derivative_i2+c2*eta*i2)*prom2
 return a_hc+a_disp,z_hc+z_disp,diameter
def _roots(data,temperature,pressure,x):
 from scipy.optimize import brentq
 m,_,_,_=_model(data);diameter=_properties(data,temperature,x,1e-12)[2];segment_volume=sum(x[i]*m[i]*diameter[i]**3 for i in range(2))
 def residual(eta):
  density=6*eta/(math.pi*segment_volume);_,z_res,_=_properties(data,temperature,x,density);return (1+z_res)*_KB*temperature*density*1e30-pressure
 upper=.72;start=-10;stop=math.log10(upper);grid=tuple(10**(start+(stop-start)*i/320) for i in range(321));roots=[];previous=grid[0];f_previous=residual(previous)
 for value in grid[1:]:
  f_value=residual(value)
  if f_previous*f_value<0:
   root=brentq(residual,previous,value,xtol=1e-14)
   if not roots or abs(root-roots[-1])>1e-9:roots.append(root)
  previous=value;f_previous=f_value
 if not roots:raise ValidationError("PC-SAFT density solve did not find a physical root")
 return tuple(roots),segment_volume
def pcsaft_fugacity_coefficients(data,compound_ids,composition,temperature_k,pressure_pa,phase):
 x=_composition(data,compound_ids,composition);_state(temperature_k,pressure_pa)
 if phase not in {"liquid","vapor"}:raise ValidationError("PC-SAFT phase must be liquid or vapor")
 roots,segment_volume=_roots(data,temperature_k,pressure_pa,x);eta=roots[-1] if phase=="liquid" else roots[0];density=6*eta/(math.pi*segment_volume);_,z_res,_=_properties(data,temperature_k,x,density);z_factor=1+z_res;step=1e-6;chemical=[]
 for index in range(2):
  def total(delta):
   amounts=list(x);amounts[index]+=delta;total_amount=sum(amounts);perturbed=tuple(value/total_amount for value in amounts);a_res,_,_=_properties(data,temperature_k,perturbed,density*total_amount);return total_amount*a_res
  chemical.append((total(step)-total(-step))/(2*step))
 return tuple(math.exp(value-math.log(z_factor)) for value in chemical)
def pcsaft_tp_flash(data,compound_ids,composition,temperature_k,pressure_pa):
 from scipy.optimize import least_squares
 z=_composition(data,compound_ids,composition);_state(temperature_k,pressure_pa);calls=[0]
 def residual(values):
  calls[0]+=1;beta,x0,y0=values;x=(x0,1-x0);y=(y0,1-y0);phil=pcsaft_fugacity_coefficients(data,compound_ids,x,temperature_k,pressure_pa,"liquid");phiv=pcsaft_fugacity_coefficients(data,compound_ids,y,temperature_k,pressure_pa,"vapor")
  return (z[0]-(1-beta)*x0-beta*y0,math.log(x0*phil[0]/(y0*phiv[0])),math.log((1-x0)*phil[1]/((1-y0)*phiv[1])))
 result=least_squares(residual,(.75,max(.05,z[0]*.5),min(.95,z[0]*1.16)),bounds=((1e-9,1e-9,1e-9),(1-1e-9,1-1e-9,1-1e-9)),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=120)
 norm=max(abs(float(x)) for x in result.fun)
 if not result.success or norm>1e-7:raise ValidationError("PC-SAFT TP flash did not converge")
 beta,x0,y0=(float(x) for x in result.x);x=(x0,1-x0);y=(y0,1-y0);return PCSAFTTPFlashResult(1-beta,beta,x,y,tuple(y[i]/x[i] for i in range(2)),calls[0],norm)
__all__=("PCSAFTParameter","PCSAFTInteraction","PCSAFTData","PCSAFTTPFlashResult","load_pcsaft_data","pcsaft_fugacity_coefficients","pcsaft_tp_flash")
