"""DWSIM Modified UNIFAC data and activity equation."""

from dataclasses import dataclass
import json
import math
from pathlib import Path
from ..errors import MissingCompoundData, ValidationError

@dataclass(frozen=True, slots=True)
class ModfacGroup:
    primary_id: int; secondary_id: int; primary_name: str; group_name: str; r: float; q: float
@dataclass(frozen=True, slots=True)
class ModfacSurfaceFraction:
    secondary_id: int; value: float
@dataclass(frozen=True, slots=True)
class ModfacCompoundBasis:
    compound_id: str; q: float; r: float; group_surface_fractions: tuple[ModfacSurfaceFraction, ...]
@dataclass(frozen=True, slots=True)
class ModfacInteractionPair:
    first_primary_id: int; second_primary_id: int
    a12: float; b12: float; c12: float; a21: float; b21: float; c21: float
@dataclass(frozen=True, slots=True)
class ModfacData:
    model: str; source_revision: str; groups_sha256: str; interactions_sha256: str
    compound_basis: tuple[ModfacCompoundBasis, ...]
    groups: tuple[ModfacGroup, ...]; interaction_pairs: tuple[ModfacInteractionPair, ...]
    def compound(self, compound_id: str) -> ModfacCompoundBasis:
        for record in self.compound_basis:
            if record.compound_id == compound_id: return record
        raise MissingCompoundData(f"missing {self.model} compound basis: {compound_id}")
    def group(self, secondary_id: int) -> ModfacGroup:
        for record in self.groups:
            if record.secondary_id == secondary_id: return record
        raise ValidationError(f"missing {self.model} subgroup: {secondary_id}")
    def coefficients(self, first: int, second: int) -> tuple[float, float, float]:
        if first == second: return (0.0, 0.0, 0.0)
        for record in self.interaction_pairs:
            if record.first_primary_id == first and record.second_primary_id == second:
                return (record.a12, record.b12, record.c12)
        for record in self.interaction_pairs:
            if record.first_primary_id == second and record.second_primary_id == first:
                return (record.a21, record.b21, record.c21)
        return (0.0, 0.0, 0.0)

def load_modfac_data(path: str | Path) -> ModfacData:
    try:
        doc=json.loads(Path(path).read_text(encoding="utf-8-sig"))
        expected_counts={"Modified UNIFAC (Dortmund)":(108,1167),"Modified UNIFAC (NIST)":(201,1969)}
        if doc["schema_version"]!="dwsim-modfac-data-1" or doc["model"] not in expected_counts:
            raise ValidationError("unsupported Modified UNIFAC data schema or model")
        s=doc["source"]
        compounds=tuple(ModfacCompoundBasis(x["compound_id"],x["q"],x["r"],tuple(ModfacSurfaceFraction(v["secondary_id"],v["value"]) for v in x["group_surface_fractions"])) for x in doc["compound_basis"])
        groups=tuple(ModfacGroup(x["primary_id"],x["secondary_id"],x["primary_name"],x["group_name"],x["r"],x["q"]) for x in doc["groups"])
        pairs=tuple(ModfacInteractionPair(x["first_primary_id"],x["second_primary_id"],x["a12"],x["b12"],x["c12"],x["a21"],x["b21"],x["c21"]) for x in doc["interaction_pairs"])
        data=ModfacData(doc["model"],s["revision"],s["groups_sha256"],s["interactions_sha256"],compounds,groups,pairs)
    except ValidationError: raise
    except (KeyError,TypeError,ValueError) as error: raise ValidationError(f"invalid Modified UNIFAC data: {error}") from error
    group_count,pair_count=expected_counts[data.model]
    if not data.source_revision or len(data.groups_sha256)!=64 or len(data.interactions_sha256)!=64 or len(data.compound_basis)!=2 or len(data.groups)!=group_count or len(data.interaction_pairs)!=pair_count:
        raise ValidationError(f"{data.model} source identity or record count is invalid")
    if len({x.compound_id for x in data.compound_basis})!=2 or len({x.secondary_id for x in data.groups})!=len(data.groups) or len({(x.first_primary_id,x.second_primary_id) for x in data.interaction_pairs})!=len(data.interaction_pairs):
        raise ValidationError("Modified UNIFAC keys must be unique")
    for g in data.groups:
        if g.primary_id<=0 or g.secondary_id<=0 or not g.primary_name or not g.group_name or not math.isfinite(g.r) or g.r<=0 or not math.isfinite(g.q) or g.q<0: raise ValidationError("invalid Modified UNIFAC group")
    for p in data.interaction_pairs:
        if p.first_primary_id<=0 or p.second_primary_id<=0 or p.first_primary_id==p.second_primary_id or any(not math.isfinite(v) for v in (p.a12,p.b12,p.c12,p.a21,p.b21,p.c21)): raise ValidationError("invalid Modified UNIFAC interaction")
    for c in data.compound_basis:
        if not c.compound_id or c.q<=0 or c.r<=0 or not c.group_surface_fractions or not math.isclose(math.fsum(v.value for v in c.group_surface_fractions),1.0,abs_tol=1e-12): raise ValidationError("invalid Modified UNIFAC compound basis")
        for v in c.group_surface_fractions:
            if v.value<=0 or not math.isfinite(v.value): raise ValidationError("invalid Modified UNIFAC surface fraction")
            data.group(v.secondary_id)
    return data

def modfac_activity_coefficients(data: ModfacData, compound_ids, composition, temperature_k: float) -> tuple[float,...]:
    if not isinstance(data,ModfacData): raise ValidationError("Modified UNIFAC data is required")
    try: ids, fractions=tuple(compound_ids),tuple(composition)
    except TypeError as error: raise ValidationError("Modified UNIFAC inputs must be sequences") from error
    if len(ids)!=2 or len(fractions)!=2 or len(set(ids))!=2 or isinstance(temperature_k,bool) or not isinstance(temperature_k,(int,float)) or not math.isfinite(temperature_k) or temperature_k<=0 or any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0 for v in fractions) or not math.isclose(math.fsum(fractions),1.0,abs_tol=1e-12): raise ValidationError("invalid Modified UNIFAC state")
    basis=tuple(data.compound(i) for i in ids); subgroups=tuple(dict.fromkeys(v.secondary_id for c in basis for v in c.group_surface_fractions)); eki=tuple({v.secondary_id:v.value for v in c.group_surface_fractions} for c in basis)
    try:
        tau={}
        for first in subgroups:
            for second in subgroups:
                a,b,c=data.coefficients(data.group(first).primary_id,data.group(second).primary_id)
                tau[(first,second)]=math.exp(-(a+b*temperature_k+c*temperature_k**2)/temperature_k)
        beta=tuple({m:math.fsum(e.get(k,0)*tau[(k,m)] for k in subgroups) for m in subgroups} for e in eki)
        sum_xq=math.fsum(fractions[i]*basis[i].q for i in range(2));theta={k:math.fsum(fractions[i]*basis[i].q*eki[i].get(k,0) for i in range(2))/sum_xq for k in subgroups};s={k:math.fsum(theta[m]*tau[(m,k)] for m in subgroups) for k in subgroups}
        sum_xr=math.fsum(fractions[i]*basis[i].r for i in range(2));sum_xr_power=math.fsum(fractions[i]*basis[i].r**0.75 for i in range(2));out=[]
        for i,c in enumerate(basis):
            j=c.r/sum_xr;jp=c.r**0.75/sum_xr_power;l=c.q/sum_xq
            lngc=1-jp+math.log(jp)-5*c.q*(1-j/l+math.log(j/l))
            total=math.fsum(theta[k]*beta[i][k]/s[k]-(eki[i][k]*math.log(beta[i][k]/s[k]) if k in eki[i] else 0) for k in subgroups)
            gamma=math.exp(lngc+c.q*(1-total))
            if not math.isfinite(gamma) or gamma<=0: raise ValidationError("Modified UNIFAC activity coefficient is unrepresentable")
            out.append(gamma)
    except (OverflowError,ValueError,ZeroDivisionError) as error: raise ValidationError("Modified UNIFAC state is outside the representable range") from error
    return tuple(out)

__all__=("ModfacData","load_modfac_data","modfac_activity_coefficients")
