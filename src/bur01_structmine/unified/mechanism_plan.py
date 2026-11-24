"""Freeze candidate selection, full-construct mappings and a bounded MD budget."""
from __future__ import annotations
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from .contracts import Contract
from pydantic import Field,model_validator

class MechanismBatch(Contract):
    batch_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
    family_id: str = "bur01"
    candidates: list[dict] = Field(min_length=1)
    reaction_id: str
    replicas: int = Field(default=3,ge=1,le=20)
    production_ns: float = Field(default=10.,gt=0,le=1000)
    segment_ns: float = Field(default=10.,gt=0,le=10)
    approved_budget_ns: float = Field(default=120.,gt=0)
    seeds: list[int] = Field(default_factory=list)
    enzyme_state: str = "free"
    @model_validator(mode="after")
    def budget(self):
        if self.enzyme_state!="free":raise ValueError("Acylated-state preparation is not yet supported")
        if len(self.candidates)*self.replicas*self.production_ns>self.approved_budget_ns+1e-9:
            raise ValueError("Requested trajectories exceed the frozen calculation budget")
        ids=[r["candidate_id"] for r in self.candidates]
        import re
        if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}",x) for x in ids):
            raise ValueError("Candidate IDs must be safe file identifiers")
        if len(set(ids))!=len(ids):raise ValueError("Candidate identities must be unique")
        count=len(ids)*self.replicas
        if not self.seeds:self.seeds=[20261124+i*7919 for i in range(count)]
        if len(self.seeds)!=count or len(set(self.seeds))!=count or min(self.seeds)<=0:
            raise ValueError("Every independent replica needs a distinct positive seed")
        return self

def pareto_front(frame,directions):
    if not directions or any(v not in {"max","min"} for v in directions.values()):
        raise ValueError("Specify max or min for each screening objective")
    cols=list(directions)
    values=frame[cols].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    present=np.isfinite(values).all(1)
    indices=np.flatnonzero(present)
    signed=values* np.array([1 if directions[c]=="max" else -1 for c in cols])
    selected=[]
    for i in indices:
        dominates=np.all(signed[indices]>=signed[i],axis=1)&np.any(signed[indices]>signed[i],axis=1)
        if not dominates.any():selected.append(i)
    result=frame.copy();result["pareto_eligible"]=present;result["pareto_front"]=False
    result.loc[result.index[selected],"pareto_front"]=True
    return result

def select_candidates(frame,parent,columns):
    required={"candidate_id","sequence_sha256","structure_review_pass",*columns}
    if required-set(frame.columns):raise ValueError("Missing selection fields: "+", ".join(sorted(required-set(frame.columns))))
    eligible=frame[frame.structure_review_pass.eq(True)].copy()
    chosen=[dict(parent,selection_reason="parent")]
    used={parent["sequence_sha256"]}
    for column in columns:
        candidates=eligible.copy();candidates[column]=pd.to_numeric(candidates[column],errors="coerce")
        candidates=candidates.dropna(subset=[column]).sort_values([column,"candidate_id"],ascending=[False,True])
        candidate=next((r for r in candidates.to_dict("records") if r["sequence_sha256"] not in used),None)
        if candidate is None:raise ValueError("No distinct, structurally reviewed candidate for "+column)
        candidate["selection_reason"]=column;chosen.append(candidate);used.add(candidate["sequence_sha256"])
    return chosen

def map_construct(core_parent,core_candidate,full_parent,*,offset=None):
    """Map substitutions only after locating one exact parent-core match."""
    if len(core_parent)!=len(core_candidate):raise ValueError("Indels require an explicit curated residue mapping")
    hits=[];start=full_parent.find(core_parent)
    while start>=0:
        hits.append(start);start=full_parent.find(core_parent,start+1)
    if offset is None:
        if len(hits)!=1:raise ValueError("Core must map uniquely to the confirmed full construct")
        offset=hits[0]
    elif offset not in hits:raise ValueError("Supplied construct offset does not match the parent core")
    sequence=full_parent[:offset]+core_candidate+full_parent[offset+len(core_parent):]
    mapping=[{"core_position":i+1,"construct_position":i+offset+1,"parent":a,"candidate":b}
             for i,(a,b) in enumerate(zip(core_parent,core_candidate)) if a!=b]
    return {"sequence":sequence,"sequence_sha256":hashlib.sha256(sequence.encode()).hexdigest(),
            "parent_sequence_sha256":hashlib.sha256(full_parent.encode()).hexdigest(),"mutations":mapping,
            "core_offset":offset,"original_prediction_preserved":True}

def freeze_batch(batch,output):
    batch=MechanismBatch.model_validate(batch)
    output=Path(output)
    payload=batch.model_dump_json(indent=2)
    if output.exists():
        if MechanismBatch.model_validate_json(output.read_text(encoding="utf-8"))!=batch:
            raise ValueError("Frozen batch already exists with different settings")
    else:
        output.parent.mkdir(parents=True,exist_ok=True);output.write_text(payload,encoding="utf-8")
    return {"path":str(output),"total_production_ns":len(batch.candidates)*batch.replicas*batch.production_ns,
            "replicas":len(batch.candidates)*batch.replicas,"sampling_status":"pilot_not_converged_mechanism"}
