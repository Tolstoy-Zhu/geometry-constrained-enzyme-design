"""Family-specific catalytic roles, structural mapping and protected positions."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
from .contracts import EnzymeFamilyProfile

def profile_dir():
    from ..config import project_root
    return project_root() / "data/families"

def profiles():
    return [EnzymeFamilyProfile.model_validate_json(p.read_text(encoding="utf-8"))
            for p in sorted(profile_dir().glob("*.json")) if p.stem != "sources"]

def get_profile(family="bur01"):
    if isinstance(family, EnzymeFamilyProfile):
        return family
    path=profile_dir()/(family+".json")
    if not path.is_file() or path.parent.resolve()!=profile_dir().resolve():
        raise ValueError("Unknown enzyme family; catalytic evidence is not evaluated: "+str(family))
    return EnzymeFamilyProfile.model_validate_json(path.read_text(encoding="utf-8"))

def reference(profile):
    from ..config import resolve_path
    path=resolve_path(profile.reference_pdb)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=profile.reference_sha256:
        raise ValueError("Family reference checksum mismatch: "+profile.family_id)
    return path

def assess_structure(trace, family="bur01"):
    from ..structure_search import load_ca_trace, compare_ca_traces
    profile=get_profile(family)
    ref=load_ca_trace(reference(profile),chain_id=profile.reference_chain)
    comparison=compare_ca_traces(ref,trace,catalytic_query_residues=tuple(
        (s.residue_number,s.amino_acid) for s in profile.sites))
    mapped=comparison["mapped_catalytic_residues"].split(";")
    labels={r.label:i for i,r in enumerate(trace.residues)}
    sites=[]
    for definition,label in zip(profile.sites,mapped):
        index=labels.get(label)
        residue=trace.residues[index] if index is not None else None
        missing=[a for a in definition.atoms if residue is None or a not in residue.atoms]
        sites.append({"role":definition.role,"reference_number":definition.residue_number,
                      "residue":label,"index":index,
                      "identity_pass":residue is not None and residue.amino_acid==definition.amino_acid,
                      "missing_atoms":missing})
    by_role={s["role"]:s for s in sites}
    distances=[]
    for rule in profile.geometry:
        a,b=by_role[rule.first_role],by_role[rule.second_role]
        value=None
        if a["identity_pass"] and b["identity_pass"] and not a["missing_atoms"] and not b["missing_atoms"]:
            da=next(s for s in profile.sites if s.role==a["role"])
            db=next(s for s in profile.sites if s.role==b["role"])
            xyz_a=np.array([trace.residues[a["index"]].atoms[n] for n in da.atoms])
            xyz_b=np.array([trace.residues[b["index"]].atoms[n] for n in db.atoms])
            value=float(np.linalg.norm(xyz_a[:,None]-xyz_b[None],axis=-1).min())
        distances.append({"roles":[a["role"],b["role"]],"distance_A":value,
                          "minimum_A":rule.minimum_A,"maximum_A":rule.maximum_A,
                          "status":"not_assessable" if value is None else
                                   ("compatible" if rule.minimum_A<=value<=rule.maximum_A else "out_of_range")})
    valid=len(sites)==len(profile.sites) and all(s["identity_pass"] for s in sites)
    protected=[]
    if profile.protected_residues:
        refs={r.residue_number:r for r in ref.residues}
        protected_comparison=compare_ca_traces(ref,trace,catalytic_query_residues=tuple(
            (number,refs[number].amino_acid) for number in profile.protected_residues))
        protected=[labels.get(label) for label in protected_comparison['mapped_catalytic_residues'].split(';')]
        if None in protected:raise ValueError('Family protected positions could not be mapped')
    return {"family_id":profile.family_id,"sites":sites,"geometry":distances,
            "protected_indices":protected,
            "fold_pass":min(comparison["qtmscore"],comparison["ttmscore"])>=profile.min_tm_score,
            "catalytic_identity_pass":valid,"qtmscore":comparison["qtmscore"],"ttmscore":comparison["ttmscore"],
            "validation_scope":profile.validation_scope}

def map_sequence(sequence,family="bur01",override=None):
    if family=="bur01":
        from ..prediction_sequence import map_8hfw_catalytic_sites
        return map_8hfw_catalytic_sites(sequence,override)
    from ..prediction_sequence import validate_sequence,_alignment_index_map
    from ..structure_search import load_ca_trace
    sequence=validate_sequence("query",sequence).sequence
    profile=get_profile(family)
    allowed={f'{s.amino_acid}{s.residue_number}' for s in profile.sites}
    if override and set(override)-allowed:raise ValueError('Override contains residues outside this family')
    ref=load_ca_trace(reference(profile),chain_id=profile.reference_chain)
    mapping,identity,score=_alignment_index_map(ref.sequence,sequence)
    positions={r.residue_number:i for i,r in enumerate(ref.residues)}
    entries={}; valid={}; missing=[]
    for site in profile.sites:
        key=f"{site.amino_acid}{site.residue_number}"
        idx=mapping.get(positions.get(site.residue_number))
        pos=(override or {}).get(key,idx+1 if idx is not None else 0)
        observed=sequence[pos-1] if 1<=pos<=len(sequence) else None
        ok=observed==site.amino_acid
        entries[key]={"reference_residue":key,"predicted_position":pos or None,"expected_amino_acid":site.amino_acid,
                      "observed_amino_acid":observed,"status":"mapped" if ok else "identity_mismatch",
                      "source":"manual" if key in (override or {}) else "global_alignment","role":site.role}
        if ok: valid[key]=pos
        else: missing.append(key)
    return {"family_id":family,"reference_pdb":str(reference(profile)),"reference_chain":profile.reference_chain,
            "sites":entries,"valid_positions":valid,"missing_evidence":missing,"aligned_identity":identity,
            "alignment_score":score,"complete":not missing,"reference_structure":profile.reference_pdb}

def scope_scores(frame,family):
    """WASE is never applied outside its defined family."""
    if family=="bur01":
        return frame
    out=frame.copy()
    for column in list(out.columns):
        if column.startswith("wase_") or column in {"structure_fit_score","catalytic_core_score",
                                                   "annotation_confidence_score","experimental_readiness_score"}:
            out[column]=None
    out["family_id"]=family
    out["wase_scope"]="not_applicable_to_family"
    out["recommendation"]="Review family-specific structural evidence"
    return out
