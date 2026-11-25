"""Sequential, budget-bound production replicas and per-replica analysis."""
from __future__ import annotations
import json
import hashlib
from pathlib import Path
from .mechanism_plan import MechanismBatch
from .mechanism import checked_file,inspect_prepared,run_prepared,atomic_json
from .mechanism_analysis import analyze_trajectory
from .datasets import file_hash

def run_batch(batch_manifest,output,platform_name="CUDA"):
    source=Path(batch_manifest)
    batch=MechanismBatch.model_validate_json(source.read_text(encoding="utf-8"))
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    receipt=output/"batch_identity.json"
    identity={"manifest_sha256":file_hash(source),"budget_ns":batch.approved_budget_ns}
    if receipt.exists() and json.loads(receipt.read_text())!=identity:
        raise ValueError("This output belongs to a different frozen mechanism batch")
    prepared=[]
    for row in batch.candidates:
        path=checked_file(row["prepared"],source.parent)
        check=inspect_prepared(path)
        if check["status"]!="passed":
            raise ValueError(row["candidate_id"]+": "+", ".join(check["missing"]))
        if check["contract"].get("family_id")!=batch.family_id:
            raise ValueError("Prepared family does not match the frozen batch")
        sequences=check["contract"]["construct_sequences"]
        expected=row.get("construct_sequence_sha256")
        if not expected or any(hashlib.sha256(s.encode()).hexdigest()!=expected for s in sequences.values()):
            raise ValueError("Prepared construct does not match the frozen candidate sequence")
        from .contracts import ReactionSpec
        reaction=ReactionSpec.model_validate_json(Path(check["files"]["reaction"]).read_text(encoding="utf-8"))
        if reaction.reaction_id!=batch.reaction_id:
            raise ValueError("Prepared reaction does not match the frozen batch")
        prepared.append((row,path,check))
    atomic_json(receipt,identity)
    reports=[]
    for index,(row,path,check) in enumerate(prepared):
        for repeat in range(batch.replicas):
            replica=f"{row['candidate_id']}_r{repeat+1:02d}"
            folder=output/replica
            result=run_prepared(path,folder,production_ns=batch.production_ns,segment_ns=batch.segment_ns,
                seed=batch.seeds[index*batch.replicas+repeat],platform_name=platform_name,
                approved_budget_ns=batch.production_ns,replica_id=replica)
            analyses={}
            for region,settings in row.get("analysis_regions",{}).items():
                if not region.isidentifier():raise ValueError("Invalid active-site analysis name")
                analyses[region]=analyze_trajectory(check["files"]["prmtop"],
                    [s["path"] for s in result["segments"]],folder/("analysis_"+region),
                    replica_id=replica,**settings)
            reports.append({"candidate_id":row["candidate_id"],"replica_id":replica,"expected_regions":list(row.get('analysis_regions',{})),
                            "production_ns":batch.production_ns,"trajectory":result,"regions":analyses})
            atomic_json(output/"batch_progress.json",{"replicas_completed":len(reports),
                "production_ns_completed":sum(r["production_ns"] for r in reports),"reports":reports})
    result={"status":"production_completed" if all(r["regions"] for r in reports) else "analysis_required",
            "replicas_completed":len(reports),"production_ns":sum(r["production_ns"] for r in reports),
            "sampling_assessment":"pilot; MSM/TPT convergence not established","reports":reports}
    atomic_json(output/"batch_report.json",result)
    return {"output":str(output),**{k:v for k,v in result.items() if k!="reports"}}
