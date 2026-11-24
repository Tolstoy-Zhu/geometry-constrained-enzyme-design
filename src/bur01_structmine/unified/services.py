"""Shared operations used by the CLI, persistent workers and the existing workbench."""
from __future__ import annotations
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from .contracts import EnzymeRecord, ReactionSpec, AssayRecord, TaskSpec
from .store import Store, uid
from .execution import submit, start_worker, server_profiles, RemoteExecutor

KINDS={"enzyme":EnzymeRecord,"reaction":ReactionSpec,"assay":AssayRecord}

def import_records(path,project_id="default",store=None):
    store=store or Store()
    data=json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if set(data)-set(KINDS):raise ValueError("Unknown record collection")
    validated={kind:[cls.model_validate(row) for row in data.get(kind,[])] for kind,cls in KINDS.items()}
    enzymes={r["enzyme_id"] for r in store.records("enzyme",project_id)}
    reactions={r["reaction_id"] for r in store.records("reaction",project_id)}
    enzymes.update(r.enzyme_id for r in validated["enzyme"])
    reactions.update(r.reaction_id for r in validated["reaction"])
    for row in validated["assay"]:
        if row.enzyme_id not in enzymes or row.reaction_id not in reactions:
            raise ValueError("An assay references an unknown enzyme or reaction")
    # Validate all conflicts before writing any rows; one transaction prevents partial imports.
    with store.connect() as db:
        for kind,rows in validated.items():
            for row in rows:
                identifier=getattr(row,{"enzyme":"enzyme_id","reaction":"reaction_id","assay":"assay_id"}[kind])
                from .store import dumps, now
                payload=dumps(row.model_dump(mode="json"))
                prior=db.execute("SELECT payload FROM records WHERE kind=? AND id=? AND project_id=?",(kind,identifier,project_id)).fetchone()
                if prior and prior[0]!=payload:raise ValueError("Conflicting record: "+identifier)
                db.execute("INSERT OR IGNORE INTO records VALUES(?,?,?,?,?)",(kind,identifier,project_id,payload,now()))
    return {kind:len(rows) for kind,rows in validated.items()}

def export_assays(output,project_id="default",store=None):
    import pandas as pd
    store=store or Store()
    enzymes={r["enzyme_id"]:r for r in store.records("enzyme",project_id)}
    reactions={r["reaction_id"]:r for r in store.records("reaction",project_id)}
    rows=[]
    for assay in store.records("assay",project_id):
        enzyme=enzymes[assay["enzyme_id"]];reaction=reactions[assay["reaction_id"]]
        substrates=[m for m in reaction["molecules"] if m["role"]=="substrate"]
        # Multi-substrate kinetics needs explicit parameter-substrate attribution.
        smiles=assay["conditions"].get("parameter_substrate_smiles")
        if not smiles and len(substrates)==1:smiles=substrates[0]["smiles"]
        rows.append({**assay,"record_id":assay["assay_id"],"sequence":enzyme["sequence"],
                     "parent_id":enzyme["parent_id"] or enzyme["enzyme_id"],"smiles":smiles,
                     "temperature_K":reaction["temperature_K"],"pH":reaction["pH"],
                     "construct_id":enzyme["construct_id"]})
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(output,index=False)
    return str(output)

def extract_features(dataset,output,device="cpu"):
    import pandas as pd
    from ..property_assets import asset_root,runtime_python
    path=Path(dataset);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    frame=pd.read_csv(path/"measurements.csv")
    sequences=dict(zip(frame.sequence_sha256,frame.sequence))
    input_path=output/"sequences.json"
    input_path.write_text(json.dumps(sequences),encoding="utf-8")
    interpreter=os.environ.get("AQUA_FEATURE_PYTHON",str(runtime_python()))
    if not Path(interpreter).is_file():interpreter=sys.executable
    script=Path(__file__).with_name("features_worker.py")
    subprocess.run([interpreter,str(script),str(input_path),str(output),str(asset_root()/"esm"),
                    "--device",device],check=True)
    return {"output":str(output),"manifest":str(output/"manifest.json")}

def reconcile_task(task_id,store=None):
    store=store or Store()
    row=store.task(task_id);spec=TaskSpec.model_validate_json(row["spec"])
    if spec.executor!="local":return RemoteExecutor(server_profiles()[spec.executor]).reconcile(task_id,store)
    start_worker(store)
    return store.task(task_id)

def cancel_task(task_id,project_id,store=None):
    store=store or Store();row=store.task(task_id)
    if row["project_id"]!=project_id:raise ValueError("Task belongs to a different project")
    spec=TaskSpec.model_validate_json(row["spec"])
    if spec.executor!="local":RemoteExecutor(server_profiles()[spec.executor]).cancel(task_id)
    store.cancel(task_id)

def retry_task(task_id,store=None):
    store=store or Store();row=store.task(task_id);spec=TaskSpec.model_validate_json(row["spec"])
    if spec.executor!="local":
        # Never retry a remotely unconfirmed task under a new identifier.
        RemoteExecutor(server_profiles()[spec.executor]).reconcile(task_id,store)
        row=store.task(task_id)
    if row["state"]=="blocked":raise ValueError("Reconcile the original task before creating a retry")
    if spec.executor!="local":
        remote=RemoteExecutor(server_profiles()[spec.executor])
        result=json.loads(remote.command("retry","--root",remote.profile.remote_root,"--task-id",task_id))
        new_id=result["task_id"]
        if new_id==task_id:store.retry(task_id)
        else:
            spec.parent_task_id=spec.task_id;spec.task_id=new_id;store.submit(spec)
        return new_id
    new_id=store.retry(task_id)
    start_worker(store)
    return new_id

def property_reports(store=None,project_id="default"):
    store=store or Store();reports={}
    for row in store.tasks(project_id):
        spec=json.loads(row["spec"])
        if spec["kind"]!="properties" or not row.get("result"):continue
        result=json.loads(row["result"])
        root=Path(result["local_output"])/"properties" if result.get("local_output") else Path(result.get("output",""))
        if (root/"index.html").is_file():reports[row["id"]]=root
    return reports


def kinetic_reports(store=None,project_id="default"):
    store=store or Store();reports={}
    for row in store.tasks(project_id):
        if json.loads(row['spec'])['kind']!='kinetics' or row['state']!='completed' or not row.get('result'):continue
        result=json.loads(row['result'])
        root=Path(result['local_output'])/'kinetics' if result.get('local_output') else Path(result.get('output',''))
        if (root/'predictions.json').is_file():reports[row['id']]=root
    return reports

class DurablePropertyJobs:
    def context(self):
        import streamlit as st
        return {"project_id":st.session_state.get("aqua_project","default"),
                "family_id":st.session_state.get("aqua_family","bur01"),
                "executor":st.session_state.get("aqua_executor","local"),
                "reaction_id":st.session_state.get("aqua_reaction") or None}
    def submit(self,title,fn,**inputs):
        mode="fasta" if fn.__name__=="profile_fasta" else "runs"
        inputs.setdefault("device","cuda" if self.context()["executor"]!="local" else "cpu")
        return submit("properties",{"mode":mode,"title":title,**inputs},**self.context())
    def snapshot(self):
        from .polling import schedule_sync
        schedule_sync(self.context()["project_id"])
        items=[]
        for row in Store().tasks(self.context()["project_id"]):
            spec=json.loads(row["spec"])
            if spec["kind"]!="properties":continue
            items.append(SimpleNamespace(key=row["id"],state="complete" if row["state"]=="completed" else row["state"],
                         message=row["message"] or "已保存，等待计算",title=spec["inputs"].get("title","性质预测")))
        return items

def register_model(model_dir,project_id="default",store=None):
    from .datasets import file_hash
    store=store or Store();folder=Path(model_dir)
    card=json.loads((folder/"model_card.json").read_text(encoding="utf-8"))
    if file_hash(folder/"model.joblib")!=card["model_sha256"]:raise ValueError("Model checksum mismatch")
    return _register_card(card,str(folder.resolve()),"local",project_id,store)


def register_task_model(task_id,project_id="default",store=None):
    store=store or Store();row=store.task(task_id)
    if row["project_id"]!=project_id:raise ValueError("Training belongs to a different project")
    spec=TaskSpec.model_validate_json(row["spec"])
    if spec.kind!="train" or row["state"]!="completed":raise ValueError("Select a completed training task")
    if spec.executor=="local":return register_model(json.loads(row["result"])["output"],project_id,store)
    remote=RemoteExecutor(server_profiles()[spec.executor])
    result=json.loads(remote.command("model-card","--root",remote.profile.remote_root,"--task-id",task_id))
    return _register_card(result["card"],"remote:"+result["path"],spec.executor,project_id,store)


def _register_card(card,path,executor,project_id,store):
    store.record("model",card["model_id"],project_id,{"path":path,"executor":executor,"card":card})
    if card["status"]=="validated_upgrade" and card["promotion"].get("promote"):
        # Explicit registration is the promotion step, separate from training.
        key="default_model_"+card["endpoint"]
        with store.connect() as db:
            db.execute("INSERT OR REPLACE INTO metadata VALUES(?,?)",(project_id+":"+key,card["model_id"]))
    return {"model_id":card["model_id"],"status":card["status"]}
