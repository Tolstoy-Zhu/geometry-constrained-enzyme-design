"""Worker entrypoint. Job kinds are explicit; no arbitrary shell command execution."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import os
import signal
import subprocess
import sys
import tarfile
import time
import threading
from pathlib import Path
from .contracts import TaskSpec
from .store import Store,state_root
from .execution import worker_lock,safe_extract,start_worker

def task_output(spec,store):
    path=store.path.parent/"tasks"/spec.task_id/"output"
    path.mkdir(parents=True,exist_ok=True)
    return path

def run_task(spec,store):
    inputs=dict(spec.inputs);out=task_output(spec,store)
    if spec.reaction:
        (out/"reaction.json").write_text(spec.reaction.model_dump_json(indent=2),encoding="utf-8")
        store.record("reaction",spec.reaction.reaction_id,spec.project_id,spec.reaction.model_dump(mode="json"))
    if spec.kind=="properties":
        from ..property_screen import profile_fasta,profile_runs
        mode=inputs.pop("mode","fasta")
        inputs.pop("title",None)
        fn=profile_fasta if mode=="fasta" else profile_runs
        result=fn(**inputs,output=out/"properties")
        return {"output":str(result),"report":str(Path(result)/"index.html")}
    if spec.kind=="design":
        from ..design import run_design
        inputs.setdefault("family_id",spec.family_id)
        result=run_design(**inputs,output=out/"design")
        return {"output":str(result)}
    if spec.kind=="predict":
        from ..models import PredictionRequest
        from ..prediction import run_prediction,run_prediction_and_mine
        mine_after=inputs.pop('mine_after',False)
        if os.name!='nt':
            inputs['runtime_backend']='native'
            inputs['colabfold_executable']=os.environ.get('AQUA_COLABFOLD_EXECUTABLE','colabfold_batch')
        inputs["family_id"]=spec.family_id;inputs["output_root"]=str(out/"predictions")
        if mine_after:
            prediction,mining=run_prediction_and_mine(PredictionRequest.model_validate(inputs),
                          mining_output_root=str(out/'runs'))
            return {'prediction':prediction.model_dump(mode='json'),'mining':mining.model_dump(mode='json') if mining else None}
        return run_prediction(PredictionRequest.model_validate(inputs)).model_dump(mode="json")
    if spec.kind=="mine":
        from ..pipeline import run_pipeline
        inputs["family_id"]=spec.family_id;inputs["output_root"]=str(out/"runs")
        return run_pipeline(**inputs).model_dump(mode="json")
    if spec.kind=="nr":
        from ..cli import main as legacy_main
        command=list(inputs["arguments"])
        allowed={"nr-expand","nr-export-job","nr-import-job","nr-demo-build","nr-novelty"}
        if not command or command[0] not in allowed:raise ValueError("Unsupported NR operation")
        for flag,key in {"--config":"config_path","--bundle":"bundle","--query-fasta":"query_fasta","--fasta":"fasta_path"}.items():
            if key in inputs:
                if flag in command:command[command.index(flag)+1]=inputs[key]
                else:command += [flag,inputs[key]]
        if "source_run_dir" in inputs:
            command += ["--structure-run-root",str(Path(inputs["source_run_dir"]).parent)]
        if command[0] in {"nr-expand","nr-export-job"}:
            command += ["--family",spec.family_id]
            if os.name!="nt" and "--config" in command:
                import yaml
                index=command.index("--config")+1
                cfg=yaml.safe_load(Path(command[index]).read_text(encoding="utf-8"))
                cfg.setdefault("runtime",{})["backend"]="native"
                config=out/"nr.runtime.yaml"
                config.write_text(yaml.safe_dump(cfg,allow_unicode=True),encoding="utf-8")
                command[index]=str(config)
        if command[0]!="nr-demo-build":command += ["--output-root",str(out/"nr_runs")]
        code=legacy_main(command)
        if code:raise RuntimeError("NR operation failed with exit code "+str(code))
        return {"output":str(out),"operation":command[0]}
    if spec.kind=="index":
        from .store import index_history
        return {"indexed":index_history(inputs["root"],store,spec.project_id)}
    if spec.kind=="dataset":
        from .datasets import build_dataset
        return build_dataset(output=out/"dataset",**inputs)
    if spec.kind=="train":
        from .learning import train_models
        return train_models(output=out/"model",**inputs)
    if spec.kind=="features":
        from .services import extract_features
        return extract_features(output=out/'features',**inputs)
    if spec.kind=="kinetics":
        from .inference import predict_kinetics
        return predict_kinetics(output=out/'kinetics',**inputs)
    if spec.kind=="mechanism":
        from .mechanism import execute_mechanism
        return execute_mechanism(output=out/"mechanism",**inputs)
    raise ValueError("Unsupported task kind")

def consume(store):
    try:
        with worker_lock(store.path.parent):
            store.recover()
            idle=0
            while idle<10:
                spec=store.claim()
                if spec is None:
                    idle+=1;time.sleep(1);continue
                idle=0
                directory=store.path.parent/"tasks"/spec.task_id
                directory.mkdir(parents=True,exist_ok=True)
                env=os.environ.copy();env["PYTHONIOENCODING"]="utf-8"
                with (directory/"process.log").open("ab") as log:
                    options={"start_new_session":True} if os.name!="nt" else {"creationflags":subprocess.CREATE_NO_WINDOW}
                    interpreter=os.environ.get('AQUA_LEARNING_PYTHON',sys.executable) if spec.kind in {'train','dataset','kinetics'} else sys.executable
                    if spec.kind=='mechanism':interpreter=os.environ.get('AQUA_MECHANISM_PYTHON',sys.executable)
                    process=subprocess.Popen([interpreter,"-m",__package__+".worker","execute",
                                              "--database",str(store.path),"--task-id",spec.task_id],
                                             stdout=log,stderr=log,env=env,**options)
                    while process.poll() is None:
                        if store.task(spec.task_id)["cancel_requested"]:
                            # This PID is an owned child, never a discovered machine-wide process.
                            if os.name!="nt": os.killpg(process.pid,signal.SIGTERM)
                            else: subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],capture_output=True)
                            try: process.wait(timeout=15)
                            except subprocess.TimeoutExpired: process.kill();process.wait()
                            store.update(spec.task_id,"cancelled","已停止本任务")
                            break
                        time.sleep(.5)
                row=store.task(spec.task_id)
                if row["state"]=="running":
                    store.update(spec.task_id,"failed",f"计算进程退出，代码 {process.returncode}；见任务日志")
    except RuntimeError as exc:
        if "Another worker" not in str(exc): raise

def restore_paths(value,base):
    if isinstance(value,dict):
        if set(value)=={"__aqua_remote_path__"}:
            # References are confined to this deployment's research workspace.
            root=base.resolve().parents[3]
            path=Path(value["__aqua_remote_path__"]).resolve()
            if not path.is_relative_to(root) or not path.exists() or path.is_symlink():
                raise ValueError("Remote input is missing or outside the deployment directory")
            return str(path)
        if set(value)=={"__aqua_path__"}:
            path=(base/value["__aqua_path__"]).resolve()
            if not path.is_relative_to(base.resolve()): raise ValueError("Invalid transferred path")
            return str(path)
        return {k:restore_paths(v,base) for k,v in value.items()}
    if isinstance(value,list): return [restore_paths(v,base) for v in value]
    return value

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("action",choices=["consume","execute","probe","ingest","status","cancel","export","exists","retry","model-card"])
    parser.add_argument("--database");parser.add_argument("--root");parser.add_argument("--task-id")
    args=parser.parse_args()
    if args.action=="probe":
        import platform,importlib.util,shutil
        print(json.dumps({"python":platform.python_version(),"platform":platform.system(),
                          "package":str(Path(__file__).resolve()),
                          "modules":{m:importlib.util.find_spec(m) is not None for m in ["torch","openmm","rdkit","sklearn"]},
                          "executables":{m:shutil.which(m) for m in ["foldseek","mmseqs","colabfold_batch"]}}));return
    if args.root: os.environ["AQUA_STATE_ROOT"]=str(Path(args.root)/"state")
    store=Store(args.database)
    if args.action=="consume": consume(store);return
    if args.action=="ingest":
        import tempfile
        blob=sys.stdin.buffer.read(2_000_000_001)
        if len(blob)>2_000_000_000: raise ValueError("Archive too large")
        # First extraction validates the identifier before selecting the persistent directory.
        with tempfile.TemporaryDirectory(dir=store.path.parent) as temp:
            safe_extract(blob,temp)
            raw=json.loads((Path(temp)/"task.json").read_text(encoding="utf-8"))
            # A marker-free temporary version validates identity and task kind.
            spec=TaskSpec.model_validate(raw)
            base=store.path.parent/"tasks"/spec.task_id/"input"
            safe_extract(blob,base)
            spec.inputs=restore_paths(raw["inputs"],base)
            spec.executor="local"
            store.submit(spec)
        start_worker(store)
        print(json.dumps({"task_id":spec.task_id,"state":store.task(spec.task_id)["state"]}));return
    if not args.task_id: raise ValueError("task-id is required")
    if args.action=='exists':
        try:store.task(args.task_id);exists=True
        except ValueError:exists=False
        print(json.dumps({'exists':exists}));return
    row=store.task(args.task_id)
    if args.action=="retry":
        new_id=store.retry(args.task_id);start_worker(store)
        print(json.dumps({"task_id":new_id}));return
    if args.action=="model-card":
        spec=TaskSpec.model_validate_json(row["spec"])
        if spec.kind!="train" or row["state"]!="completed":raise ValueError("Select a completed training task")
        folder=store.path.parent/"tasks"/args.task_id/"output/model"
        from .datasets import file_hash
        card=json.loads((folder/"model_card.json").read_text(encoding="utf-8"))
        if file_hash(folder/"model.joblib")!=card["model_sha256"]:raise ValueError("Model checksum mismatch")
        print(json.dumps({"path":str(folder),"card":card}));return
    if args.action=="execute":
        spec=TaskSpec.model_validate_json(row["spec"])
        stopped=threading.Event()
        def monitor_cancel():
            while not stopped.wait(.5):
                if store.task(spec.task_id)['cancel_requested']:
                    store.update(spec.task_id,'cancelled','已取消本任务及其子进程')
                    if os.name=='nt':subprocess.run(['taskkill','/PID',str(os.getpid()),'/T','/F'],capture_output=True)
                    else:os.killpg(os.getpgrp(),signal.SIGTERM)
        try:
            with worker_lock(store.path.parent/'tasks'/spec.task_id), worker_lock(store.path.parent/'resources'):
                if row['cancel_requested']:
                    store.update(spec.task_id,'cancelled','执行前已取消');return
                threading.Thread(target=monitor_cancel,daemon=True).start()
                result=run_task(spec,store)
                if not store.task(spec.task_id)['cancel_requested']:
                    store.update(spec.task_id,"completed","任务已完成",result)
        except Exception as exc:
            if store.task(spec.task_id)['state']!='cancelled':store.update(spec.task_id,"failed",str(exc))
            raise
        finally:
            stopped.set()
            start_worker(store)
    elif args.action=="status":
        spec=TaskSpec.model_validate_json(row['spec'])
        if spec.kind=='mechanism' and spec.inputs.get('action')=='equilibrate':
            progress=store.path.parent/'tasks'/args.task_id/'output/mechanism/progress.json'
            if progress.is_file():
                data=json.loads(progress.read_text())
                row['mechanism_progress']={'completed_ps':data.get('completed_ps',0),
                    'completed_chunks':data.get('completed_chunks',0),
                    'last_sample':data.get('chunks',[{}])[-1].get('last_sample') if data.get('chunks') else None}
        log=store.path.parent/"tasks"/args.task_id/"process.log"
        if log.exists():
            with log.open("rb") as handle:
                handle.seek(max(0,log.stat().st_size-64000))
                row["log_tail"]=handle.read().decode("utf-8",errors="replace")
        print(json.dumps(row,ensure_ascii=False))
    elif args.action=="cancel":
        store.cancel(args.task_id);print(json.dumps({"cancel_requested":True}))
    elif args.action=="export":
        from .export import export_task
        export_task(store,args.task_id,sys.stdout.buffer)

if __name__=="__main__": main()
