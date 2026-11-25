"""Detached local workers and SSH transport with immutable, idempotent task IDs."""
from __future__ import annotations
import io
import json
import os
import re
import shlex
import subprocess
import sys
import tarfile
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from pydantic import Field, field_validator
from .contracts import Contract, TaskSpec
from .store import Store, state_root, uid

class ServerProfile(Contract):
    name: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")
    host: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9.-]*$")
    port: int = Field(default=22,ge=1,le=65535)
    user: str = Field(default="root",pattern=r"^[A-Za-z_][A-Za-z0-9_.-]*$")
    identity_file: str | None = None
    remote_root: str = "/root/autodl-tmp/aquaesternet"
    python: str = "/root/miniconda3/bin/python"
    @field_validator("remote_root","python")
    @classmethod
    def absolute_posix(cls,value):
        p=PurePosixPath(value)
        if not p.is_absolute() or ".." in p.parts or any(c in value for c in "\r\n\0"):
            raise ValueError("Use an absolute POSIX path without parent traversal")
        return value

def server_profiles():
    path=state_root()/"servers.json"
    if not path.exists(): return {}
    rows=json.loads(path.read_text(encoding="utf-8"))
    return {p.name:p for p in (ServerProfile.model_validate(r) for r in rows)}

def save_server(profile):
    values=server_profiles();values[profile.name]=profile
    path=state_root()/"servers.json";path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps([p.model_dump() for p in values.values()],indent=2),encoding="utf-8")
    os.replace(tmp,path)

@contextmanager
def worker_lock(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    handle=(directory/"worker.lock").open("a+b")
    try:
        handle.seek(0)
        if os.name=="nt":
            import msvcrt
            if handle.read(1)==b"": handle.write(b"0");handle.flush()
            handle.seek(0);msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except (OSError,BlockingIOError):
        handle.close()
        raise RuntimeError("Another worker owns this resource")
    try:
        yield
    finally:
        handle.close()

def start_worker(store=None):
    store=store or Store()
    log=store.path.parent/"worker.log"
    env=os.environ.copy()
    env["PYTHONIOENCODING"]="utf-8"
    env["AQUA_STATE_ROOT"]=str(store.path.parent)
    # Ensure this checkout, not an unrelated editable install, is used.
    env["PYTHONPATH"]=str(Path(__file__).resolve().parents[2])
    options={"start_new_session":True} if os.name!="nt" else {
        "creationflags":subprocess.CREATE_NO_WINDOW|subprocess.DETACHED_PROCESS|subprocess.CREATE_NEW_PROCESS_GROUP}
    with log.open("ab") as out:
        subprocess.Popen([sys.executable,"-m","bur01_structmine.unified.worker","consume",
                          "--database",str(store.path)],stdout=out,stderr=out,stdin=subprocess.DEVNULL,env=env,**options)

def submit(kind,inputs,*,project_id="default",family_id="bur01",executor="local",store=None,reaction_config=None,reaction_id=None):
    store=store or Store()
    from .contracts import ReactionSpec
    reaction=None
    if reaction_config:
        reaction=ReactionSpec.model_validate_json(Path(reaction_config).read_text(encoding="utf-8-sig"))
        store.record("reaction",reaction.reaction_id,project_id,reaction.model_dump(mode="json"))
    elif reaction_id:
        found=next((r for r in store.records("reaction",project_id) if r["reaction_id"]==reaction_id),None)
        if found is None:raise ValueError("Unknown project reaction")
        reaction=ReactionSpec.model_validate(found)
    spec=TaskSpec(task_id=uid("task"),project_id=project_id,kind=kind,family_id=family_id,
                  executor=executor,inputs=json.loads(json.dumps(inputs,default=str)),reaction=reaction)
    if kind=='mechanism' and inputs.get('action')=='assess':spec.resource='cpu'
    if executor=="local":
        # Uploads and selected input files remain valid after the web request ends.
        from .worker import restore_paths
        base=store.path.parent/'tasks'/spec.task_id/'input'
        safe_extract(task_bundle(spec),base)
        raw=json.loads((base/'task.json').read_text(encoding='utf-8'))
        spec.inputs=restore_paths(raw['inputs'],base)
        store.submit(spec)
        start_worker(store)
    else:
        store.submit(spec)
        try: RemoteExecutor(server_profiles()[executor]).submit(spec)
        except Exception as exc:
            store.update(spec.task_id,"blocked","远端提交未确认；重连时先核对同一任务编号："+str(exc))
    return spec.task_id

PATH_KEYS={"parent","fasta_path","local_msa_path","parameter_dir","config_path","reaction_config","dataset","features","cluster_file",
           "prepared_manifest","batch_manifest","msa","structure","run","input","protocol","reference",
           'query_structure','input_tsv','topology','trajectories','channel_manifest','model_dir','source_run_dir','bundle','query_fasta',
           'assembly_report','protonation_reference','system_dir','prerequisite_dir','quality_report'}
def task_bundle(spec):
    """Copy only task-selected inputs. Existing result trees and private projects are never swept."""
    files={};blobs={};copied={};counter=0
    def copy_file(source,key=""):
        nonlocal counter
        if Path(source).is_symlink():raise ValueError("Task inputs must not be symbolic links")
        source=Path(source).resolve()
        if source in copied:return copied[source]
        if not source.exists():raise FileNotFoundError(source)
        if source.is_symlink():raise ValueError("Task inputs must not be symbolic links")
        name=f"inputs/{counter:04d}/{source.name}"
        counter+=1;copied[source]=name
        if source.is_dir():
            if not any(source.rglob("*")):raise ValueError("Selected input directory is empty")
            for p in sorted(source.rglob("*")):
                if p.is_symlink():raise ValueError("Input directories cannot contain symbolic links")
                if p.is_file():files[name+"/"+p.relative_to(source).as_posix()]=p
        else:
            files[name]=source
            if key=='local_msa_path':
                evidence=source.with_suffix('.provenance.json')
                if evidence.is_file():files[str(PurePosixPath(name).with_suffix('.provenance.json'))]=evidence
            if key=="cluster_file":
                evidence=source.with_suffix(".json")
                if not evidence.is_file():raise ValueError("Cluster mapping needs its provenance JSON")
                files[str(PurePosixPath(name).with_suffix(".json"))]=evidence
            if key in {"prepared_manifest","batch_manifest","channel_manifest"} or (key=="manifest_dependency" and source.suffix==".json"):
                from .transfer import manifest_dependencies,receipt_bytes
                deps=[(raw,copy_file(path,"manifest_dependency")) for raw,path in manifest_dependencies(source)]
                if deps:blobs[str(PurePosixPath(name).parent/".aqua-transfer.json")]=receipt_bytes(name,deps)
        return name
    def visit(value,key=""):
        nonlocal counter
        if isinstance(value,dict): return {k:visit(v,k) for k,v in value.items()}
        if isinstance(value,list): return [visit(v,"run" if key=="directories" else key) for v in value]
        if isinstance(value,str) and key in PATH_KEYS:
            if value.startswith("remote:"):
                if spec.executor=="local": raise ValueError("Remote inputs require a remote executor")
                return {"__aqua_remote_path__":value[7:]}
            if spec.executor=='local' and spec.kind=='mechanism' and spec.inputs.get('action')=='assess' and key in {'run','batch_manifest','quality_report'}:
                # Existing MD is read-only. The assessment verifies its frozen source
                # hashes before analysis; staging a multi-GB trajectory tree is wasteful.
                path=Path(value)
                if path.is_symlink() or not path.exists():raise ValueError('Missing or symbolic analysis source')
                if key=='run' and not path.is_dir():raise ValueError('Production directory required')
                if key!='run' and not path.is_file():raise ValueError('Analysis manifest file required')
                return str(path.resolve())
            return {"__aqua_path__":copy_file(value,key)}
        return value
    payload=spec.model_dump();payload["inputs"]=visit(spec.inputs)
    if sum(p.stat().st_size for p in files.values())>2_000_000_000:
        raise ValueError("Task inputs exceed 2 GB; use a remote reference for large stored data")
    out=io.BytesIO()
    with tarfile.open(fileobj=out,mode="w:gz") as tar:
        blob=json.dumps(payload,ensure_ascii=False).encode()
        info=tarfile.TarInfo("task.json");info.size=len(blob);tar.addfile(info,io.BytesIO(blob))
        for name,blob in blobs.items():
            info=tarfile.TarInfo(name);info.size=len(blob);tar.addfile(info,io.BytesIO(blob))
        for name,path in files.items():
            tar.add(path,arcname=name,recursive=False)
    return out.getvalue()

def safe_extract(data,destination):
    destination=Path(destination).resolve();destination.mkdir(parents=True,exist_ok=True)
    total=0
    with tarfile.open(fileobj=io.BytesIO(data),mode="r:gz") as tar:
        for member in tar:
            rel=PurePosixPath(member.name)
            target=destination.joinpath(*rel.parts).resolve()
            if (rel.is_absolute() or ".." in rel.parts or "\\" in member.name or ":" in member.name
                    or not target.is_relative_to(destination) or not member.isfile()):
                raise ValueError("Unsafe task archive entry")
            if os.name=='nt' and not str(target).startswith('\\\\?\\'):
                # Extended paths avoid the Windows 260-character legacy limit;
                # all containment checks above still use the ordinary resolved path.
                target=Path('\\\\?\\'+str(target))
            total+=member.size
            if total>2_000_000_000: raise ValueError("Task archive exceeds 2 GB")
            if target.exists():
                with tar.extractfile(member) as source:
                    content=source.read()
                if target.read_bytes()!=content: raise ValueError("Existing task input differs")
            else:
                target.parent.mkdir(parents=True,exist_ok=True)
                with tar.extractfile(member) as source:
                    target.write_bytes(source.read())

class RemoteExecutor:
    def __init__(self,profile): self.profile=profile

    def command(self,*args,stdin=None,timeout=90):
        p=self.profile
        command=["ssh","-T","-o","BatchMode=yes","-o","StrictHostKeyChecking=yes",
                 "-o","ConnectTimeout=15","-p",str(p.port)]
        if p.identity_file: command+=["-i",str(Path(p.identity_file).expanduser())]
        command += [f"{p.user}@{p.host}",shlex.join([p.python,"-m","bur01_structmine.unified.worker",*args])]
        result=subprocess.run(command,input=stdin,capture_output=True,timeout=timeout)
        if result.returncode:
            raise RuntimeError(result.stderr.decode("utf-8",errors="replace")[-2000:])
        return result.stdout

    def check(self):
        return json.loads(self.command("probe"))

    def submit(self,spec):
        return json.loads(self.command("ingest","--root",self.profile.remote_root,stdin=task_bundle(spec),timeout=300))

    def reconcile(self,task_id,store=None):
        store=store or Store()
        result=json.loads(self.command('exists','--root',self.profile.remote_root,'--task-id',task_id))
        if result['exists']:
            return self.sync(task_id,store)
        spec=TaskSpec.model_validate_json(store.task(task_id)['spec'])
        self.submit(spec)
        return self.sync(task_id,store,download=False)

    def sync(self,task_id,store=None,download=True):
        store=store or Store()
        # Validate locally before incorporating the ID in a remote command.
        previous=store.task(task_id)
        TaskSpec.model_validate_json(previous["spec"])
        prior_result=json.loads(previous["result"]) if previous.get("result") else {}
        row=json.loads(self.command("status","--root",self.profile.remote_root,"--task-id",task_id))
        log_dir=store.path.parent/"tasks"/task_id
        log_dir.mkdir(parents=True,exist_ok=True)
        if "log_tail" in row:
            (log_dir/"process.log").write_text(row["log_tail"],encoding="utf-8")
        remote_result=json.loads(row["result"]) if row.get("result") else None
        result={"remote":remote_result,"remote_updated":row["updated"]}
        if row.get('mechanism_progress'):result['progress']=row['mechanism_progress']
        if prior_result.get("remote_updated")==row["updated"] and prior_result.get("local_output"):
            if (Path(prior_result["local_output"])/".aqua-download.json").is_file():
                result["local_output"]=prior_result["local_output"]
        # Persist authoritative remote state even when transfer is interrupted.
        store.update(task_id,row["state"],row["message"],result)
        if download and "local_output" not in result and row["state"] in {"completed","failed","cancelled","interrupted"}:
            import hashlib
            snapshot=hashlib.sha256(row["updated"].encode()).hexdigest()[:16]
            dest=store.path.parent/"tasks"/task_id/"download"/snapshot
            try:
                data=self.command("export","--root",self.profile.remote_root,"--task-id",task_id,timeout=300)
                safe_extract(data,dest)
                result["local_output"]=str(dest)
            except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as exc:
                result["download_error"]=str(exc)
        store.update(task_id,row["state"],row["message"],result)
        return store.task(task_id)

    def cancel(self,task_id):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}",task_id): raise ValueError("Invalid task")
        return json.loads(self.command("cancel","--root",self.profile.remote_root,"--task-id",task_id))
