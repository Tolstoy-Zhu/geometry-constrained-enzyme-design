"""Read historical and durable-task artifacts without rewriting the source manifests."""
from __future__ import annotations
import json
from pathlib import Path
from .store import Store

def local_artifact(value,base):
    if not isinstance(value,str):return value
    base=Path(base).resolve()
    for parent in [base,*list(base.parents)[:5]]:
        receipt=parent/".aqua-download.json"
        if receipt.is_file():
            mapping=json.loads(receipt.read_text(encoding="utf-8"))
            if value in mapping:
                candidate=(parent/mapping[value]).resolve()
                if candidate.is_relative_to(parent) and candidate.is_file():return str(candidate)
            break
    return value

def result_roots(project_id="default", store=None):
    store=store or Store()
    roots=[]
    for row in store.tasks(project_id):
        if row.get("result"):
            result=json.loads(row["result"])
            downloaded=result.get("local_output")
            if downloaded and Path(downloaded).is_dir():
                roots.append(Path(downloaded))
        local=store.path.parent/"tasks"/row["id"]/"output"
        if local.is_dir(): roots.append(local)
    return list(dict.fromkeys(p.resolve() for p in roots))

def manifest_paths(name, legacy_root=None, project_id="default", store=None):
    store=store or Store(); paths=[]
    if project_id=="default" and legacy_root:
        paths.extend(Path(legacy_root).glob("*/"+name))
    for record in store.records("history",project_id):
        path=Path(record["path"])
        if path.name==name and path.is_file(): paths.append(path)
    for root in result_roots(project_id,store): paths.extend(root.rglob(name))
    return sorted(set(p.resolve() for p in paths))

def read_manifest(path, directory_key):
    path=Path(path).resolve()
    data=json.loads(path.read_text(encoding="utf-8"))
    # Preserve the remote/original manifest on disk. Only the viewer's copy is rebased.
    def rebase(value):
        if isinstance(value,dict): return {k:rebase(v) for k,v in value.items()}
        if isinstance(value,list): return [rebase(v) for v in value]
        if isinstance(value,str):
            resolved=local_artifact(value,path.parent)
            if resolved!=value:return resolved
            parts=value.replace("\\","/").split("/")
            for i,part in enumerate(parts):
                if part==path.parent.name and i<len(parts)-1:
                    rel=Path(*parts[i+1:])
                    candidate=(path.parent/rel).resolve()
                    if candidate.is_relative_to(path.parent) and candidate.is_file():
                        return str(candidate)
        return value
    result=rebase(data)
    result["_manifest_path"]=str(path)
    result[directory_key]=str(path.parent)
    return result

def scan_manifests(name, directory_key, legacy_root=None, project_id="default", store=None):
    rows=[]
    for path in manifest_paths(name,legacy_root,project_id,store):
        try: rows.append(read_manifest(path,directory_key))
        except (OSError,ValueError): continue
    return sorted(rows,key=lambda r:r.get("created_at",""),reverse=True)

def resolve_table_paths(frame,base):
    result=frame.copy()
    for column in result.columns:
        if column.endswith(("_path","_file")):
            result[column]=result[column].map(lambda value:local_artifact(value,base))
    return result
