import io,json,tarfile,hashlib

def test_remote_sync_preserves_state_and_retries_interrupted_download(tmp_path, monkeypatch):
    from bur01_structmine.unified.execution import RemoteExecutor, ServerProfile
    store=Store(tmp_path/"client/projects.sqlite")
    spec=TaskSpec(task_id="sync-test",project_id="default",kind="properties",executor="server")
    store.submit(spec)
    remote=RemoteExecutor(ServerProfile(name="server",host="example.org"))
    row={"state":"completed","message":"任务已完成","updated":"frozen-time",
         "result":json.dumps({"output":"/remote/result"}),"log_tail":"中文日志"}
    calls=[]
    def broken(action,*args,**kwargs):
        calls.append(action)
        if action=="status":return json.dumps(row).encode()
        raise TimeoutError("connection interrupted")
    monkeypatch.setattr(remote,"command",broken)
    remote.sync(spec.task_id,store)
    assert store.task(spec.task_id)["state"]=="completed"
    assert "download_error" in json.loads(store.task(spec.task_id)["result"])
    buf=io.BytesIO()
    with tarfile.open(fileobj=buf,mode="w:gz") as tar:
        data=b"{}";entry=tarfile.TarInfo(".aqua-download.json");entry.size=len(data);tar.addfile(entry,io.BytesIO(data))
    def recovered(action,*args,**kwargs):
        calls.append(action)
        return json.dumps(row).encode() if action=="status" else buf.getvalue()
    monkeypatch.setattr(remote,"command",recovered)
    remote.sync(spec.task_id,store)
    remote.sync(spec.task_id,store,download=False)
    remote.sync(spec.task_id,store)
    assert calls.count("export")==2
    assert json.loads(store.task(spec.task_id)["result"])["local_output"]
    assert (store.path.parent/"tasks"/spec.task_id/"process.log").read_text(encoding="utf-8")=="中文日志"
from pathlib import Path
import pandas as pd
import pytest
from bur01_structmine.unified.contracts import TaskSpec
from bur01_structmine.unified.execution import task_bundle,safe_extract
from bur01_structmine.unified.transfer import resolve_reference
from bur01_structmine.unified.datasets import file_hash
from bur01_structmine.unified.worker import restore_paths
from bur01_structmine.unified.export import export_task
from bur01_structmine.unified.store import Store
from bur01_structmine.unified.artifacts import read_manifest,resolve_table_paths


def test_long_prediction_filenames_survive_download(tmp_path):
    import os
    name='predictions/'+('candidate_'+'a'*80)+'/'+('result_'+'b'*150)+'.json'
    blob=io.BytesIO()
    with tarfile.open(fileobj=blob,mode='w:gz') as archive:
        info=tarfile.TarInfo(name);info.size=2;archive.addfile(info,io.BytesIO(b'{}'))
    destination=tmp_path/'download'
    safe_extract(blob.getvalue(),destination)
    target=destination/Path(name)
    if os.name=='nt':target=Path('\\\\?\\'+str(target.resolve()))
    assert target.read_bytes()==b'{}'

def test_nested_prepared_dependencies_preserve_hashes(tmp_path):
    source=tmp_path/"source";source.mkdir()
    coords=source/"protein.pdb";coords.write_text("coordinate fixture")
    prepared=source/"prepared.json"
    prepared.write_text(json.dumps({"protein_pdb":{"path":str(coords),"sha256":file_hash(coords)}}))
    batch=source/"batch.json"
    batch.write_text(json.dumps({"prepared":{"path":str(prepared),"sha256":file_hash(prepared)}}))
    expected=file_hash(prepared)
    spec=TaskSpec(task_id="transfer",project_id="default",kind="mechanism",inputs={"batch_manifest":str(batch)})
    blob=task_bundle(spec)
    coords.unlink();prepared.unlink();batch.unlink()
    destination=tmp_path/"copied";safe_extract(blob,destination)
    raw=json.loads((destination/"task.json").read_text())
    restored=restore_paths(raw["inputs"],destination)
    transferred=Path(restored["batch_manifest"])
    a=json.loads(transferred.read_text())["prepared"]
    prepared=resolve_reference(a["path"],transferred.parent)
    assert file_hash(prepared)==expected
    b=json.loads(prepared.read_text())["protein_pdb"]
    assert file_hash(resolve_reference(b["path"],prepared.parent))==b["sha256"]

def test_remote_download_rebases_coordinates_and_keeps_original_manifest(tmp_path):
    root=tmp_path/"deployment";store=Store(root/"state/projects.sqlite")
    spec=TaskSpec(task_id="remote",project_id="default",kind="mine")
    store.submit(spec)
    source=root/"state/tasks/remote/output/runs/run1";source.mkdir(parents=True)
    coordinate=root/"inputs/candidate.pdb";coordinate.parent.mkdir()
    coordinate.write_text("coordinate fixture")
    table=source/"candidates.tsv"
    pd.DataFrame([{"structure_file":str(coordinate)}]).to_csv(table,sep="\t",index=False)
    manifest=source/"run_manifest.json"
    manifest.write_text(json.dumps({"run_id":"run1","query_path":str(coordinate),"artifacts":{"annotated_tsv":str(table)}}))
    original=manifest.read_bytes()
    buffer=io.BytesIO();export_task(store,"remote",buffer)
    dest=tmp_path/"client";safe_extract(buffer.getvalue(),dest)
    # Server paths deliberately become unavailable.
    coordinate.unlink();table.unlink()
    read=read_manifest(dest/"runs/run1/run_manifest.json","_run_dir")
    assert Path(read["query_path"]).read_text()=="coordinate fixture"
    frame=pd.read_csv(read["artifacts"]["annotated_tsv"],sep="\t")
    rebased=resolve_table_paths(frame,Path(read["artifacts"]["annotated_tsv"]).parent)
    assert Path(rebased.iloc[0].structure_file).is_file()
    assert manifest.read_bytes()==original

def test_sabio_mutants_require_explicit_construct_and_parameter_mapping(tmp_path):
    from bur01_structmine.unified.sabio import map_assays
    # Small API-shape fixture; no fixture is a released measurement dataset.
    entry={"id":123,"enzyme_description":{"wildtype":"Mutant","mutant_spec":"A2G"},
           "general":{"experiment_type":"InVitro"},
           "reaction":{"species":[{"role":"Substrate","compound":{"id":1,"name":"ethanol"}},
                                  {"role":"Substrate","compound":{"id":2,"name":"acid"}}]},
           "external_links":{"compound":[{"key":"SMILES","id":1,"value":"CCO"},{"key":"SMILES","id":2,"value":"CC(=O)O"}]},
           "experimental_conditions":{"envvar_ph":{"start_value":7},"envvar_temperature":{"start_value":25,"unit":"°C"}},
           "kineticlaw":{"parameter":[{"name":"cat","parameter_type":{"name":"kcat"},"unit":{"n_name":"s^(-1)"},
                                      "n_start_value":2,"species":{}},
                                     {"name":"km","parameter_type":{"name":"Km"},"unit":{"n_name":"M"},
                                      "n_start_value":.01,"species":{"species_key":"1 | ethanol | Substrate"}}]}}
    page=tmp_path/"page.json";page.write_text(json.dumps({"data":[entry]}))
    (tmp_path/"manifest.json").write_text(json.dumps({"pages":[{"path":"page.json","sha256":file_hash(page)}]}))
    mapping=tmp_path/"map.json"
    curated={"sequence":"AGDE","sequence_source":"fixture","construct_confirmed":True}
    mapping.write_text(json.dumps({"123":curated}))
    assert map_assays(tmp_path,mapping,tmp_path/"unconfirmed")["rows"]==0
    curated["mutations_confirmed"]=True;mapping.write_text(json.dumps({"123":curated}))
    result=map_assays(tmp_path,mapping,tmp_path/"confirmed")
    assert result["rows"]==1 and result["unmapped"]==1
    frame=pd.read_csv(tmp_path/"confirmed/mapped.csv")
    assert frame.iloc[0].temperature_K==298.15 and frame.iloc[0].endpoint=="Km"
