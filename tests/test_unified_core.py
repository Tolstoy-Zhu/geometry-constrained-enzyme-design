import io,json,tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pytest
from bur01_structmine.unified.contracts import EnzymeRecord,AssayRecord,TaskSpec
from bur01_structmine.unified.store import Store,uid,index_history
from bur01_structmine.unified.execution import worker_lock,task_bundle,safe_extract
from bur01_structmine.unified.services import cancel_task,import_records

def task(**kwargs):
    return TaskSpec(task_id=uid("test"),project_id="default",kind="index",inputs={},**kwargs)

def test_contract_identity_and_schema():
    row=EnzymeRecord(enzyme_id="e",sequence="ACD",name="core",construct_id="core_3aa",source="natural")
    with pytest.raises(ValueError):EnzymeRecord(**{**row.model_dump(),"sequence":"ACE"})
    with pytest.raises(ValueError):EnzymeRecord(**{**row.model_dump(),"schema_version":"2"})
    with pytest.raises(ValueError):EnzymeRecord(**{**row.model_dump(),"sequence":"ACX","sequence_sha256":""})

def test_censor_and_unmeasured_are_not_exact_zeros():
    raw=dict(assay_id="a",enzyme_id="e",reaction_id="r",endpoint="kcat",unit="s^-1",
             method="HPLC",replicate_id="1",measurement_group="g",source="paper")
    with pytest.raises(ValueError):AssayRecord(**raw,qualifier="<")
    with pytest.raises(ValueError):AssayRecord(**raw,qualifier="not_measured",value=0)
    assert AssayRecord(**raw,qualifier="<",detection_limit=.01).value is None

def test_atomic_claim_and_idempotent_submission(tmp_path):
    store=Store(tmp_path/"state.sqlite")
    specs=[task() for _ in range(75)]
    with ThreadPoolExecutor(max_workers=8) as executor:list(executor.map(store.submit,specs))
    assert len(store.tasks())==75
    store.submit(specs[0])
    with ThreadPoolExecutor(max_workers=8) as executor:claims=list(executor.map(lambda _:store.claim(),range(8)))
    assert sum(item is not None for item in claims)==1
    spec=specs[0].model_copy(deep=True);spec.inputs["different"]=True
    with pytest.raises(ValueError,match="different inputs"):store.submit(spec)

def test_remote_jobs_are_never_claimed_locally(tmp_path):
    store=Store(tmp_path/"state.sqlite");store.submit(task(executor="autodl"))
    assert store.claim() is None
    store.recover()
    assert store.tasks()[0]["state"]=="queued"


def test_late_poll_preserves_completed_download_and_newer_state(tmp_path):
    store=Store(tmp_path/'state.sqlite');spec=task(executor='autodl');store.submit(spec)
    out=tmp_path/'download';out.mkdir();(out/'.aqua-download.json').write_text('{}')
    stamp='2026-11-24T13:00:00Z'
    downloaded={'remote_updated':stamp,'remote':{'done':True},'local_output':str(out)}
    polling={'remote_updated':stamp,'remote':{'done':True}}
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda r:store.update(spec.task_id,'completed',result=r),[downloaded,polling]*20))
    assert json.loads(store.task(spec.task_id)['result'])['local_output']==str(out)
    store.update(spec.task_id,'running',result={'remote_updated':'2026-11-24T12:59:00Z'})
    assert store.task(spec.task_id)['state']=='completed'
    store.update(spec.task_id,'running',result={'remote_updated':'2026-11-24T13:01:00Z'})
    assert 'local_output' not in json.loads(store.task(spec.task_id)['result'])

def test_recovery_does_not_duplicate_live_orphan(tmp_path):
    store=Store(tmp_path/"state.sqlite");spec=task();store.submit(spec);store.claim()
    with worker_lock(store.path.parent/"tasks"/spec.task_id):
        store.recover()
        assert store.task(spec.task_id)["state"]=="running"
    store.recover()
    assert store.task(spec.task_id)["state"]=="interrupted"

def test_cancellation_project_scope(tmp_path):
    store=Store(tmp_path/"state.sqlite");spec=task();store.submit(spec)
    with pytest.raises(ValueError,match="different project"):cancel_task(spec.task_id,"another",store)
    assert not store.task(spec.task_id)["cancel_requested"]
    cancel_task(spec.task_id,"default",store)
    assert store.task(spec.task_id)["state"]=="cancelled"

def test_design_history_is_project_scoped(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from bur01_structmine import design_view
    legacy=tmp_path/'legacy';(legacy/'old').mkdir(parents=True)
    (legacy/'old/manifest.json').write_text('{}')
    imported=tmp_path/'imported';imported.mkdir();(imported/'manifest.json').write_text('{}')
    monkeypatch.setattr(design_view,'design_root',lambda:legacy)
    monkeypatch.setattr(design_view,'st',SimpleNamespace(session_state={'aqua_project':'current'}))
    store=Store();store.record('history','imported','current',{'area':'designs','path':str(imported/'manifest.json')})
    assert design_view.design_runs()=={'imported':imported}
    design_view.st.session_state['aqua_project']='another'
    assert design_view.design_runs()=={}
    design_view.st.session_state['aqua_project']='default'
    assert design_view.design_runs()=={'old':legacy/'old'}

@pytest.mark.parametrize("kind",["properties","mechanism","features"])
def test_concurrent_retry_is_single_job(tmp_path,kind):
    store=Store(tmp_path/"state.sqlite")
    spec=TaskSpec(task_id=uid("test"),project_id="default",kind=kind,inputs={})
    store.submit(spec);store.update(spec.task_id,"interrupted")
    with ThreadPoolExecutor(max_workers=8) as pool:
        retried=list(pool.map(lambda _:store.retry(spec.task_id),range(20)))
    assert len(set(retried))==1
    assert sum(t["state"]=="queued" for t in store.tasks())==1
    assert (retried[0]==spec.task_id)==(kind in {"mechanism","features"})

@pytest.mark.parametrize("name",["../escape","/absolute","C:/escape","bad\\file"])
def test_archive_path_rejection(tmp_path,name):
    blob=io.BytesIO()
    with tarfile.open(fileobj=blob,mode="w:gz") as archive:
        info=tarfile.TarInfo(name);info.size=1;archive.addfile(info,io.BytesIO(b"x"))
    with pytest.raises(ValueError):safe_extract(blob.getvalue(),tmp_path/"out")

def test_staged_inputs_survive_original_removal(tmp_path):
    from bur01_structmine.unified.worker import restore_paths
    source=tmp_path/"input.fasta";source.write_text(">p\nACD")
    spec=task();spec.inputs={"fasta_path":str(source),"fasta_text":">q\nACDE"}
    blob=task_bundle(spec);source.unlink()
    destination=tmp_path/"copied";safe_extract(blob,destination)
    data=json.loads((destination/"task.json").read_text())
    restored=restore_paths(data["inputs"],destination)
    assert Path(restored["fasta_path"]).read_text()==">p\nACD"
    assert restored["fasta_text"]==">q\nACDE"
    safe_extract(blob,destination)
    Path(restored["fasta_path"]).write_text("changed")
    with pytest.raises(ValueError,match="differs"):safe_extract(blob,destination)

def test_history_indexes_snapshots_without_mutating_files(tmp_path):
    root=tmp_path/"repo";path=root/"results/nr_runs/old/nr_run_manifest.json"
    path.parent.mkdir(parents=True);path.write_text('{"status":"running"}')
    store=Store(tmp_path/"state.sqlite")
    assert index_history(root,store)==1
    path.write_text('{"status":"completed"}')
    assert index_history(root,store)==1
    assert len(store.records("history"))==2
    assert path.read_text()=='{"status":"completed"}'

def test_experimental_import_rolls_back_conflicts(tmp_path):
    store=Store(tmp_path/"state.sqlite")
    enzyme=dict(enzyme_id="e",sequence="ACD",name="e",source="natural",construct_id="core")
    path=tmp_path/"input.json";path.write_text(json.dumps({"enzyme":[enzyme]}))
    import_records(path,store=store)
    path.write_text(json.dumps({"enzyme":[{**enzyme,"enzyme_id":"new"},{**enzyme,"name":"conflict"}]}))
    with pytest.raises(ValueError):import_records(path,store=store)
    assert len(store.records("enzyme"))==1
