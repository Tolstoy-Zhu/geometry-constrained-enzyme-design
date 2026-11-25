import json
import io
import tarfile
from pathlib import Path
import pytest
from bur01_structmine.unified.store import Store
from bur01_structmine.unified.contracts import TaskSpec
from bur01_structmine.unified.datasets import file_hash
from bur01_structmine.unified.assessment_results import load_assessment,assessment_tasks,submit_assessment
from bur01_structmine.unified.acceptance import feedback_roundtrip,run_acceptance
from bur01_structmine.unified.mechanism_assessment import AssessmentProtocol


def test_feedback_uses_isolated_records_and_excludes_unmeasured_and_simulated(tmp_path):
    pytest.importorskip('rdkit')
    store=Store();before=store.records('assay')
    result=feedback_roundtrip(tmp_path/'feedback')
    assert result['training_filter']['accepted']==0
    assert result['training_filter']['rejected']==2
    assert result['conflict_preserved'] and result['idempotent_import']
    assert store.records('assay')==before


def assessment_files(tmp_path):
    source=tmp_path/'production';source.mkdir()
    (source/'batch_report.json').write_text('{"production_ns":120}')
    (source/'batch_identity.json').write_text('{"manifest_sha256":"batch"}')
    out=tmp_path/'assessment';out.mkdir();(out/'assessment_identity.json').write_text('{}')
    report={'production_task_id':'production','identity':{'batch_sha256':'batch',
        'batch_report_sha256':file_hash(source/'batch_report.json'),
        'analysis_identity_sha256':file_hash(out/'assessment_identity.json')}}
    path=out/'assessment_report.json';path.write_text(json.dumps(report))
    task={'result':json.dumps({'output':str(out),'report_sha256':file_hash(path)})}
    return task,source,out


def test_assessment_is_bound_to_actual_original_report(tmp_path):
    task,source,out=assessment_files(tmp_path)
    assert load_assessment(task,'production',source/'batch_report.json')[1]==out
    (source/'batch_report.json').write_text('{"production_ns":240}')
    with pytest.raises(ValueError,match='another production'):load_assessment(task,'production',source/'batch_report.json')


def test_changed_analysis_policy_or_task_cannot_reuse_pass_status(tmp_path):
    task,source,out=assessment_files(tmp_path)
    with pytest.raises(ValueError,match='another task'):load_assessment(task,'another',source/'batch_report.json')
    (out/'assessment_identity.json').write_text('{"changed":true}')
    with pytest.raises(ValueError,match='policy'):load_assessment(task,'production',source/'batch_report.json')


def test_assessment_submission_and_visibility_are_project_scoped(tmp_path):
    store=Store(tmp_path/'state.sqlite')
    batch=TaskSpec(task_id='production',project_id='default',kind='mechanism',inputs={'action':'batch'})
    store.submit(batch);store.update(batch.task_id,'completed')
    analysis=TaskSpec(task_id='analysis',project_id='default',kind='mechanism',resource='cpu',
                      inputs={'action':'assess','production_task_id':'production',
                              'protocol':AssessmentProtocol().model_dump()})
    store.submit(analysis)
    assert submit_assessment('default','production',store)=='analysis'
    assert assessment_tasks('another','production',store)==[]
    with pytest.raises(ValueError):submit_assessment('another','production',store)


def test_empty_project_does_not_pass_acceptance(tmp_path,monkeypatch):
    from bur01_structmine.unified import acceptance
    monkeypatch.setattr(acceptance,'feedback_roundtrip',lambda *args,**kwargs:{'scope':'test_only'})
    store=Store(tmp_path/'state.sqlite');report=run_acceptance('default',tmp_path/'report',store)
    assert report['status']=='requires_attention'
    assert next(r for r in report['checks'] if r['id']=='production')['status']=='requires_attention'
    assert report['research_validation']['experimental_assays']==0
    assert report['research_validation']['full_NR_acceptance']=='pending_separate_validation'
    assert not store.records('assay')
    with pytest.raises(ValueError,match='new output'):run_acceptance('default',tmp_path/'report',store)


def test_local_assessment_references_existing_trajectories_without_copying(tmp_path):
    from bur01_structmine.unified.execution import task_bundle
    source=tmp_path/'production';source.mkdir();(source/'original.dcd').write_bytes(b'trajectory fixture')
    manifest=tmp_path/'batch.json';manifest.write_text('{}')
    quality=tmp_path/'quality.json';quality.write_text('{}')
    spec=TaskSpec(task_id='assess',project_id='default',kind='mechanism',inputs={
        'action':'assess','run':str(source),'batch_manifest':str(manifest),'quality_report':str(quality)})
    with tarfile.open(fileobj=io.BytesIO(task_bundle(spec)),mode='r:gz') as archive:
        assert archive.getnames()==['task.json']
        data=json.load(archive.extractfile('task.json'))
    assert data['inputs']['run']==str(source.resolve())
    assert (source/'original.dcd').read_bytes()==b'trajectory fixture'
