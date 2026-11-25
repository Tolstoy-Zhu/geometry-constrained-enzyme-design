"""Bind analysis reports to the selected project, production task and original batch."""
from __future__ import annotations
import json
from pathlib import Path
from .datasets import file_hash
from .store import Store


def assessment_tasks(project_id,production_task_id,store=None):
    store=store or Store();rows=[]
    for row in store.tasks(project_id):
        spec=json.loads(row['spec'])
        if spec['kind']=='mechanism' and spec['inputs'].get('action')=='assess' and spec['inputs'].get('production_task_id')==production_task_id:
            rows.append(row)
    return rows


def load_assessment(task,production_task_id,production_report):
    result=json.loads(task.get('result') or '{}');metadata=result.get('remote') or result
    root=Path(result['local_output'])/'mechanism' if result.get('local_output') else Path(metadata['output'])
    path=root/'assessment_report.json'
    if file_hash(path)!=metadata.get('report_sha256'):raise ValueError('Assessment report checksum mismatch')
    report=json.loads(path.read_text(encoding='utf-8'))
    if report['production_task_id']!=production_task_id:raise ValueError('Assessment belongs to another task')
    production_report=Path(production_report)
    identity=json.loads((production_report.parent/'batch_identity.json').read_text(encoding='utf-8'))
    if report['identity']['batch_sha256']!=identity['manifest_sha256'] or report['identity']['batch_report_sha256']!=file_hash(production_report):
        raise ValueError('Assessment belongs to another production result')
    if file_hash(root/'assessment_identity.json')!=report['identity']['analysis_identity_sha256']:
        raise ValueError('Assessment policy or source inventory changed')
    return report,root


def submit_assessment(project_id,production_task_id,store=None):
    from .execution import submit
    from .mechanism_assessment import AssessmentProtocol
    store=store or Store();row=store.task(production_task_id)
    spec=json.loads(row['spec'])
    if row['project_id']!=project_id or row['state']!='completed' or spec['kind']!='mechanism' or spec['inputs'].get('action')!='batch':
        raise ValueError('Select a completed production batch in this project')
    policy=AssessmentProtocol().model_dump()
    existing=[r for r in assessment_tasks(project_id,production_task_id,store)
              if json.loads(r['spec'])['inputs'].get('protocol')==policy]
    live=next((r for r in existing if r['state'] in {'queued','running','completed','blocked'}),None)
    if live:return live['id']
    reviews=[r for r in store.records('mechanism_quality',project_id) if r['task_id']==production_task_id]
    if not reviews:raise ValueError('A reviewed production batch is required')
    review=reviews[-1]
    if file_hash(review['report_path'])!=review['sha256']:raise ValueError('Quality report checksum mismatch')
    result=json.loads(row['result']);metadata=result.get('remote') or result
    run=('remote:' if spec['executor']!='local' else '')+metadata['output']
    return submit('mechanism',{'action':'assess','batch_manifest':spec['inputs']['batch_manifest'],
        'run':run,'quality_report':review['report_path'],'production_task_id':production_task_id,'workers':4,'protocol':policy},
        project_id=project_id,family_id=spec['family_id'],executor=spec['executor'],store=store)
