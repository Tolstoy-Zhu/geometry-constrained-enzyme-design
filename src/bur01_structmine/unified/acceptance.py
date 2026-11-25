"""Point-in-time platform acceptance with separate software, model and research states."""
from __future__ import annotations
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from .contracts import EnzymeRecord,ReactionSpec,AssayRecord
from .store import Store,now,uid
from .datasets import file_hash
from .mechanism import atomic_json


def feedback_roundtrip(output,learning_python=None):
    """Exercise import/export in a separate index; never invent a research measurement."""
    from .services import import_records,export_assays
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    enzyme={'enzyme_id':'acceptance_example','name':'Software fixture','sequence':'ACDEFGHIK',
            'source':'designed','construct_id':'fixture','provenance':{'data_scope':'software_fixture'}}
    reaction={'reaction_id':'fixture_reaction','name':'Software fixture','molecules':[{'name':'ethanol','smiles':'CCO'}]}
    assay={'assay_id':'unmeasured','enzyme_id':enzyme['enzyme_id'],'reaction_id':reaction['reaction_id'],
        'endpoint':'kcat','value':None,'unit':'s^-1','method':'software_fixture','replicate_id':'r1',
        'measurement_group':'fixture','qualifier':'not_measured','source':'isolated software acceptance; no experimental observation'}
    simulated={**assay,'assay_id':'simulated_control','value':1.,'qualifier':'=','data_kind':'simulated'}
    payload={'enzyme':[enzyme],'reaction':[reaction],'assay':[assay,simulated]}
    source=output/'fixture_records.json';atomic_json(source,payload)
    with tempfile.TemporaryDirectory(prefix='aqua-acceptance-') as temp:
        store=Store(Path(temp)/'projects.sqlite')
        imported=import_records(source,store=store)
        assert import_records(source,store=store)==imported
        assert len(store.records('assay'))==2
        exported=output/'fixture_assays.csv';export_assays(exported,store=store)
        conflict={**payload,'enzyme':[{**enzyme,'name':'Conflicting name'}]}
        atomic_json(output/'conflict.json',conflict)
        try:import_records(output/'conflict.json',store=store)
        except ValueError:pass
        else:raise ValueError('Conflicting scientific record was overwritten')
        if store.records('enzyme')[0]['name']!=enzyme['name']:raise ValueError('Conflict changed original record')
    interpreter=learning_python or os.environ.get('AQUA_LEARNING_PYTHON',sys.executable)
    receipt=output/'training_filter.json'
    subprocess.run([interpreter,'-m','bur01_structmine.unified.acceptance','--filter',str(exported),
                    '--output',str(receipt)],check=True,capture_output=True,timeout=60)
    filtered=json.loads(receipt.read_text(encoding='utf-8'))
    if filtered['accepted']!=0 or set(filtered['reasons'])!={'not_experimental','censored_or_unmeasured'}:
        raise ValueError('Simulation or unmeasured value reached experimental training input')
    return {'scope':'isolated_software_fixture_never_research_dataset','imported':imported,
            'idempotent_import':True,'conflict_preserved':True,'training_filter':filtered}


def run_acceptance(project_id,output,store=None,learning_python=None,ui_evidence=None,recovery_evidence=None):
    from .families import profiles,reference
    from .services import property_reports,kinetic_reports
    from .mechanism_results import load_result
    from .assessment_results import assessment_tasks,load_assessment
    store=store or Store();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if (output/'platform_acceptance.json').exists():raise ValueError('Use a new output directory for a new acceptance snapshot')
    project=next((r for r in store.projects() if r['id']==project_id),None)
    if not project:raise ValueError('Unknown project')
    checks=[];evidence={};research={'mechanism_rates':'not_established','full_NR_acceptance':'pending_separate_validation'}
    def keep(path):
        path=Path(path).resolve();evidence[str(path)]=file_hash(path);return path
    def check(key,title,operation):
        try:detail=operation();checks.append({'id':key,'title':title,'status':'passed','detail':detail})
        except (OSError,ValueError,KeyError,TypeError,AssertionError,subprocess.SubprocessError) as exc:
            checks.append({'id':key,'title':title,'status':'requires_attention','detail':str(exc)})
    def families():
        rows=[]
        for profile in profiles():
            keep(reference(profile))
            rows.append({'family_id':profile.family_id,'catalytic_roles':len(profile.sites),'wase_screen':profile.wase_screen})
        if len(rows)<3:raise ValueError('Expected three preview family configurations')
        return rows
    check('families','家族配置与参考结构文件',families)
    enzymes=store.records('enzyme',project_id);reactions=store.records('reaction',project_id);assays=store.records('assay',project_id)
    def records():
        for row in enzymes:EnzymeRecord.model_validate(row)
        for row in reactions:ReactionSpec.model_validate(row)
        for row in assays:
            record=AssayRecord.model_validate(row)
            if record.enzyme_id not in {e['enzyme_id'] for e in enzymes} or record.reaction_id not in {r['reaction_id'] for r in reactions}:
                raise ValueError('Assay references another or unknown project record')
        if not enzymes or not reactions:raise ValueError('Enzyme and reaction records are required')
        return {'enzymes':len(enzymes),'reactions':len(reactions),'assays':len(assays)}
    check('records','序列身份、构建和反应记录',records)
    research['experimental_assays']=sum(r['data_kind']=='experimental' and r['qualifier']=='=' for r in assays)
    def structures():
        rows=store.records('structure_batch_review',project_id)
        if not rows:raise ValueError('No complete-construct review report')
        row=rows[-1];path=keep(row['report_path'])
        if file_hash(path)!=row['sha256']:raise ValueError('Structure-review source changed')
        report=json.loads(path.read_text(encoding='utf-8'))
        if not report['rows'] or not all(r['structure_review_pass'] for r in report['rows']):raise ValueError('Some complete constructs require review')
        return {'reviewed':len(report['rows']),'passed':sum(r['structure_review_pass'] for r in report['rows'])}
    check('structures','完整构建结构复核报告',structures)
    def available_reports(getter,filename):
        reports=getter(store,project_id)
        if not reports:raise ValueError('No readable completed prediction report')
        for path in reports.values():keep(path/filename)
        return {'readable_reports':len(reports),'scope':'result availability; predictive validity evaluated separately'}
    check('properties','性质预测报告读取',lambda:available_reports(property_reports,'index.html'))
    check('kinetics','自主动力学模型结果读取',lambda:available_reports(kinetic_reports,'predictions.json'))
    def models():
        rows=store.records('model',project_id)
        if not rows:raise ValueError('No registered trained models')
        for row in rows:
            card=row['card']
            if not all(k in card for k in ['model_id','endpoint','model_sha256','metrics','status']):raise ValueError('Incomplete model card')
        return [{'endpoint':r['card']['endpoint'],'status':r['card']['status']} for r in rows]
    check('models','模型版本与分项评估记录',models)
    tasks=store.tasks(project_id)
    production=[r for r in tasks if json.loads(r['spec'])['kind']=='mechanism' and json.loads(r['spec'])['inputs'].get('action')=='batch' and r['state']=='completed']
    paths={}
    def mechanism():
        if not production:raise ValueError('No completed production batch')
        reviews={r['task_id']:r for r in store.records('mechanism_quality',project_id)};rows=[]
        for task in production:
            result=json.loads(task['result']);meta=result.get('remote') or result
            root=Path(result['local_output'])/'mechanism' if result.get('local_output') else Path(meta['output'])
            path=keep(root/'batch_report.json');paths[task['id']]=path
            if task['id'] not in reviews:raise ValueError('Production quality report is missing')
            keep(reviews[task['id']]['report_path'])
            data=load_result(path,reviews[task['id']])
            if data['analyzed_replicas']!=data['replicas'] or data['quality_passed']!=data['replicas']:raise ValueError('Production analysis or numerical quality incomplete')
            rows.append({k:data[k] for k in ['replicas','production_ns','analyzed_replicas','quality_passed']})
        return rows
    check('production','生产轨迹、逐位点分析和数值质检',mechanism)
    def assessment():
        if not paths:raise ValueError('Production result not verified')
        summaries=[]
        for task_id,path in paths.items():
            rows=assessment_tasks(project_id,task_id,store)
            done=rows[0] if rows and rows[0]['state']=='completed' else None
            if not done:raise ValueError('Parameter and sampling assessment is not complete')
            report,folder=load_assessment(done,task_id,path);keep(folder/'assessment_report.json')
            summaries.append({'replicas':report['replicas'],'parameter_evaluations':report['channel']['frame_parameter_evaluations'],
                              'sampling_groups':len(report['sampling']['groups'])})
        return summaries
    check('mechanism_assessment','通道参数与采样评估执行及报告绑定',assessment)
    def history():
        rows=store.records('history',project_id)
        if not rows:raise ValueError('No indexed history available')
        for row in rows:
            if file_hash(keep(row['path']))!=row['manifest_sha256']:raise ValueError('Historical manifest changed: '+row['run_id'])
        return {'unchanged_history_manifests':len(rows)}
    check('history','历史结果原始文件保持不变',history)
    check('feedback','实验回填接口与训练数据排除规则',lambda:feedback_roundtrip(output/'feedback_fixture',learning_python))
    def ui():
        if not ui_evidence:raise ValueError('Web acceptance evidence has not been supplied')
        receipt=json.loads(keep(ui_evidence).read_text(encoding='utf-8'))
        if receipt.get('status')!='passed' or len(receipt.get('tabs',[]))<13:raise ValueError('Web pages have not all been verified')
        for path,sha in receipt.get('source_sha256',{}).items():
            if file_hash(keep(path))!=sha:raise ValueError('UI changed after validation')
        if not receipt.get('source_sha256'):raise ValueError('UI source identity missing')
        return {'pages':len(receipt['tabs']),'checks':receipt.get('checks',[])}
    check('web','网页逐页、中文、筛选和下载验收',ui)
    def recovery():
        if not recovery_evidence:raise ValueError('Remote recovery evidence has not been supplied')
        receipt=json.loads(keep(recovery_evidence).read_text(encoding='utf-8'))
        row=store.task(receipt['task_id'])
        if row['project_id']!=project_id or row['state']!='completed':raise ValueError('Recovered task is not complete in this project')
        if receipt.get('status')!='passed' or receipt.get('unchanged_saved_frames',0)<1 or not receipt.get('same_task_id'):
            raise ValueError('Checkpoint preservation was not established')
        return {'unchanged_saved_frames':receipt['unchanged_saved_frames'],'same_task_id':True}
    check('remote_recovery','远端分析中断恢复与已有结果保留',recovery)
    report={'schema_version':'1.0.0','project_id':project_id,'created_at':now(),
        'status':'engineering_preview_checks_passed' if all(r['status']=='passed' for r in checks) else 'requires_attention',
        'checks':checks,'evidence_sha256':evidence,'research_validation':research,
        'task_states':dict(Counter(r['state'] for r in tasks)),
        'scope':'Point-in-time engineering acceptance; model accuracy, full NR and experimental mechanism validation remain separate.'}
    path=output/'platform_acceptance.json';atomic_json(path,report)
    store.record('platform_acceptance',uid('acceptance'),project_id,{'report_path':str(path.resolve()),'sha256':file_hash(path)})
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--filter',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args()
    import pandas as pd
    from .datasets import canonicalize
    accepted,rejected=canonicalize(pd.read_csv(args.filter))
    atomic_json(args.output,{'accepted':len(accepted),'rejected':len(rejected),'reasons':rejected.reason.tolist()})
