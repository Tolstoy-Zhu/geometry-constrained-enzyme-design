"""Display measured trajectory descriptors without treating chains as replicates."""
from __future__ import annotations
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from .datasets import file_hash


def finite(value, name, minimum=None, maximum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Invalid finite value: {name}")
    if minimum is not None and value < minimum or maximum is not None and value > maximum:
        raise ValueError(f"Out-of-range value: {name}")
    return float(value)


def summarize_batch(report, quality=None, expected_regions=None):
    """First average sites within each replica, then summarize independent replicas."""
    trajectories=report.get('reports',[])
    expected_regions=expected_regions or {}
    if not trajectories:raise ValueError('No production replicas in this report')
    ids=[r['replica_id'] for r in trajectories]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate replica identity')
    quality_rows={r['replica_id']:r for r in (quality or {}).get('reports',[])}
    if len(quality_rows)!=len((quality or {}).get('reports',[])):
        raise ValueError('Duplicate quality-review identity')
    if set(quality_rows)-set(ids):raise ValueError('Quality review contains another batch')
    site_rows=[];replica_rows=[];quality_table=[];complete=0;production=0.;durations_ns=set()
    for trajectory in trajectories:
        replica=trajectory['replica_id'];candidate=trajectory['candidate_id']
        duration=finite(trajectory['production_ns'],'production_ns',minimum=0.)
        if not duration:raise ValueError('Production duration must be positive')
        production+=duration
        durations_ns.add(duration)
        regions=trajectory.get('regions',{})
        expected=trajectory.get('expected_regions',expected_regions.get(replica))
        if expected is not None and (not isinstance(expected,list) or any(not isinstance(s,str) for s in expected) or len(expected)!=len(set(expected))):
            raise ValueError('Invalid expected analysis regions')
        if replica in expected_regions and trajectory.get('expected_regions') is not None and set(expected)!=set(expected_regions[replica]):
            raise ValueError('Conflicting expected analysis regions')
        analysis_complete=bool(expected) and set(regions)==set(expected)
        complete+=int(analysis_complete)
        q=quality_rows.get(replica)
        status='unreviewed'
        if q:
            status='passed' if q.get('status') in {'passed','numerical_checks_passed_descriptors_available'} and not q.get('issues') else 'requires_review'
        quality_row={'candidate_id':candidate,'replica_id':replica,'production_ns':duration,
            'analysis_complete':analysis_complete,'quality_status':status,'issues':'; '.join((q or {}).get('issues',[]))}
        for key,value in (q or {}).get('metrics',{}).items():
            if key=='core_rmsd_A':
                quality_row['core_rmsd_mean_A']=finite(value['mean'],key,minimum=0.)
                quality_row['core_rmsd_max_A']=finite(value['max'],key,minimum=0.)
        quality_table.append(quality_row)
        per_solvent=defaultdict(list)
        for region,analysis in regions.items():
            frames=analysis.get('frames',0)
            if isinstance(frames,bool) or not isinstance(frames,int) or frames<2:raise ValueError('Analysis needs at least two frames')
            radius=finite(analysis['radius_A'],'site radius',minimum=0.)
            paths=analysis.get('channel',{}).get('frame_results',[])
            widths=[finite(p['bottleneck_A'],'channel bottleneck',minimum=0.) for p in paths if p.get('status')=='path_found']
            for solvent,values in analysis.get('observables',{}).items():
                events=values.get('residence_events',[])
                durations=[finite(e['observed_duration_ps'],'observed residence',minimum=0.) for e in events]
                fully_observed=[d for d,e in zip(durations,events) if not e['left_censored'] and not e['right_censored']]
                row={'candidate_id':candidate,'replica_id':replica,'site':region,'solvent':solvent,
                    'frames':frames,'radius_A':radius,
                    'mean_molecules':finite(values['mean_molecules_in_region'],'mean occupancy',minimum=0.),
                    'occupied_fraction':finite(values['fraction_frames_occupied'],'occupied fraction',minimum=0.,maximum=1.),
                    'entry_events':finite(values['entry_events'],'entry count',minimum=0.),
                    'complete_residence_events':len(fully_observed),
                    'complete_event_observed_median_ps':median(fully_observed) if fully_observed else None,
                    'left_censored_events':sum(bool(e['left_censored']) for e in events),
                    'right_censored_events':sum(bool(e['right_censored']) for e in events),
                    'channel_sample_frames':len(paths),'path_found_fraction':len(widths)/len(paths) if paths else None,
                    'path_bottleneck_mean_A':mean(widths) if widths else None}
                site_rows.append(row);per_solvent[solvent].append(row)
        for solvent,rows in per_solvent.items():
            # Missing sites stay visible individually, but cannot enter a candidate average.
            if not analysis_complete or {r['site'] for r in rows}!=set(expected):continue
            if len({r['radius_A'] for r in rows})!=1:raise ValueError('Cannot average different site radii')
            replica_rows.append({'candidate_id':candidate,'replica_id':replica,'solvent':solvent,
                'sites':len(rows),'radius_A':rows[0]['radius_A'],
                **{key:mean(r[key] for r in rows) for key in ('mean_molecules','occupied_fraction','entry_events')}})
    if 'production_ns' in report and not math.isclose(production,finite(report['production_ns'],'reported total'),abs_tol=1e-8):
        raise ValueError('Production total disagrees with replica records')
    if report.get('replicas_completed',len(ids))!=len(ids):raise ValueError('Replica count disagrees with records')
    if len(durations_ns)!=1:raise ValueError('Different production durations require separate comparisons')
    for solvent in {r['solvent'] for r in replica_rows}:
        if len({r['radius_A'] for r in replica_rows if r['solvent']==solvent})!=1:
            raise ValueError('Cannot compare different site radii across candidates')
    groups=defaultdict(list)
    for row in replica_rows:groups[(row['candidate_id'],row['solvent'])].append(row)
    candidates=[]
    for (candidate,solvent),rows in groups.items():
        if len({r['radius_A'] for r in rows})!=1:raise ValueError('Cannot compare different site radii')
        entry={'candidate_id':candidate,'solvent':solvent,'replicates':len(rows),'radius_A':rows[0]['radius_A']}
        for key in ('mean_molecules','occupied_fraction','entry_events'):
            values=[r[key] for r in rows]
            entry.update({key+'_mean':mean(values),key+'_min':min(values),key+'_max':max(values)})
        candidates.append(entry)
    return {'production_ns':production,'replicas':len(ids),'analyzed_replicas':complete,
        'quality_passed':sum(r['quality_status']=='passed' for r in quality_table),
        'quality_available':len(quality_rows),'sampling_status':report.get('sampling_assessment','not_assessed'),
        'candidates':candidates,'replica_values':replica_rows,'sites':site_rows,'quality':quality_table}


def load_result(path, quality_record=None):
    """Check local report provenance, retaining the original downloaded manifest bytes."""
    path=Path(path);report=json.loads(path.read_text(encoding='utf-8'));expected={};assumptions=[]
    for row in report.get('reports',[]):
        replica=row['replica_id']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]{0,127}',replica):raise ValueError('Unsafe replica identifier')
        preflight=path.parent/replica/'preflight.json'
        if preflight.is_file():
            cfg=json.loads(preflight.read_text(encoding='utf-8')).get('contract',{})
            assumptions.extend(cfg.get('protonation_assumptions',[]))
    quality=None
    if quality_record:
        qpath=Path(quality_record['report_path'])
        if file_hash(qpath)!=quality_record['sha256']:raise ValueError('Quality report checksum mismatch')
        quality=json.loads(qpath.read_text(encoding='utf-8'))
        identity=json.loads((path.parent/'batch_identity.json').read_text(encoding='utf-8'))
        if quality.get('identity',{}).get('batch_sha256')!=identity['manifest_sha256']:
            raise ValueError('Quality report belongs to a different frozen batch')
        if quality_record.get('expected_regions'):
            if quality_record.get('batch_manifest_sha256')!=identity['manifest_sha256']:
                raise ValueError('Expected analysis regions belong to a different batch')
            expected=quality_record['expected_regions']
        replica_ids={r['replica_id'] for r in report.get('reports',[])}
        for row in quality.get('reports',[]):
            replica=row['replica_id']
            if replica not in replica_ids or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]{0,127}',replica):
                raise ValueError('Quality review has no matching prepared replica')
            original=path.parent/replica/'trajectory_report.json'
            if file_hash(original)!=row.get('identity',{}).get('trajectory_report_sha256'):
                raise ValueError('Quality review does not match this trajectory report')
    result=summarize_batch(report,quality,expected)
    result.update(assumptions=list(dict.fromkeys(assumptions)),quality_report_path=str(qpath) if quality_record else None)
    return result
