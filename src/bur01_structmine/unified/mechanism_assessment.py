"""Bounded, resumable analysis of completed trajectories; this module never runs MD."""
from __future__ import annotations
import hashlib
import json
import platform
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Literal
import numpy as np
import pandas as pd
from pydantic import Field, model_validator
from .contracts import Contract
from .datasets import file_hash
from .mechanism import atomic_json
from .mechanism_plan import MechanismBatch
from .sampling import assess_occupancy


class AssessmentProtocol(Contract):
    channel_method: Literal['point_grid_v1','continuous_edge_v2'] = 'continuous_edge_v2'
    spacings_A: list[float] = Field(default_factory=lambda:[.6,.8,1.])
    probes_A: list[float] = Field(default_factory=lambda:[1.,1.2,1.4,1.8])
    origin_fractions: list[float] = Field(default_factory=lambda:[0.,.5])
    channel_stride: int = Field(default=100,ge=1)
    include_last_frame: bool = True
    max_voxels: int = Field(default=4_000_000,ge=1000,le=8_000_000)
    grid_disagreement_fraction: float = Field(default=.1,ge=0,le=1)
    grid_bottleneck_range_A: float = Field(default=.3,gt=0)
    lags_frames: list[int] = Field(default_factory=lambda:[1,2,5,10,20,50])
    base_lag_frames: int = Field(default=5,ge=1)
    plateau_lags_frames: list[int] = Field(default_factory=lambda:[5,10,20])
    ck_multipliers: list[int] = Field(default_factory=lambda:[1,2,4,10])
    blocks: int = Field(default=5,ge=2,le=20)
    minimum_replicas: int = Field(default=3,ge=3)
    minimum_bidirectional_events: int = Field(default=5,ge=1)
    minimum_state_frames: int = Field(default=20,ge=2)
    minimum_timescale_lag_ratio: float = Field(default=2.,ge=1.)
    its_relative_range: float = Field(default=.2,gt=0,le=1)
    ck_absolute_error: float = Field(default=.1,gt=0,le=1)

    @model_validator(mode='after')
    def check(self):
        for key,low,high,count in [('spacings_A',.3,2.,2),('probes_A',.1,3.,2),('origin_fractions',0.,.999,2)]:
            values=getattr(self,key)
            if len(values)<count or len(set(values))!=len(values) or any(not low<=v<=high for v in values):
                raise ValueError('Invalid parameter sweep: '+key)
        for key in ['lags_frames','plateau_lags_frames','ck_multipliers']:
            values=getattr(self,key)
            if len(values)<2 or values!=sorted(set(values)) or any(not isinstance(v,int) or v<1 for v in values):
                raise ValueError('Ordered distinct positive lags required')
        if self.base_lag_frames not in self.lags_frames or not set(self.plateau_lags_frames)<=set(self.lags_frames):
            raise ValueError('Base and plateau lags must be evaluated')
        return self


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def source_inventory(batch_manifest,run,quality_report):
    batch_manifest=Path(batch_manifest);run=Path(run);quality_report=Path(quality_report)
    batch=MechanismBatch.model_validate_json(batch_manifest.read_text(encoding='utf-8'))
    identity=json.loads((run/'batch_identity.json').read_text(encoding='utf-8'))
    if file_hash(batch_manifest)!=identity['manifest_sha256']:raise ValueError('Wrong production batch')
    report=json.loads((run/'batch_report.json').read_text(encoding='utf-8'))
    quality=json.loads(quality_report.read_text(encoding='utf-8'))
    if quality['identity']['batch_sha256']!=identity['manifest_sha256']:raise ValueError('Wrong quality batch')
    sources={};items=[];assumptions=[]
    def checked(path,expected=None):
        path=Path(path).resolve()
        if path.is_symlink() or not path.is_file():raise ValueError('Missing source: '+str(path))
        sha=sources.setdefault(str(path),file_hash(path))
        if expected and sha!=expected:raise ValueError('Source checksum mismatch: '+str(path))
        return str(path)
    checked(batch_manifest);checked(run/'batch_identity.json');checked(run/'batch_report.json');checked(quality_report)
    qrows={r['replica_id']:r for r in quality['reports']}
    records={r['replica_id']:r for r in report['reports']}
    expected_ids={f"{c['candidate_id']}_r{i+1:02d}" for c in batch.candidates for i in range(batch.replicas)}
    if len(records)!=len(report['reports']) or set(records)!=expected_ids or set(qrows)!=expected_ids:
        raise ValueError('Incomplete or duplicate production/quality records')
    for candidate in batch.candidates:
        for i in range(batch.replicas):
            replica=f"{candidate['candidate_id']}_r{i+1:02d}";folder=run/replica
            row=records[replica]
            if row['candidate_id']!=candidate['candidate_id'] or row['production_ns']!=batch.production_ns:
                raise ValueError('Production record differs from frozen batch')
            checked(folder/'trajectory_report.json',qrows[replica]['identity']['trajectory_report_sha256'])
            checked(folder/'preflight.json')
            trajectory=json.loads((folder/'trajectory_report.json').read_text(encoding='utf-8'))
            check=json.loads((folder/'preflight.json').read_text(encoding='utf-8'))
            if check['manifest_sha256']!=candidate['prepared']['sha256']:raise ValueError('Prepared construct mismatch')
            topology=checked(check['files']['prmtop'],check['contract']['prmtop']['sha256'])
            trajectories=[checked(p['path'],p['sha256']) for p in trajectory['segments']]
            regions=candidate['analysis_regions']
            if not regions or set(regions)!=set(row['regions']):raise ValueError('Missing expected active-site analysis')
            if any(not key.isidentifier() for key in regions):raise ValueError('Invalid site identifier')
            for region in regions:
                checked(folder/('analysis_'+region)/'analysis.json')
                checked(folder/('analysis_'+region)/'occupancy.npz')
            assumptions.extend(check['contract'].get('protonation_assumptions',[]))
            items.append({'candidate_id':candidate['candidate_id'],'replica_id':replica,'topology':topology,
                'trajectories':trajectories,'regions':regions,'folder':str(folder),
                'numerical_quality_pass':qrows[replica].get('status') in {'passed','numerical_checks_passed_descriptors_available'} and not qrows[replica].get('issues')})
    return batch,items,sources,list(dict.fromkeys(assumptions))


def sampling_report(items,protocol):
    groups=defaultdict(list);dt_by_group={}
    for item in items:
        for region in item['regions']:
            root=Path(item['folder'])/('analysis_'+region)
            analysis=json.loads((root/'analysis.json').read_text(encoding='utf-8'))
            with np.load(root/'occupancy.npz',allow_pickle=False) as data:
                times=data['time_ps'];delta=np.diff(times)
                if len(times)<2 or not np.isfinite(times).all() or (delta<=0).any() or not np.allclose(delta,delta[0],rtol=1e-5):
                    raise ValueError('Uniform positive sampling interval required')
                if len(times)!=analysis['frames']:raise ValueError('Occupancy frame count mismatch')
                for solvent in item['regions'][region]['solvent_selections']:
                    occupied=data[solvent]
                    if occupied.ndim!=2 or occupied.shape[0]!=len(times) or occupied.dtype!=bool:
                        raise ValueError('Boolean molecule-by-frame occupancy required')
                    counts=occupied.sum(1)
                    if not np.isclose(counts.mean(),analysis['observables'][solvent]['mean_molecules_in_region']):
                        raise ValueError('Occupancy archive disagrees with original analysis')
                    key=(item['candidate_id'],region,solvent)
                    if key in dt_by_group and not np.isclose(dt_by_group[key],delta[0]):raise ValueError('Sampling intervals differ')
                    dt_by_group[key]=float(delta[0])
                    groups[key].append({'replica_id':item['replica_id'],'counts':counts})
    rows=[]
    for (candidate,region,solvent),replicas in groups.items():
        rows.append({'candidate_id':candidate,'site':region,'solvent':solvent,
            **assess_occupancy(replicas,dt_by_group[(candidate,region,solvent)],protocol)})
    return {'groups':rows,'state_scope':'site_occupancy_not_molecular_transport',
        'rates':None,'tpt_flux':None,'independence':'Chains assessed separately; replicates share one equilibrated starting structure per construct.'}


def assess_replica(item,protocol,output,identity_hash):
    import MDAnalysis as mda
    from MDAnalysis.transformations import unwrap
    from .mechanism_analysis import minimum_image
    from .channels import clearance_grid,grid_path,continuous_path_clearance
    universe=mda.Universe(item['topology'],item['trajectories'])
    protein=universe.select_atoms('protein')
    if not len(protein) or not len(protein.bonds):raise ValueError('Bonded protein topology required')
    universe.trajectory.add_transformations(unwrap(protein))
    sites={k:universe.select_atoms(v['site_selection']) for k,v in item['regions'].items()}
    if any(not len(s) for s in sites.values()):raise ValueError('Empty site selection')
    try:elements=list(protein.elements)
    except mda.exceptions.NoDataError:elements=[a.name.lstrip('0123456789')[0] for a in protein]
    indices=list(range(0,len(universe.trajectory),protocol['channel_stride']))
    if protocol['include_last_frame'] and indices[-1]!=len(universe.trajectory)-1:indices.append(len(universe.trajectory)-1)
    output=Path(output)/'frames'/item['replica_id'];output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for index in indices:
        path=output/f'{index:06d}.json'
        if path.is_file():
            saved=json.loads(path.read_text(encoding='utf-8'))
            if saved['identity_hash']!=identity_hash or saved['frame']!=index:raise ValueError('Changed analysis checkpoint')
            rows.extend(saved['rows']);continue
        ts=universe.trajectory[index]
        if ts.dimensions is None or not np.isfinite(ts.positions).all():raise ValueError('Invalid trajectory frame')
        fragments=protein.fragments;anchor=fragments[0].center_of_geometry()
        for fragment in fragments[1:]:
            offset=fragment.center_of_geometry()-anchor
            fragment.positions+=minimum_image(offset,ts.triclinic_dimensions)-offset
        frame_rows=[]
        for spacing in protocol['spacings_A']:
            for origin in protocol['origin_fractions']:
                grid=clearance_grid(protein.positions,elements,spacing_A=spacing,
                    max_voxels=protocol['max_voxels'],origin_fraction=origin)
                for region,site in sites.items():
                    for probe in protocol['probes_A']:
                        result=grid_path(grid,site.center_of_geometry(),probe,continuous_edges=protocol['channel_method']=='continuous_edge_v2')
                        points=result.pop('centerline_A',None)
                        continuous=continuous_path_clearance(points,grid['xyz'],grid['elements']) if points else None
                        valid=continuous is not None and continuous>=probe-1e-6
                        if points and protocol['channel_method']=='continuous_edge_v2' and not valid:
                            raise ValueError('Continuous-edge path failed independent clearance verification')
                        frame_rows.append({'candidate_id':item['candidate_id'],'replica_id':item['replica_id'],
                            'site':region,'frame':index,'time_ps':float(ts.time),'spacing_A':spacing,
                            'probe_radius_A':probe,'origin_fraction':origin,**result,
                            'continuous_bottleneck_A':continuous,'edge_clearance_pass':valid})
                del grid
        atomic_json(path,{'identity_hash':identity_hash,'frame':index,'rows':frame_rows})
        rows.extend(frame_rows)
    return rows


def summarize_channels(rows,protocol):
    groups=defaultdict(list)
    for row in rows:groups[(row['candidate_id'],row['replica_id'],row['site'],row['probe_radius_A'])].append(row)
    summary=[]
    expected={(s,o) for s in protocol['spacings_A'] for o in protocol['origin_fractions']}
    for (candidate,replica,site,probe),variants in groups.items():
        frames=defaultdict(list)
        for row in variants:frames[row['frame']].append(row)
        ranges=[];disagree=0;valid_total=0;edge_fail=0
        for values in frames.values():
            if len(values)!=len(expected) or {(v['spacing_A'],v['origin_fraction']) for v in values}!=expected:
                raise ValueError('Incomplete channel parameter sweep')
            valid=[v['edge_clearance_pass'] for v in values];valid_total+=sum(valid)
            disagree+=int(any(valid) and not all(valid))
            edge_fail+=sum(v['status']=='path_found' and not v['edge_clearance_pass'] for v in values)
            if all(valid):ranges.append(float(np.ptp([v['continuous_bottleneck_A'] for v in values])))
        fraction=disagree/len(frames)
        sensitive=fraction>protocol['grid_disagreement_fraction'] or bool(ranges and max(ranges)>protocol['grid_bottleneck_range_A'])
        status='grid_sensitive' if sensitive else ('no_paths_on_sampled_frames' if not valid_total else 'grid_agreement_on_sampled_frames')
        if edge_fail:status='edge_clearance_requires_review'
        summary.append({'candidate_id':candidate,'replica_id':replica,'site':site,'probe_radius_A':probe,
            'sampled_frames':len(frames),'settings_per_frame':len(expected),'grid_disagreement_fraction':fraction,
            'valid_path_fraction_across_settings':valid_total/len(variants),'matched_open_frames':len(ranges),
            'matched_bottleneck_range_max_A':max(ranges) if ranges else None,
            'edge_clearance_failures':edge_fail,'status':status})
    return summary


def assess_batch(batch_manifest,run,quality_report,output,protocol=None,workers=4,production_task_id=None):
    if not isinstance(workers,int) or not 1<=workers<=8:raise ValueError('Use 1 to 8 analysis workers')
    if isinstance(protocol,(str,Path)):protocol=json.loads(Path(protocol).read_text(encoding='utf-8'))
    protocol=AssessmentProtocol.model_validate(protocol or {}).model_dump()
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    batch,items,sources,assumptions=source_inventory(batch_manifest,run,quality_report)
    import MDAnalysis,scipy
    from . import channels,sampling
    identity={'sources':sources,'protocol':protocol,'code_sha256':{m.__name__:file_hash(m.__file__) for m in [channels,sampling]},
              'assessment_code_sha256':file_hash(__file__),'runtime':{'python':platform.python_version(),'numpy':np.__version__,
              'MDAnalysis':MDAnalysis.__version__,'scipy':scipy.__version__}}
    identity_hash=digest(identity);identity_path=output/'assessment_identity.json'
    if identity_path.exists() and json.loads(identity_path.read_text(encoding='utf-8'))!=identity:
        raise ValueError('Analysis inputs, policy, code or runtime changed; use a new task')
    atomic_json(identity_path,identity)
    sampling_data=sampling_report(items,protocol);atomic_json(output/'sampling_report.json',sampling_data)
    def progress(completed):
        atomic_json(output/'assessment_progress.json',{'status':'running','replicas_completed':completed,
                    'total':len(items),'identity_hash':identity_hash})
    rows=[];progress(0)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(assess_replica,item,protocol,str(output),identity_hash) for item in items]
        for completed,future in enumerate(as_completed(futures),1):
            rows.extend(future.result());progress(completed)
    rows.sort(key=lambda r:(r['candidate_id'],r['replica_id'],r['site'],r['frame'],r['spacing_A'],r['origin_fraction'],r['probe_radius_A']))
    summaries=summarize_channels(rows,protocol)
    pd.DataFrame(rows).to_csv(output/'channel_parameter_frames.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(summaries).to_csv(output/'channel_parameter_summary.csv',index=False,encoding='utf-8-sig')
    sample_rows=[{k:v for k,v in r.items() if k in {'candidate_id','site','solvent','independent_replicas','dt_ps','its_pass','ck_pass','occupancy_diagnostics_pass','status'}}|{'reasons':'; '.join(r['reasons'])} for r in sampling_data['groups']]
    pd.DataFrame(sample_rows).to_csv(output/'sampling_summary.csv',index=False,encoding='utf-8-sig')
    status_counts=dict(pd.Series([r['status'] for r in summaries]).value_counts().items())
    result={'schema_version':'1.0.0','status':'assessment_completed','production_task_id':production_task_id,
        'identity':{'batch_sha256':file_hash(batch_manifest),'batch_report_sha256':file_hash(Path(run)/'batch_report.json'),
                    'analysis_identity_sha256':file_hash(identity_path)},'protocol':protocol,
        'production_ns':len(items)*batch.production_ns,'replicas':len(items),'numerical_quality_passed':sum(r['numerical_quality_pass'] for r in items),
        'channel':{'rows':summaries,'frame_parameter_evaluations':len(rows),'status_counts':status_counts,
            'method':protocol['channel_method'],
            'scope':'Geometric probe sweep; probe radii are not calibrated substrate sizes. Grid agreement does not establish dynamical convergence.'},
        'sampling':sampling_data,'assumptions':assumptions,'rates':None,'tpt_flux':None,
        'scope':'Completed audit of a fixed pilot batch; no additional production simulation.'}
    atomic_json(output/'assessment_report.json',result)
    atomic_json(output/'assessment_progress.json',{'status':'completed','replicas_completed':len(items),'total':len(items),'identity_hash':identity_hash})
    return {'output':str(output),'status':'assessment_completed','replicas':len(items),
            'production_task_id':production_task_id,'report_sha256':file_hash(output/'assessment_report.json')}
