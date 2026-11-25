"""Replica-aware diagnostics of a two-state *occupancy* model, never transport rates."""
from __future__ import annotations
import numpy as np


def transition_counts(series,lag):
    """Count within each trajectory; never create transitions at replica boundaries."""
    if not isinstance(lag,int) or lag<1:raise ValueError('Positive integer lag required')
    counts=np.zeros((2,2),dtype=np.int64)
    for values in series:
        values=np.asarray(values)
        if values.ndim!=1 or not np.isin(values,[0,1]).all():raise ValueError('Binary state trajectory required')
        if len(values)<=lag:continue
        np.add.at(counts,(values[:-lag].astype(int),values[lag:].astype(int)),1)
    return counts


def transition_model(counts,lag_ps):
    counts=np.asarray(counts,float)
    if counts.shape!=(2,2) or not np.isfinite(counts).all() or (counts<0).any():
        raise ValueError('Invalid transition counts')
    rows=counts.sum(1)
    if (rows==0).any():return {'status':'unobserved_state','matrix':None,'timescale_ps':None}
    matrix=counts/rows[:,None];eigenvalue=float(np.trace(matrix)-1)
    status='estimated' if 0<eigenvalue<1 else ('disconnected' if eigenvalue>=1 else 'unresolved_or_oscillatory')
    return {'status':status,'matrix':matrix.tolist(),'eigenvalue':eigenvalue,
            'timescale_ps':float(-lag_ps/np.log(eigenvalue)) if status=='estimated' else None}


def assess_occupancy(replicas,dt_ps,protocol):
    """Preserve per-replica evidence and held-out CK errors; no pseudo-independent chains."""
    if not np.isfinite(dt_ps) or dt_ps<=0:raise ValueError('Positive sampling interval required')
    if len({r['replica_id'] for r in replicas})!=len(replicas):raise ValueError('Duplicate independent replica')
    series=[np.asarray(r['counts']) for r in replicas]
    if any(a.ndim!=1 or len(a)<2 or not np.isfinite(a).all() or (a<0).any() for a in series):
        raise ValueError('Finite molecule count series required')
    binary=[(a>0).astype(np.int8) for a in series]
    stats=[]
    for replica,counts,states in zip(replicas,series,binary):
        changes=transition_counts([states],1)
        blocks=[float(b.mean()) for b in np.array_split(counts,protocol['blocks']) if len(b)]
        stats.append({'replica_id':replica['replica_id'],'frames':len(states),
            'mean_molecules':float(counts.mean()),'occupied_fraction':float(states.mean()),
            'empty_frames':int((states==0).sum()),'occupied_frames':int(states.sum()),
            'empty_to_occupied':int(changes[0,1]),'occupied_to_empty':int(changes[1,0]),
            'block_mean_molecules':blocks,'first_half_mean':float(counts[:len(counts)//2].mean()),
            'second_half_mean':float(counts[len(counts)//2:].mean())})
    lags=protocol['lags_frames'];base=protocol['base_lag_frames']
    models=[]
    for lag in lags:
        counts=transition_counts(binary,lag)
        models.append({'lag_frames':lag,'lag_ps':lag*dt_ps,'counts':counts.tolist(),
                       **transition_model(counts,lag*dt_ps)})
    by_lag={r['lag_frames']:r for r in models}
    plateau=[by_lag[k]['timescale_ps'] for k in protocol['plateau_lags_frames']]
    resolved=all(t is not None and t>=protocol['minimum_timescale_lag_ratio']*lag*dt_ps
                 for t,lag in zip(plateau,protocol['plateau_lags_frames']))
    spread=float(np.ptp(plateau)/np.mean(plateau)) if resolved else None
    plateau_pass=resolved and spread<=protocol['its_relative_range']
    ck=[]
    for held,replica in enumerate(replicas):
        train=[s for i,s in enumerate(binary) if i!=held]
        training_counts=transition_counts(train,base)
        model=transition_model(training_counts,base*dt_ps)
        for multiple in protocol['ck_multipliers']:
            observed_counts=transition_counts([binary[held]],base*multiple)
            observed=transition_model(observed_counts,base*multiple*dt_ps)
            sufficient=(training_counts.sum(1)>=protocol['minimum_state_frames']).all() and (
                observed_counts.sum(1)>=protocol['minimum_state_frames']).all()
            eligible=sufficient and model['matrix'] is not None and observed['matrix'] is not None
            predicted=np.linalg.matrix_power(np.asarray(model['matrix']),multiple) if eligible else None
            error=float(np.max(np.abs(predicted-np.asarray(observed['matrix'])))) if eligible else None
            ck.append({'held_out_replica':replica['replica_id'],'multiple':multiple,
                'lag_ps':base*multiple*dt_ps,'status':'evaluated' if eligible else 'insufficient_state_counts',
                'max_absolute_error':error,'predicted':predicted.tolist() if eligible else None,
                'observed':observed['matrix'],'held_out_counts':observed_counts.tolist()})
    ck_pass=bool(ck) and all(r['max_absolute_error'] is not None and r['max_absolute_error']<=protocol['ck_absolute_error'] for r in ck)
    enough_replicas=len(replicas)>=protocol['minimum_replicas']
    transitions_pass=all(min(r['empty_to_occupied'],r['occupied_to_empty'])>=protocol['minimum_bidirectional_events'] for r in stats)
    states_pass=all(min(r['empty_frames'],r['occupied_frames'])>=protocol['minimum_state_frames'] for r in stats)
    reasons=[]
    if not enough_replicas:reasons.append('too_few_independent_replicas')
    if not transitions_pass:reasons.append('too_few_bidirectional_occupancy_transitions')
    if not states_pass:reasons.append('insufficient_empty_or_occupied_state_frames')
    if not plateau_pass:reasons.append('no_resolved_implied_timescale_plateau')
    if not ck_pass:reasons.append('held_out_chapman_kolmogorov_not_passed')
    # Passing diagnostics for whether any molecule is present does not validate a solute transport model.
    reasons.append('occupancy_states_not_validated_for_transport_TPT')
    return {'state_definition':'0=no_marker_in_site; 1=one_or_more_markers_in_site',
        'independent_replicas':len(replicas),'dt_ps':float(dt_ps),'replicas':stats,'implied_timescales':models,
        'its_relative_range':spread,'its_pass':bool(plateau_pass),'held_out_ck':ck,'ck_pass':bool(ck_pass),
        'occupancy_diagnostics_pass':len(reasons)==1,'status':'transport_sampling_not_established',
        'reasons':reasons,'rates':None,'tpt_flux':None}
