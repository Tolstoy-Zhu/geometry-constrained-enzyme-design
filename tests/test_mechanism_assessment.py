"""Geometric and stochastic fixtures test algorithms, not biochemical validation."""
import numpy as np
import pytest
from bur01_structmine.unified.channels import clearance_grid,grid_path,widest_path,continuous_path_clearance
from bur01_structmine.unified.sampling import transition_counts,transition_model,assess_occupancy
from bur01_structmine.unified.mechanism_assessment import AssessmentProtocol,summarize_channels


def test_continuous_edge_path_matches_reported_clearance():
    pytest.importorskip('scipy')
    xyz=np.array([[0.,0.,0.],[8.,0.,0.]])
    grid=clearance_grid(xyz,['C','C'],spacing_A=1.)
    original=widest_path(xyz,['C','C'],[4,0,0],spacing_A=1.,probe_radius_A=.5)
    assert grid_path(grid,[4,0,0],.5,continuous_edges=True)==original
    assert original['status']=='path_found'
    assert continuous_path_clearance(original['centerline_A'],xyz,['C','C'])==pytest.approx(original['bottleneck_A'])
    assert grid_path(grid,[4,0,0],8.)['status']=='no_accessible_seed'


def test_segment_clearance_detects_collision_missed_at_nodes():
    # Endpoints clear a unit-radius probe, but the segment crosses the carbon sphere.
    assert continuous_path_clearance([[-3,0,0],[3,0,0]],[[0,0,0]],['C'])==pytest.approx(-1.7)
    assert continuous_path_clearance([[-3,3,0],[3,3,0]],[[0,0,0]],['C'])==pytest.approx(1.3)


def test_path_search_rejects_an_edge_crossing_an_atom_between_clear_nodes():
    pytest.importorskip('scipy')
    grid=clearance_grid([[0,0,0]],['C'],spacing_A=2.)
    # A masked synthetic corridor isolates the edge tangent to a carbon sphere.
    shape=np.array([3,3,3]);ijk=np.stack(np.unravel_index(np.arange(27),shape),axis=1)
    points=np.array([-3.,-.5,-2.])+ijk*2
    clearance=np.full(27,-10.)
    for node in [13,22]:clearance[node]=np.linalg.norm(points[node])-1.7
    grid.update(shape=shape,ijk=ijk,points=points,clearance=clearance)
    assert grid_path(grid,points[13],.1)['status']=='path_found'
    assert grid_path(grid,points[13],.1,continuous_edges=True)['status']=='no_open_path_at_probe_radius'


@pytest.mark.parametrize('kwargs',[{'spacing_A':0.1},{'origin_fraction':1.},{'max_voxels':1}])
def test_channel_parameter_and_memory_limits(kwargs):
    with pytest.raises(ValueError):clearance_grid([[0,0,0]],['C'],**kwargs)


def test_no_transition_is_created_between_independent_replicas():
    counts=transition_counts([np.zeros(5),np.ones(5)],1)
    assert counts.tolist()==[[4,0],[0,4]]
    assert transition_model(counts,10)['status']=='disconnected'


def test_analytic_two_state_implied_timescale():
    model=transition_model([[90,10],[20,80]],10)
    assert model['eigenvalue']==pytest.approx(.7)
    assert model['timescale_ps']==pytest.approx(-10/np.log(.7))


def test_unvisited_state_does_not_produce_rate_or_ck_pass():
    replicas=[{'replica_id':str(i),'counts':np.zeros(1000)} for i in range(3)]
    report=assess_occupancy(replicas,10,AssessmentProtocol().model_dump())
    assert not report['its_pass'] and not report['ck_pass']
    assert report['rates'] is None and report['tpt_flux'] is None
    assert all(r['matrix'] is None for r in report['implied_timescales'])
    assert all(r['status']=='insufficient_state_counts' for r in report['held_out_ck'])


def test_sampled_markov_chain_can_pass_diagnostics_but_not_transport_validation():
    rng=np.random.default_rng(407);replicas=[]
    for i in range(3):
        flip=rng.random(50000)<.005
        states=np.logical_xor.accumulate(flip).astype(int)
        replicas.append({'replica_id':str(i),'counts':states})
    report=assess_occupancy(replicas,10,AssessmentProtocol().model_dump())
    assert report['its_pass'] and report['ck_pass'] and report['occupancy_diagnostics_pass']
    assert report['reasons']==['occupancy_states_not_validated_for_transport_TPT']
    assert report['rates'] is None
    assert {r['held_out_replica'] for r in report['held_out_ck']}=={'0','1','2'}


def sweep():
    protocol=AssessmentProtocol().model_dump()
    rows=[{'candidate_id':'example','replica_id':'r1','site':'A','frame':frame,'probe_radius_A':1.,
           'spacing_A':spacing,'origin_fraction':origin,'status':'path_found',
           'edge_clearance_pass':True,'continuous_bottleneck_A':1.5}
          for frame in [0,100] for spacing in protocol['spacings_A'] for origin in protocol['origin_fractions']]
    return protocol,rows


def test_grid_disagreement_is_reported_instead_of_averaged_away():
    protocol,rows=sweep()
    assert summarize_channels(rows,protocol)[0]['status']=='grid_agreement_on_sampled_frames'
    rows[0].update(status='no_open_path_at_probe_radius',edge_clearance_pass=False,continuous_bottleneck_A=None)
    result=summarize_channels(rows,protocol)[0]
    assert result['status']=='grid_sensitive' and result['grid_disagreement_fraction']==.5
    assert result['matched_open_frames']==1


def test_collision_and_missing_sweep_cannot_be_marked_passed():
    protocol,rows=sweep();rows[0]['edge_clearance_pass']=False
    assert summarize_channels(rows,protocol)[0]['status']=='edge_clearance_requires_review'
    with pytest.raises(ValueError,match='Incomplete'):summarize_channels(rows[:-1],protocol)


def test_no_path_keeps_missing_bottleneck():
    protocol,rows=sweep()
    for row in rows:row.update(status='no_accessible_seed',edge_clearance_pass=False,continuous_bottleneck_A=None)
    result=summarize_channels(rows,protocol)[0]
    assert result['status']=='no_paths_on_sampled_frames'
    assert result['matched_bottleneck_range_max_A'] is None


@pytest.mark.parametrize('changes',[{'spacings_A':[.8,.8]}, {'origin_fractions':[0.]},
    {'base_lag_frames':3},{'plateau_lags_frames':[100,200]}, {'ck_multipliers':[2,1]}])
def test_invalid_protocol_rejected_before_computing(changes):
    with pytest.raises(ValueError):AssessmentProtocol(**changes)
