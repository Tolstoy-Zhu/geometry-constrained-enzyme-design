import copy,json
import pytest
from bur01_structmine.unified.mechanism_results import summarize_batch,load_result
from bur01_structmine.unified.datasets import file_hash


def analysis(value):
    return {'frames':1000,'radius_A':5.,'channel':{'frame_results':[
        {'status':'no_accessible_seed','bottleneck_A':None},
        {'status':'path_found','bottleneck_A':1.5}]},
        'observables':{'water':{'mean_molecules_in_region':value,'fraction_frames_occupied':.8,'entry_events':10,
            'residence_events':[{'observed_duration_ps':500.,'left_censored':True,'right_censored':True},
                {'observed_duration_ps':10.,'left_censored':False,'right_censored':False}]}}}


def batch():
    return {'replicas_completed':3,'production_ns':30.,'reports':[
        {'candidate_id':'example','replica_id':f'example_r{i}','production_ns':10.,'expected_regions':['A','B'],
         'regions':{'A':analysis(i),'B':analysis(i+2)}} for i in range(1,4)]}


def test_sites_are_averaged_within_each_independent_replica():
    result=summarize_batch(batch())
    assert len(result['sites'])==6
    assert result['candidates'][0]['replicates']==3
    assert result['candidates'][0]['mean_molecules_mean']==3
    assert result['candidates'][0]['mean_molecules_min']==2
    assert result['candidates'][0]['mean_molecules_max']==4
    assert result['analyzed_replicas']==3
    assert result['quality_passed']==0


def test_incomplete_site_analysis_cannot_enter_candidate_average():
    data=batch();del data['reports'][1]['regions']['B']
    result=summarize_batch(data)
    assert result['analyzed_replicas']==2
    assert result['candidates'][0]['replicates']==2
    assert len(result['sites'])==5


def test_completion_without_expected_sites_is_unverified():
    data=batch()
    for row in data['reports']:del row['expected_regions']
    result=summarize_batch(data)
    assert result['analyzed_replicas']==0
    assert result['candidates']==[]


def test_censored_events_and_missing_channels_are_not_zero_measurements():
    data=batch();result=summarize_batch(data);row=result['sites'][0]
    assert row['complete_event_observed_median_ps']==10
    assert row['left_censored_events']==row['right_censored_events']==1
    assert row['path_found_fraction']==.5
    assert row['path_bottleneck_mean_A']==1.5
    data['reports'][0]['regions']['A']['channel']['frame_results']=[{'status':'no_accessible_seed','bottleneck_A':None}]
    assert summarize_batch(data)['sites'][0]['path_bottleneck_mean_A'] is None


def test_failed_quality_record_remains_failed_even_with_passed_label():
    quality={'reports':[{'replica_id':'example_r1','status':'passed','issues':['broken_bond']}]}
    result=summarize_batch(batch(),quality)
    assert result['quality_passed']==0
    assert result['quality'][0]['quality_status']=='requires_review'
    assert result['quality'][1]['quality_status']=='unreviewed'


@pytest.mark.parametrize('change',[lambda b:b['reports'].append(copy.deepcopy(b['reports'][0])),
    lambda b:b.update(production_ns=120),
    lambda b:b['reports'][0]['regions']['A']['observables']['water'].update(mean_molecules_in_region=float('nan')),
    lambda b:b['reports'][0]['regions']['A']['observables']['water'].update(fraction_frames_occupied=1.2),
    lambda b:b['reports'][0]['regions']['A'].update(radius_A=4.)])
def test_invalid_or_incommensurate_reports_fail_closed(change):
    data=batch();change(data)
    with pytest.raises(ValueError):summarize_batch(data)


def test_quality_is_bound_to_batch_and_each_trajectory_report(tmp_path):
    data=batch();path=tmp_path/'batch_report.json';path.write_text(json.dumps(data))
    (tmp_path/'batch_identity.json').write_text(json.dumps({'manifest_sha256':'frozen-batch'}))
    rows=[]
    for replica in data['reports']:
        folder=tmp_path/replica['replica_id'];folder.mkdir()
        source=folder/'trajectory_report.json';source.write_text(json.dumps({'replica':replica['replica_id']}))
        (folder/'preflight.json').write_text(json.dumps({'contract':{'construct_sequences':{'A':'AC','B':'AC'}}}))
        rows.append({'replica_id':replica['replica_id'],'status':'numerical_checks_passed_descriptors_available','issues':[],
                     'identity':{'trajectory_report_sha256':file_hash(source)}})
    qpath=tmp_path/'quality.json';qpath.write_text(json.dumps({'identity':{'batch_sha256':'frozen-batch'},'reports':rows}))
    record={'report_path':str(qpath),'sha256':file_hash(qpath)}
    assert load_result(path,record)['quality_passed']==3
    source.write_text('changed')
    with pytest.raises(ValueError,match='does not match'):load_result(path,record)
    qpath.write_text('{}')
    with pytest.raises(ValueError,match='checksum'):load_result(path,record)


def test_unsafe_replica_directory_cannot_be_read(tmp_path):
    data=batch();data['reports'][0]['replica_id']='../elsewhere'
    path=tmp_path/'batch_report.json';path.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='Unsafe'):load_result(path)


def test_different_candidate_site_radii_cannot_be_compared():
    data=batch();second=copy.deepcopy(data['reports'])
    for row in second:
        row['candidate_id']='second';row['replica_id']='second_'+row['replica_id']
        for region in row['regions'].values():region['radius_A']=6.
    data['reports']+=second;data.update(production_ns=60.,replicas_completed=6)
    with pytest.raises(ValueError,match='across candidates'):summarize_batch(data)


def test_entry_counts_are_not_compared_across_unequal_durations():
    data=batch();data['reports'][0]['production_ns']=20.;data['production_ns']=40.
    with pytest.raises(ValueError,match='durations'):summarize_batch(data)
