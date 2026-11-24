import json
import numpy as np
import pandas as pd
import pytest
from bur01_structmine.unified.families import get_profile,reference,assess_structure,map_sequence,scope_scores
from bur01_structmine.structure_search import load_ca_trace
from bur01_structmine.unified.mechanism_plan import MechanismBatch,map_construct,select_candidates,pareto_front
from bur01_structmine.unified.mechanism_analysis import minimum_image,residence_events
from bur01_structmine.unified.mechanism import inspect_prepared

@pytest.mark.parametrize("family,count",[("bur01",3),("calb",3),("gh1_bgl",2)])
def test_reference_family_catalytic_roles(family,count):
    p=get_profile(family);trace=load_ca_trace(reference(p),chain_id=p.reference_chain)
    checked=assess_structure(trace,p)
    assert checked["catalytic_identity_pass"] and checked["fold_pass"]
    assert len(checked["sites"])==count
    mapped=map_sequence(trace.sequence,family)
    assert mapped["complete"] and len(mapped["sites"])==count
    if family=="gh1_bgl":
        assert set(mapped["sites"])=={"E166","E355"}
        assert not p.protect_serine_motif

def test_non_bur01_never_receives_wase():
    frame=pd.DataFrame({"wase_screen_score":[99.],"wase_tier":["A"],"target":["one"]})
    result=scope_scores(frame,"gh1_bgl")
    assert result.wase_screen_score.isna().all()
    assert result.wase_tier.isna().all()
    with pytest.raises(ValueError,match="Unknown enzyme family"):get_profile("unknown")

def test_construct_mapping_is_separate_and_unambiguous():
    result=map_construct("ACDE","ACFE","MMACDEKK")
    assert result["sequence"]=="MMACFEKK"
    assert result["mutations"]==[{"core_position":3,"construct_position":5,"parent":"D","candidate":"F"}]
    with pytest.raises(ValueError):map_construct("AC","AE","ACAC")
    with pytest.raises(ValueError):map_construct("AC","ACE","MMAC")

def test_120ns_budget_and_independent_seeds():
    batch=MechanismBatch(batch_id="pilot",candidates=[{"candidate_id":str(i)} for i in range(4)],reaction_id="r")
    assert batch.production_ns*batch.replicas*len(batch.candidates)==120
    assert len(set(batch.seeds))==12
    with pytest.raises(ValueError,match="budget"):MechanismBatch(**{**batch.model_dump(),"production_ns":20})

def test_shortlisting_uses_distinct_structurally_reviewed_sequences():
    frame=pd.DataFrame({"candidate_id":["a","b","c","d"],"sequence_sha256":["a","b","c","d"],
                        "structure_review_pass":[True,True,True,False],"sol":[4,3,2,99],"tm":[4,3,2,99],"pair":[4,3,2,99]})
    selected=select_candidates(frame,{"candidate_id":"parent","sequence_sha256":"p"},["sol","tm","pair"])
    assert [r["candidate_id"] for r in selected]==["parent","a","b","c"]

def test_pareto_respects_units_and_missing_values():
    frame=pd.DataFrame({"sol":[1,2,1,np.nan],"tm":[400,300,350,999]})
    result=pareto_front(frame,{"sol":"max","tm":"max"})
    assert result.pareto_front.tolist()==[True,True,False,False]

def test_periodic_box_and_censored_residences():
    assert np.allclose(minimum_image([9.5,0,0],np.eye(3)*10),[-.5,0,0])
    events=residence_events([True,True,False,True,True],[0,2,4,6,8])
    assert events[0]["left_censored"] and not events[0]["right_censored"]
    assert events[1]["right_censored"] and not events[1]["left_censored"]
    assert events[0]["observed_duration_ps"]==2

def test_old_benchmark_pass_is_not_production_release(tmp_path):
    path=tmp_path/"old.json"
    path.write_text(json.dumps({"identity_gate":"PASS","assembly_validation":"PASS","parameter_validation":"PASS",
                                 "equilibration_validation":"PASS"}))
    result=inspect_prepared(path)
    assert result["status"]=="blocked"
    assert "confirmed_full_construct_sequences" in result["missing"]
    assert "construct_validation" in result["missing"]
