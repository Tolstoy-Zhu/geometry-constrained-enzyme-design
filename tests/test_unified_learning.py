"""Numerical fixtures validate algorithms only; they are never exported as measured training data."""
import numpy as np
import pandas as pd
import pytest
pytest.importorskip("rdkit")
pytest.importorskip("sklearn")
from bur01_structmine.unified.datasets import canonicalize,grouped_split
from bur01_structmine.unified.learning import validate_splits,grouped_improvement
from bur01_structmine.unified.channels import widest_path


def test_kinetic_export_keeps_candidate_names_and_duplicate_sequence_aliases(tmp_path,monkeypatch):
    import json
    from bur01_structmine.unified import inference,learning
    monkeypatch.setattr(inference.subprocess,'run',lambda *args,**kwargs:None)
    def predict(folder,frame,features):
        return {'predictions':[{'enzyme_id':r.sequence_sha256,'reaction_id':r.substrate_id,
                               'endpoint':'Km','unit':'M','value':.001} for _,r in frame.iterrows()], 'excluded':[]}
    monkeypatch.setattr(learning,'predict_model',predict)
    text='>parent\nACDE\n>same_sequence_alias\nACDE\n>variant\nACDF\n'
    output=tmp_path/'result'
    inference.predict_kinetics(tmp_path/'model',text,[('ethanol','CCO')],output)
    report=json.loads((output/'predictions.json').read_text())
    assert len(report['predictions'])==2
    assert report['predictions'][0]['candidate_ids']==['parent','same_sequence_alias']
    assert report['predictions'][1]['candidate_ids']==['variant']
    assert all(r['substrate_names']==['ethanol'] and r['data_kind']=='prediction' for r in report['predictions'])
    assert (output/'inputs.fasta').read_text()==text

def test_promotion_also_requires_lower_primary_test_mae():
    from bur01_structmine.unified.learning import promotion_decision
    y=np.zeros(105);baseline=np.ones(105)
    candidate=np.array([1.2]*100+[0.0]*5)
    groups=np.array([0]*100+list(range(1,6)))
    report=promotion_decision("ExtraTrees",y,baseline,candidate,groups)
    assert report["sample_weighted_delta_MAE"]>0
    assert not report["promote"]

def test_installed_openmm_patch_version_is_used(tmp_path):
    pytest.importorskip("openmm")
    from bur01_structmine.unified.mechanism import run_prepared
    manifest=tmp_path/"unprepared.json";manifest.write_text("{}")
    with pytest.raises(ValueError,match="Prepared system is blocked"):
        run_prepared(manifest,tmp_path/"run",production_ns=.001)

def test_new_substrate_challenge_is_distinct_and_never_used_for_selection():
    from bur01_structmine.unified.learning import evaluation_masks
    frame=pd.DataFrame([{"sequence_sha256":split+str(i),"substrate_id":"s"+str(j),"split":split,"log10_value":i+j}
        for split in ("train","validation","calibration","test") for i in range(3) for j in range(20)])
    masks,protocol=evaluation_masks(frame)
    held=set(protocol["reserved_substrates"])
    for key in ("train","validation","calibration"):
        assert not set(frame.loc[masks[key],"substrate_id"]) & held
    assert masks["new_substrate"].any() and not (masks["new_substrate"] & masks["test"]).any()
    assert set(frame.loc[masks["new_substrate"],"sequence_sha256"]) <= set(frame.loc[masks["train"],"sequence_sha256"])
    changed=frame.copy();changed["log10_value"]=999
    assert evaluation_masks(changed)[1]==protocol

def test_units_censor_conflicts_and_simulations():
    base={"record_id":"one","sequence":"ACDE","smiles":"CCO","endpoint":"kcat","value":60.,
          "unit":"min^-1","source":"fixture","data_kind":"experimental","qualifier":"="}
    frame=pd.DataFrame([base,{**base,"record_id":"two","data_kind":"simulated"},
                        {**base,"record_id":"three","qualifier":"<","detection_limit":1},
                        {**base,"value":120}])
    clean,rejected=canonicalize(frame)
    assert clean.value.tolist()==[1,2]
    assert len(rejected)==2

def test_parent_and_repeated_assays_stay_together():
    frame=pd.DataFrame({"sequence_sha256":[str(i) for i in range(40)],"parent_id":["shared"]*3+[f"p{i}" for i in range(37)],
                        "measurement_group":[f"r{i}" for i in range(40)],"data_kind":"experimental"})
    groups=pd.DataFrame({"sequence_sha256":frame.sequence_sha256,"group_id":frame.sequence_sha256})
    split=grouped_split(frame,groups)
    validate_splits(split)
    assert split.iloc[:3].split.nunique()==1
    assert set(split.split)=={"train","validation","calibration","test"}
    changed=split.copy();changed.loc[1,"split"]="test" if split.loc[0,"split"]!="test" else "train"
    with pytest.raises(ValueError,match="Leakage"):validate_splits(changed)

def test_promotion_requires_grouped_improvement():
    y=np.arange(10.)
    assert grouped_improvement(y,y+2,y+.2,np.arange(10))["promote"]
    assert not grouped_improvement(y,y+.1,y+2,np.arange(10))["promote"]
    assert not grouped_improvement(y,y+2,y+.1,np.zeros(10))["promote"]

def test_channels_are_computed_from_current_coordinates():
    # Simple geometric fixture, not a protein or mechanistic validation.
    result=widest_path(np.array([[0,0,0],[8,0,0]]),["C","C"],[4,0,0],spacing_A=1.,probe_radius_A=.5)
    assert result["status"]=="path_found"
    assert result["bottleneck_A"]>0
    assert len(result["centerline_A"])>1
