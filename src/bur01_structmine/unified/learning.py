"""Independent endpoint regression, held-out calibration and grouped model promotion."""
from __future__ import annotations
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
import pandas as pd
from .datasets import file_hash
from ..property_assets import ASSETS

def metrics(y,pred):
    from scipy.stats import spearmanr
    error=np.asarray(pred)-np.asarray(y)
    correlation=float(spearmanr(y,pred).statistic) if len(y)>=3 and np.ptp(y)>0 and np.ptp(pred)>0 else None
    return {"n":len(y),"MAE":float(np.mean(np.abs(error))),"RMSE":float(np.sqrt(np.mean(error**2))),
            "Spearman":correlation if correlation is None or math.isfinite(correlation) else None}

def feature_matrix(frame,feature_dir,conditions=False):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    feature_dir=Path(feature_dir)
    manifest=json.loads((feature_dir/"manifest.json").read_text())
    if manifest.get("model_sha256")!=ASSETS["esm/model.safetensors"]:
        raise ValueError("Features must use the repository's pinned ESM-2 weights")
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    cache={};out=[];retained=[];excluded=[]
    for index,row in frame.iterrows():
        key=row.sequence_sha256
        path=feature_dir/(key+".npy")
        if key not in manifest["files"] or not path.exists():
            excluded.append({"record_id":row.record_id,"reason":"embedding_unavailable"});continue
        if conditions and (pd.isna(row.get("pH")) or pd.isna(row.get("temperature_K"))):
            excluded.append({"record_id":row.record_id,"reason":"conditions_unreported"});continue
        if key not in cache:
            if file_hash(path)!=manifest["files"][key]:raise ValueError("Embedding checksum mismatch")
            vector=np.load(path,allow_pickle=False)
            if vector.shape!=(1280,) or not np.isfinite(vector).all():raise ValueError("Expected ESM-2 650M 1280-vector")
            cache[key]=vector
        mol=Chem.MolFromSmiles(row.smiles)
        if mol is None:raise ValueError("Invalid frozen dataset substrate")
        parts=[cache[key],generator.GetFingerprintAsNumPy(mol)]
        if conditions:parts.append(np.array([float(row.pH),float(row.temperature_K)]))
        out.append(np.concatenate(parts).astype(np.float32));retained.append(index)
    if not out:raise ValueError("No eligible rows with verified embeddings")
    return np.stack(out),frame.loc[retained].reset_index(drop=True),excluded,manifest

def validate_splits(frame):
    required={"train","validation","calibration","test"}
    if set(frame.split)!=required:raise ValueError("All four frozen splits are required")
    for key in ["sequence_sha256","group_id","parent_id","measurement_group"]:
        if key not in frame:continue
        sub=frame[frame[key].notna() & frame[key].astype(str).str.strip().ne("")]
        if (sub.groupby(key).split.nunique()>1).any():raise ValueError("Leakage across splits: "+key)
    if (frame.data_kind!="experimental").any():raise ValueError("Only measured labels may enter training")

def grouped_improvement(y,baseline,candidate,groups,repetitions=1000,seed=2026):
    delta=pd.DataFrame({"group":groups,"delta":np.abs(candidate-y)-np.abs(baseline-y)}).groupby("group").delta.mean().to_numpy()
    if len(delta)<5:return {"status":"insufficient_test_groups","interval_95":None,"promote":False}
    rng=np.random.default_rng(seed)
    means=np.array([rng.choice(delta,len(delta),replace=True).mean() for _ in range(repetitions)])
    bounds=np.quantile(means,[.025,.975]).tolist()
    return {"status":"assessed","delta_MAE":float(delta.mean()),"interval_95":bounds,
            "promote":bool(bounds[1]<0),"bootstrap_unit":"sequence_parent_group","repetitions":repetitions}

def promotion_decision(selected, y, baseline, candidate, groups, seed=2026):
    report=grouped_improvement(y,baseline,candidate,groups,seed=seed)
    sample_delta=float(np.mean(np.abs(candidate-y))-np.mean(np.abs(baseline-y)))
    return {**report,"sample_weighted_delta_MAE":sample_delta,
            "promote":bool(selected=="ExtraTrees" and sample_delta<0 and report["promote"])}

def evaluation_masks(frame,seed=2026):
    """Reserve substrate identities without consulting labels or test performance.

    Homology groups retain their primary 70/10/10/10 allocation. Ten percent of
    substrate identities are withheld from fitting, selection and calibration.
    Their rows in training groups form a separate seen-enzyme/new-substrate
    challenge; they never contribute to the main four-split fit or calibration.
    """
    substrates=np.array(sorted(frame.substrate_id.unique()),dtype=object)
    rng=np.random.default_rng(seed+1);rng.shuffle(substrates)
    count=max(1,round(.1*len(substrates))) if len(substrates)>1 else 0
    held=set(substrates[:count])
    reserved=frame.substrate_id.isin(held).to_numpy()
    masks={n:frame.split.eq(n).to_numpy() & ~reserved for n in ("train","validation","calibration")}
    masks["test"]=frame.split.eq("test").to_numpy()
    seen_sequences=set(frame.loc[masks["train"],"sequence_sha256"])
    masks["new_substrate"]=frame.split.eq("train").to_numpy() & reserved & frame.sequence_sha256.isin(seen_sequences).to_numpy()
    protocol={"version":"substrate_identity_holdout_v1","seed":seed+1,"fraction":.1,
              "reserved_substrates":sorted(held),"counts":{n:int(m.sum()) for n,m in masks.items()},
              "selection_uses_labels":False,"primary_homology_allocation_unchanged":True}
    return masks,protocol

def train_models(dataset,features,output,conditions=False,seed=2026):
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler,FunctionTransformer
    from sklearn.linear_model import Ridge
    from sklearn.ensemble import ExtraTreesRegressor
    import joblib
    dataset=Path(dataset);output=Path(output)
    if output.exists():raise FileExistsError(output)
    metadata=json.loads((dataset/"manifest.json").read_text())
    if metadata.get("status")!="ready_for_features":raise ValueError("Dataset has no verified homology split")
    if file_hash(dataset/"split.csv")!=metadata["split_sha256"]:raise ValueError("Frozen split checksum mismatch")
    original=pd.read_csv(dataset/"split.csv")
    validate_splits(original)
    X,frame,excluded,feature_info=feature_matrix(original,features,conditions)
    validate_splits(frame)
    endpoint=metadata["endpoint"]
    if endpoint not in {"kcat","Km"}:raise ValueError("No trained task for this endpoint")
    if set(frame.endpoint)!={endpoint}:raise ValueError("A kinetic model must contain exactly one endpoint")
    y=frame.log10_value.to_numpy(float)
    masks,challenge_protocol=evaluation_masks(frame,seed)
    if min(int(masks[n].sum()) for n in ("train","validation","calibration","test"))<2:raise ValueError("Insufficient observations in frozen splits")
    train,val,cal,test=[masks[n] for n in ["train","validation","calibration","test"]]
    output.mkdir(parents=True)
    (output/"evaluation_protocol.json").write_text(json.dumps(challenge_protocol,indent=2),encoding="utf-8")
    candidates=[];selection=[]
    for alpha in (1.,10.,100.):
        print(f"Fitting Ridge alpha={alpha}",flush=True)
        model=make_pipeline(FunctionTransformer(np.asarray,kw_args={"dtype":"float64"}),StandardScaler(),Ridge(alpha=alpha))
        model.fit(X[train],y[train]);score=metrics(y[val],model.predict(X[val]))
        candidates.append(("Ridge",{"alpha":alpha},model,score["MAE"]))
        print(f"Validation MAE={score['MAE']:.6f}",flush=True)
    for leaf in (2,5,10):
        print(f"Fitting ExtraTrees min_samples_leaf={leaf}",flush=True)
        model=ExtraTreesRegressor(n_estimators=200,min_samples_leaf=leaf,max_features=.5,
                                 n_jobs=max(1,int(os.environ.get("AQUA_CPU_THREADS","8"))),random_state=seed)
        model.fit(X[train],y[train]);score=metrics(y[val],model.predict(X[val]))
        candidates.append(("ExtraTrees",{"min_samples_leaf":leaf,"trees":200,"max_features":.5},model,score["MAE"]))
        print(f"Validation MAE={score['MAE']:.6f}",flush=True)
    ridge=min((c for c in candidates if c[0]=="Ridge"),key=lambda c:c[3])
    tree=min((c for c in candidates if c[0]=="ExtraTrees"),key=lambda c:c[3])
    selected=min((ridge,tree),key=lambda c:c[3])
    b=ridge[2].predict(X[test]);t=tree[2].predict(X[test]);p=selected[2].predict(X[test])
    improvement=promotion_decision(selected[0],y[test],b,t,frame.loc[test,"group_id"].to_numpy(),seed=seed)
    promotion=improvement["promote"]
    residual=pd.DataFrame({"group":frame.loc[cal,"group_id"].to_numpy(),
                           "residual":np.abs(selected[2].predict(X[cal])-y[cal])}).groupby("group").residual.max().to_numpy()
    radius=None
    if len(residual)>=20:
        level=min(1.,math.ceil((len(residual)+1)*.95)/len(residual))
        radius=float(np.quantile(residual,level,method="higher"))
    training_substrates=set(frame.loc[train,"substrate_id"])
    cold=~frame.loc[test,"substrate_id"].isin(training_substrates).to_numpy()
    reports={"Ridge":metrics(y[test],b),"ExtraTrees":metrics(y[test],t),
             "training_mean":metrics(y[test],np.full(test.sum(),float(y[train].mean()))),
             "new_enzyme":metrics(y[test],p),
             "new_enzyme_and_substrate":metrics(y[test][cold],p[cold]) if cold.any() else {"n":0,"status":"no_eligible_test_rows"}}
    challenge=masks["new_substrate"]
    reports["new_substrate"]={**(metrics(y[challenge],selected[2].predict(X[challenge])) if challenge.any() else {"n":0,"status":"no_eligible_test_rows"}),
                               "cohort":"seen_exact_enzyme_reserved_substrate"}
    reports["new_enzyme_seen_substrate"]=metrics(y[test][~cold],p[~cold]) if (~cold).any() else {"n":0}
    import sklearn,scipy,rdkit
    recipe={"version":"esm_morgan_ridge_extratrees_v2","seed":seed,"conditions":conditions,"ridge_precision":"float64",
            "split_sha256":metadata["split_sha256"],"features_sha256":file_hash(Path(features)/"manifest.json"),
            "evaluation_protocol_sha256":file_hash(output/"evaluation_protocol.json"),
            "versions":{"sklearn":sklearn.__version__,"numpy":np.__version__,"scipy":scipy.__version__,"rdkit":rdkit.__version__}}
    recipe_hash=hashlib.sha256(json.dumps(recipe,sort_keys=True).encode()).hexdigest()
    model_id=f"aqua_{endpoint}_{recipe_hash[:16]}"
    joblib.dump(selected[2],output/"model.joblib")
    pd.DataFrame({"record_id":frame.loc[test,"record_id"].to_numpy(),"group_id":frame.loc[test,"group_id"].to_numpy(),
                  "observed_log10":y[test],"predicted_log10":p,"ridge_log10":b,"extratrees_log10":t}).to_csv(output/"test_predictions.csv",index=False)
    report={"schema_version":"1.0.0","model_id":model_id,"endpoint":endpoint,
            "unit":"s^-1" if endpoint=="kcat" else "M","target_transform":"log10","model":selected[0],
            "parameters":selected[1],"status":"validated_upgrade" if promotion else "experimental",
            "recipe":recipe,"promotion":{**improvement,"promote":promotion},"metrics":reports,
            "hyperparameter_selection":[{"model":n,"parameters":params,"validation_MAE":score} for n,params,_,score in candidates],
            "conditions":conditions,"n":len(frame),"excluded":excluded,
            "source_sha256":metadata["source_sha256"],"split_sha256":metadata["split_sha256"],
            "feature_manifest_sha256":file_hash(Path(features)/"manifest.json"),
            "encoder_sha256":feature_info["model_sha256"],"model_sha256":file_hash(output/"model.joblib"),
            "calibration":{"radius_log10":radius,"independent_groups":len(residual),
                           "status":"calibrated_in_domain" if radius is not None else "insufficient_calibration_groups",
                           "test_coverage":float(np.mean(np.abs(p-y[test])<=radius)) if radius is not None else None},
            "applicability":{"max_sequence_length":feature_info.get("max_length",1022),
                             "training_sequence_ids":sorted(frame.loc[train,"sequence_sha256"].unique()),
                             "training_substrate_ids":sorted(training_substrates),
                             "unreported_conditions":"not_inferred","kinetic_predictions_are_not_esterification_yields":True}}
    (output/"model_card.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
    return {"output":str(output),"model_id":model_id,"status":report["status"],"metrics":reports}

def predict_model(model_dir,frame,features):
    """Only load a local model whose serialized checksum matches its model card."""
    import joblib
    from .contracts import PredictionResult
    folder=Path(model_dir);card=json.loads((folder/"model_card.json").read_text())
    if file_hash(folder/"model.joblib")!=card["model_sha256"]:raise ValueError("Model checksum mismatch")
    X,eligible,excluded,_=feature_matrix(frame,features,card["conditions"])
    model=joblib.load(folder/"model.joblib")
    values=model.predict(X);out=[]
    radius=card["calibration"]["radius_log10"]
    for (_,row),value in zip(eligible.iterrows(),values):
        known=row.substrate_id in card["applicability"]["training_substrate_ids"]
        sequence_known=row.sequence_sha256 in card['applicability']['training_sequence_ids']
        assessed=known and sequence_known
        out.append(PredictionResult(enzyme_id=row.sequence_sha256,reaction_id=row.substrate_id,
                    endpoint=card["endpoint"],unit=card["unit"],value=float(10**value),
                    model_id=card["model_id"],model_version="1",
                    status="predicted" if assessed else "out_of_domain",
                    interval=(float(10**(value-radius)),float(10**(value+radius))) if assessed and radius is not None else None,
                    applicability={"substrate_seen":known,"model_status":card["status"],
                                   "enzyme_similarity_assessment":"seen_sequence" if sequence_known else "unassessed_sequence"},
                    uncertainty_status=card["calibration"]["status"] if assessed else "domain_not_established"
                    ).model_dump(mode="json"))
    return {"predictions":out,"excluded":excluded}
