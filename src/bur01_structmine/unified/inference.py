"""Apply a registered kinetic predictor to actual input sequences and molecules."""
from __future__ import annotations
import hashlib,json,os,subprocess,sys
from pathlib import Path

def predict_kinetics(model_dir,text,substrates,output,device="cpu",temperature_K=None,pH=None):
    import pandas as pd
    from rdkit import Chem
    from ..property_fasta import collect_fasta
    from ..property_assets import asset_root,runtime_python
    from .learning import predict_model
    sequences,_,_,sources=collect_fasta(text,False)
    names={key:[] for key in sequences}
    for record in sources['FASTA']['records']:
        names[record['sequence_sha256']].append(record['id'])
    substrate_names={}
    rows=[]
    for name,smiles in substrates:
        mol=Chem.MolFromSmiles(smiles)
        if mol is None or len(Chem.GetMolFrags(mol))!=1:raise ValueError("Invalid individual substrate: "+name)
        canonical=Chem.MolToSmiles(mol,isomericSmiles=True)
        substrate_id=hashlib.sha256(canonical.encode()).hexdigest()
        substrate_names.setdefault(substrate_id,[]).append(name)
        for key,sequence in sequences.items():
            rows.append({"record_id":key+":"+name,"sequence":sequence,"sequence_sha256":key,
                         "smiles":canonical,"substrate_id":substrate_id,
                         "substrate_name":name,"temperature_K":temperature_K,"pH":pH})
    if not rows:raise ValueError("Supply at least one substrate")
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    directory=output/"features";directory.mkdir(exist_ok=True)
    inp=directory/"sequences.json";inp.write_text(json.dumps(sequences),encoding="utf-8")
    interpreter=os.environ.get("AQUA_FEATURE_PYTHON",str(runtime_python()))
    if not Path(interpreter).is_file():interpreter=sys.executable
    subprocess.run([interpreter,str(Path(__file__).with_name("features_worker.py")),str(inp),
                    str(directory),str(asset_root()/"esm"),"--device",device],check=True)
    result=predict_model(model_dir,pd.DataFrame(rows),directory)
    for prediction in result['predictions']:
        prediction['candidate_ids']=names[prediction['enzyme_id']]
        prediction['substrate_names']=substrate_names[prediction['reaction_id']]
        prediction['data_kind']='prediction'
    result['input_records']=sources['FASTA']['records']
    (output/'inputs.fasta').write_text(text,encoding='utf-8')
    (output/"predictions.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
    pd.DataFrame(result["predictions"]).to_csv(output/"predictions.csv",index=False)
    return {"output":str(output),"prediction_count":len(result["predictions"]),"data_kind":"prediction"}
