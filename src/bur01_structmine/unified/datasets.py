"""Frozen public measurements, explicit units and leakage-aware grouped splits."""
from __future__ import annotations
import hashlib
import json
import math
import os
import subprocess
import tempfile
import zipfile
from pathlib import Path
import pandas as pd
import requests
from .contracts import AA

CATAPRO_URL="https://zenodo.org/records/14894710/files/catapro_data.zip?download=1"
CATAPRO_MD5="2f8411d35fba9f5a9b3069dac59d3010"

def file_hash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()

def fetch_catapro(destination):
    path=Path(destination);path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        part=path.with_suffix(".part")
        with requests.get(CATAPRO_URL,stream=True,timeout=(20,120)) as response:
            response.raise_for_status()
            with part.open("wb") as f:
                for block in response.iter_content(1024*1024): f.write(block)
        part.replace(path)
    digest=hashlib.md5(path.read_bytes()).hexdigest()
    if digest!=CATAPRO_MD5: raise ValueError("CataPro published archive checksum mismatch")
    return {"source":"CataPro","url":CATAPRO_URL,"sha256":file_hash(path),"published_md5":digest}

def fetch_sabio(query,destination):
    """Freeze the current export API's paginated JSON and companion SBML."""
    from .sabio import fetch
    return fetch(query,destination)

def catapro_rows(path,endpoint):
    if endpoint not in {"kcat","Km"}: raise ValueError("CataPro supports separate kcat and Km tasks")
    name="kcat-data_0.4simi-10fold.csv" if endpoint=="kcat" else "Km-data_0.4simi-10fold.csv"
    with zipfile.ZipFile(path) as archive:
        frame=pd.read_csv(archive.open("catapro_data/datasets/"+name))
    value="kcat(s^-1)" if endpoint=="kcat" else "Km(M)"
    return pd.DataFrame({
        "record_id":frame.iloc[:,0].astype(str),"sequence":frame.Sequence,"smiles":frame.Smiles,
        "endpoint":endpoint,"value":frame[value],"unit":"s^-1" if endpoint=="kcat" else "M",
        "source":"https://doi.org/10.5281/zenodo.14894710","parent_id":frame.UniProtID.fillna(""),
        "family":frame.EC.fillna(""),"organism":frame.Organism.fillna(""),
        "enzyme_type":frame.EnzymeType,"source_fold":frame.fold,"data_kind":"experimental",
        "qualifier":"=","temperature_K":None,"pH":None,"conditions":"not_reported",
        "measurement_group":frame.iloc[:,0].astype(str),"method":"source_dataset",
        "provenance_resolution":"curated_dataset_row_not_original_assay"
    })

def canonicalize(raw):
    from rdkit import Chem
    accepted=[];rejected=[]
    factors={"kcat":{"s^-1":1.,"1/s":1.,"min^-1":1/60},
             "Km":{"M":1.,"mM":1e-3,"uM":1e-6,"µM":1e-6}}
    for index,row in raw.iterrows():
        item=row.where(pd.notna(row),None).to_dict()
        try:
            seq=str(item.get("sequence") or "").strip().upper()
            if not seq or set(seq)-AA: raise ValueError("invalid_sequence")
            if item.get("data_kind")!="experimental": raise ValueError("not_experimental")
            if item.get("qualifier","=")!="=": raise ValueError("censored_or_unmeasured")
            endpoint=item["endpoint"];unit=item["unit"];value=float(item["value"])
            if endpoint not in factors or unit not in factors[endpoint]: raise ValueError("unsupported_endpoint_or_unit")
            if not math.isfinite(value) or value<=0: raise ValueError("nonpositive_or_invalid_kinetic_value")
            mol=Chem.MolFromSmiles(str(item.get("smiles") or ""))
            if mol is None or not mol.GetNumAtoms() or len(Chem.GetMolFrags(mol))!=1:
                raise ValueError("invalid_or_disconnected_substrate")
            smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
            digest=hashlib.sha256(seq.encode()).hexdigest()
            item.update(sequence=seq,sequence_sha256=digest,smiles=smiles,
                        substrate_id=hashlib.sha256(smiles.encode()).hexdigest(),
                        value=value*factors[endpoint][unit],unit="s^-1" if endpoint=="kcat" else "M",
                        log10_value=math.log10(value*factors[endpoint][unit]))
            item["record_id"]=str(item.get("record_id") or f"row_{index}")
            accepted.append(item)
        except (ValueError,KeyError,TypeError) as exc:
            rejected.append({"row":int(index)+2,"record_id":item.get("record_id"),"reason":str(exc)})
    frame=pd.DataFrame(accepted)
    if frame.empty: return frame,pd.DataFrame(rejected)
    # Preserve discordant values and assay conditions. Only literal duplicate evidence is removed.
    # Equal values from distinct conditions or replicates are separate measurements.
    duplicates=frame.duplicated()
    for row in frame.loc[duplicates].to_dict("records"):
        rejected.append({"record_id":row["record_id"],"reason":"literal_duplicate_evidence"})
    frame=frame.loc[~duplicates].reset_index(drop=True)
    return frame,pd.DataFrame(rejected)

class UnionFind:
    def __init__(self,items): self.parent={x:x for x in items}
    def find(self,x):
        while self.parent[x]!=x:
            self.parent[x]=self.parent[self.parent[x]];x=self.parent[x]
        return x
    def union(self,a,b):
        ra,rb=self.find(a),self.find(b)
        if ra!=rb:self.parent[max(ra,rb)]=min(ra,rb)

def mmseqs_groups(frame,output,executable="mmseqs"):
    """All-v-all edges, both-sequence coverage, then connected components."""
    directory=Path(output);directory.mkdir(parents=True,exist_ok=True)
    sequences=frame[["sequence_sha256","sequence"]].drop_duplicates()
    fasta=directory/"sequences.fasta"
    fasta.write_text("".join(f">{r.sequence_sha256}\n{r.sequence}\n" for r in sequences.itertuples()),encoding="ascii")
    edges=directory/"homology.tsv"
    command=[str(executable),"easy-search",str(fasta),str(fasta),str(edges),str(directory/"mmseqs_tmp"),
             "--min-seq-id","0.4","-c","0.8","--cov-mode","0","--max-seqs",str(len(sequences)),
             "-s","7.5","--threads",str(max(1,int(os.environ.get("AQUA_CPU_THREADS","8")))),
             "--format-output","query,target,fident,qcov,tcov"]
    result=subprocess.run(command,capture_output=True)
    if result.returncode:
        raise RuntimeError("MMseqs clustering failed: "+result.stderr.decode("utf-8",errors="replace")[-1500:])
    uf=UnionFind(sequences.sequence_sha256)
    for line in edges.read_text().splitlines():
        a,b,identity,qcov,tcov=line.split("\t")
        if float(identity)>=.4 and min(float(qcov),float(tcov))>=.8: uf.union(a,b)
    grouping=pd.DataFrame([{"sequence_sha256":key,"group_id":uf.find(key)} for key in uf.parent])
    grouping.to_csv(directory/"groups.tsv",sep="\t",index=False)
    info={"algorithm":"MMseqs2 all-v-all connected components","identity":.4,"coverage":.8,
          "command":command,"edges_sha256":file_hash(edges),'groups_sha256':file_hash(directory/'groups.tsv')}
    (directory/'groups.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    return grouping,info

def grouped_split(frame,groups,seed=2026):
    import numpy as np
    mapping=dict(zip(groups.sequence_sha256,groups.group_id.astype(str)))
    if set(frame.sequence_sha256)-set(mapping): raise ValueError("Cluster assignment missing sequences")
    uf=UnionFind(set(mapping.values()))
    # Keep WT/mutants and source-associated repetitions together in addition to homology.
    for column in ("parent_id","measurement_group"):
        if column not in frame: continue
        for value,rows in frame.groupby(column,dropna=True):
            if not str(value).strip() or str(value).lower() in {"nan","none","unknown"}:continue
            values=list(dict.fromkeys(mapping[x] for x in rows.sequence_sha256))
            for other in values[1:]:uf.union(values[0],other)
    out=frame.copy()
    out["group_id"]=[uf.find(mapping[x]) for x in out.sequence_sha256]
    unique=sorted(out.group_id.unique())
    if len(unique)<10: raise ValueError("At least 10 independent groups are needed for the frozen 70/10/10/10 split")
    rng=np.random.default_rng(seed);rng.shuffle(unique)
    n=len(unique);a=max(1,round(.7*n));b=max(a+1,round(.8*n));c=min(n-1,max(b+1,round(.9*n)))
    assignment={g:("train" if i<a else "validation" if i<b else "calibration" if i<c else "test") for i,g in enumerate(unique)}
    out["split"]=out.group_id.map(assignment)
    return out

def build_dataset(input,output,source="catapro",endpoint="kcat",cluster_file=None,mmseqs=None,seed=2026):
    output=Path(output)
    if output.exists(): raise FileExistsError(output)
    output.mkdir(parents=True)
    if endpoint not in {"kcat", "Km"}:
        raise ValueError("Select one endpoint: kcat or Km")
    raw=catapro_rows(input,endpoint) if source=="catapro" else pd.read_csv(input)
    if "endpoint" not in raw:
        raise ValueError("Measurements require an explicit endpoint column")
    other=raw.loc[raw.endpoint!=endpoint]
    raw_for_task=raw.loc[raw.endpoint==endpoint]
    clean,rejected=canonicalize(raw_for_task)
    excluded_other=pd.DataFrame([{"row":int(index)+2,"record_id":row.get("record_id"),
                                 "reason":"different_endpoint"} for index,row in other.iterrows()])
    rejected=pd.concat([rejected,excluded_other],ignore_index=True)
    rejected.to_csv(output/"excluded.csv",index=False)
    clean.to_csv(output/"measurements.csv",index=False)
    manifest={"schema_version":"1.0.0","source":source,"source_sha256":file_hash(input),
              "endpoint":endpoint,"raw_rows":len(raw),"accepted_rows":len(clean),"excluded_rows":len(rejected),
              "status":"cleaned_awaiting_homology_groups","split_seed":seed,
              "split_fractions":{"train":.7,"validation":.1,"calibration":.1,"test":.1},
              "conditions_policy":"unreported_conditions_are_not_imputed"}
    if cluster_file or mmseqs:
        if cluster_file:
            groups=pd.read_csv(cluster_file,sep="\t")
            provenance=Path(cluster_file).with_suffix(".json")
            if not provenance.exists(): raise ValueError("Imported groups require a matching .json provenance file")
            info=json.loads(provenance.read_text())
            if info.get("identity")!=.4 or info.get("coverage")!=.8:
                raise ValueError("Clustering thresholds do not match the frozen contract")
            if info.get("groups_sha256")!=file_hash(cluster_file):
                raise ValueError("Cluster mapping checksum mismatch")
        else: groups,info=mmseqs_groups(clean,output/"clustering",mmseqs)
        split=grouped_split(clean,groups,seed)
        split.to_csv(output/"split.csv",index=False)
        manifest.update(status="ready_for_features",groups=int(split.group_id.nunique()),clustering=info,
                        split_counts=split.split.value_counts().to_dict(),split_sha256=file_hash(output/"split.csv"))
    manifest["measurements_sha256"]=file_hash(output/"measurements.csv")
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    return {"output":str(output),**manifest}
