"""Standalone frozen ESM-2 extractor, usable in the existing property runtime."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

def extract(input,output,model,device="cpu",max_length=1022):
    import torch
    from transformers import AutoTokenizer,AutoModel
    rows=json.loads(Path(input).read_text(encoding="utf-8"))
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    model=Path(model)
    h=hashlib.sha256()
    with (model/'model.safetensors').open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    identity={'model_sha256':h.hexdigest(),'max_length':max_length,'method':'ESM2_t33_mean_residue_embedding'}
    receipt=out/'encoder.json'
    if receipt.exists():
        if json.loads(receipt.read_text())!=identity:raise ValueError('Embedding cache belongs to a different encoder')
    elif any(out.glob('*.npy')):
        raise ValueError('Existing embeddings have no verified encoder receipt')
    receipt.write_text(json.dumps(identity),encoding='utf-8')
    tokenizer=AutoTokenizer.from_pretrained(model,local_files_only=True)
    encoder=AutoModel.from_pretrained(model,local_files_only=True,add_pooling_layer=False).to(device).eval()
    torch.set_num_threads(min(6,__import__("os").cpu_count() or 1))
    manifest={"method":"ESM2_t33_mean_residue_embedding","device":device,"dimension":encoder.config.hidden_size,
              "max_length":max_length,"excluded":[],"files":{}}
    pending=[]
    for key,sequence in rows.items():
        if len(sequence)>max_length:
            manifest["excluded"].append({"sequence_sha256":key,"length":len(sequence),"reason":"above_model_limit_no_truncation"})
            continue
        if hashlib.sha256(sequence.encode()).hexdigest()!=key: raise ValueError("Sequence hash mismatch")
        path=out/(key+".npy")
        if path.exists():
            vector=np.load(path,allow_pickle=False)
            if vector.shape!=(encoder.config.hidden_size,) or not np.isfinite(vector).all():raise ValueError("Invalid embedding")
            manifest["files"][key]=hashlib.sha256(path.read_bytes()).hexdigest()
        else:pending.append((key,sequence))
    pending.sort(key=lambda row:len(row[1]))
    while pending:
        size=min(16 if device=="cuda" else 1,len(pending))
        while size>1 and max(len(s)+2 for _,s in pending[:size])*size>8192:size//=2
        batch=pending[:size]
        tokens=tokenizer([s for _,s in batch],return_tensors="pt",padding=True,truncation=False).to(device)
        with torch.inference_mode():
            hidden=encoder(**tokens).last_hidden_state
            vectors=[hidden[i,1:len(s)+1].mean(0).float().cpu().numpy() for i,(_,s) in enumerate(batch)]
        for (key,_),vector in zip(batch,vectors):
            if not np.isfinite(vector).all():raise ValueError("Invalid embedding")
            path=out/(key+".npy");temp=path.with_suffix(".npy.tmp")
            with temp.open("wb") as handle:np.save(handle,vector,allow_pickle=False)
            temp.replace(path)
            manifest["files"][key]=hashlib.sha256(path.read_bytes()).hexdigest()
        del hidden,tokens,vectors
        pending=pending[size:]
        print(f"embeddings_saved={len(manifest['files'])} remaining={len(pending)}",flush=True)
    manifest["model_sha256"]=identity['model_sha256']
    temp=out/"manifest.json.tmp"
    temp.write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    temp.replace(out/"manifest.json")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("input");p.add_argument("output");p.add_argument("model")
    p.add_argument("--device",default="cpu",choices=["cpu","cuda"])
    args=p.parse_args();extract(**vars(args))
