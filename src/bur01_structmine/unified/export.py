"""Small result transfer, including explicitly referenced structure coordinates."""
from __future__ import annotations
import csv,hashlib,io,json,tarfile
from pathlib import Path
EXTENSIONS={".json",".jsonl",".csv",".tsv",".txt",".html",".pdb",".cif",".fasta",".xlsx",".log",".yaml",".zip"}
COORDINATES={".pdb",".cif",".fasta"}
def export_task(store,task_id,stream):
    base=store.path.parent/"tasks"/task_id/"output"
    deployment=store.path.parent.parent.resolve()
    files={};mapping={}
    for path in sorted(base.rglob("*")):
        if path.is_file() and not path.is_symlink() and path.suffix.lower() in EXTENSIONS:
            name=path.relative_to(base).as_posix()
            files[name]=path;mapping[str(path.resolve())]=name
    def add(value):
        if isinstance(value,dict):
            for child in value.values():add(child)
        elif isinstance(value,list):
            for child in value:add(child)
        elif isinstance(value,str) and Path(value).suffix.lower() in COORDINATES:
            path=Path(value).resolve()
            if path.is_relative_to(deployment) and path.is_file() and not path.is_symlink():
                name=mapping.get(str(path))
                if not name:
                    name="referenced/"+hashlib.sha256(str(path).encode()).hexdigest()[:16]+"/"+path.name
                    files[name]=path;mapping[str(path)]=name
    for path in list(files.values()):
        try:
            if path.suffix==".json":add(json.loads(path.read_text(encoding="utf-8")))
            elif path.suffix in {".csv",".tsv"}:
                with path.open(encoding="utf-8-sig",newline="") as handle:
                    for row in csv.DictReader(handle,delimiter="\t" if path.suffix==".tsv" else ","):
                        add({k:v for k,v in row.items() if k.endswith(("_path","_file"))})
        except (OSError,ValueError):continue
    with tarfile.open(fileobj=stream,mode="w|gz") as tar:
        if sum(p.stat().st_size for p in files.values())>1_900_000_000:
            raise ValueError("Result tables exceed the automatic download limit; fetch selected files separately")
        data=json.dumps(mapping,ensure_ascii=False).encode("utf-8")
        entry=tarfile.TarInfo(".aqua-download.json");entry.size=len(data)
        tar.addfile(entry,io.BytesIO(data))
        for name,path in files.items():tar.add(path,arcname=name,recursive=False)
