"""Package manifest dependencies while preserving the scientific manifest's bytes."""
from __future__ import annotations
import hashlib,json,posixpath
from pathlib import Path

def resolve_reference(value,base):
    base=Path(base).resolve()
    receipt=base/".aqua-transfer.json"
    if receipt.is_file():
        mapping=json.loads(receipt.read_text(encoding="utf-8"))
        if str(value) in mapping:
            # Receipts only point inside the enclosing task inputs directory.
            candidate=(base/mapping[str(value)]).resolve()
            task_inputs=base.parent
            if not candidate.is_relative_to(task_inputs):
                raise ValueError("Invalid transferred dependency")
            return candidate
    path=Path(value)
    return path.resolve() if path.is_absolute() else (base/path).resolve()

def manifest_dependencies(path):
    """Only explicit path+checksum entries are followed; directories are not swept."""
    path=Path(path)
    data=json.loads(path.read_text(encoding="utf-8"))
    entries=[]
    def visit(value):
        if isinstance(value,dict):
            if "path" in value and "sha256" in value:
                dependency=resolve_reference(value["path"],path.parent)
                if not dependency.is_file() or dependency.is_symlink():
                    raise ValueError("Missing manifest dependency: "+str(dependency))
                from .datasets import file_hash
                if file_hash(dependency)!=value["sha256"]:
                    raise ValueError("Manifest dependency checksum mismatch: "+str(dependency))
                entries.append((str(value["path"]),dependency))
            for child in value.values(): visit(child)
        elif isinstance(value,list):
            for child in value: visit(child)
    visit(data)
    return entries

def receipt_bytes(archive_name,dependencies):
    directory=posixpath.dirname(archive_name)
    return json.dumps({raw:posixpath.relpath(target,directory) for raw,target in dependencies},
                      ensure_ascii=False,sort_keys=True).encode("utf-8")
