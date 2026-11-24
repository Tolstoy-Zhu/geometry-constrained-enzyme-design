"""SABIO-RK export API v1: freeze pages and map explicitly identified assays."""
from __future__ import annotations
import json,time
from pathlib import Path
import pandas as pd
import requests
from .datasets import file_hash
BASE="https://sabiork.h-its.org/export-api/sabio/kinlaw-entry"

def fetch(query,destination,page_size=100):
    if not query or not query.strip():raise ValueError("Provide an explicit SABIO-RK query")
    root=Path(destination)
    if root.exists():raise FileExistsError(root)
    root.mkdir(parents=True)
    pages=[];total=None
    for page in range(1,10001):
        params={"q":query,"page":page,"pageSize":page_size}
        response=requests.get(BASE+"/json",params=params,timeout=(20,120))
        response.raise_for_status();data=response.json()
        if not isinstance(data.get("data"),list) or "meta" not in data:raise ValueError("Unexpected SABIO export schema")
        meta=data["meta"]
        if total is None:total=meta["total_count"]
        if meta["total_count"]!=total:raise ValueError("Database changed during export; create a new snapshot")
        raw=root/f"page_{page:05d}.json";raw.write_bytes(response.content)
        xml=requests.get(BASE+"/sbml",params=params,timeout=(20,120));xml.raise_for_status()
        import xml.etree.ElementTree as ET
        if ET.fromstring(xml.content).tag.split("}")[-1]!="sbml":raise ValueError("Expected SBML")
        sbml=raw.with_suffix(".sbml");sbml.write_bytes(xml.content)
        pages.append({"path":raw.name,"sha256":file_hash(raw),"sbml":{"path":sbml.name,"sha256":file_hash(sbml)},"rows":len(data["data"])})
        if page>=meta["total_pages"]:break
        time.sleep(2.1)  # Both export formats share the official 60 requests/minute quota.
    if sum(p["rows"] for p in pages)!=total:raise ValueError("Incomplete SABIO export")
    manifest={"schema_version":"1.0.0","source":"SABIO-RK","api":BASE,"api_version":"1.0",
              "query":query,"total_rows":total,"pages":pages,"training_status":"requires_confirmed_construct_mapping",
              "terms":"https://sabiork.h-its.org/ui/about"}
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    return {"output":str(root),**manifest}

def map_assays(snapshot,mapping,output):
    """Mapping keys are entry IDs; mutant sequences are never replaced by WT UniProt sequences."""
    root=Path(snapshot);manifest=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    curated=json.loads(Path(mapping).read_text(encoding="utf-8"))
    rows=[];excluded=[]
    for page in manifest["pages"]:
        path=root/page["path"]
        if file_hash(path)!=page["sha256"]:raise ValueError("Frozen SABIO page checksum mismatch")
        for entry in json.loads(path.read_text(encoding="utf-8"))["data"]:
            entry_id=str(entry["id"]);enzyme=entry["enzyme_description"]
            defined=curated.get(entry_id,{})
            sequence=defined.get("sequence");source=defined.get("sequence_source")
            if not sequence or not source or not defined.get("construct_confirmed"):
                excluded.append({"entry_id":entry_id,"reason":"confirmed_construct_mapping_missing"});continue
            if enzyme.get("wildtype")=="Mutant" and not defined.get("mutations_confirmed"):
                excluded.append({"entry_id":entry_id,"reason":"mutant_sequence_not_confirmed"});continue
            substrate_species=[s["compound"] for s in entry["reaction"]["species"] if s["role"]=="Substrate"]
            compounds={int(x["id"]):x for x in substrate_species}
            smiles={}
            for link in entry.get("external_links",{}).get("compound",[]):
                if link["key"]=="SMILES":smiles.setdefault(int(link["id"]),set()).add(link["value"])
            conditions=entry.get("experimental_conditions") or {}
            ph=conditions.get("envvar_ph") or {};temp=conditions.get("envvar_temperature") or {}
            temperature=temp.get("start_value")
            if temp.get("end_value") is not None:temperature=None
            if temperature is not None:
                if temp.get("unit") in {"°C","C"}:temperature+=273.15
                elif temp.get("unit")!="K":temperature=None
            for i,parameter in enumerate(entry["kineticlaw"]["parameter"]):
                endpoint=(parameter.get("parameter_type") or {}).get("name")
                if endpoint not in {"kcat","Km"}:continue
                identifier=f"sabio_{entry_id}_{i}"
                assigned=defined.get("parameter_substrates",{}).get(parameter["name"])
                if assigned is None:
                    key=(parameter.get("species") or {}).get("species_key") or ""
                    parts=[s.strip() for s in key.split("|")]
                    found=[c["id"] for c in substrate_species if len(parts)>=2 and c["name"]==parts[1]]
                    if len(found)==1:assigned=found[0]
                    elif len(compounds)==1:assigned=next(iter(compounds))
                if assigned is None or int(assigned) not in compounds or len(smiles.get(int(assigned),set()))!=1:
                    excluded.append({"entry_id":entry_id,"record_id":identifier,"reason":"parameter_substrate_ambiguous"});continue
                if parameter.get("n_end_value") is not None:
                    excluded.append({"entry_id":entry_id,"record_id":identifier,"reason":"interval_value_requires_curated_qualifier"});continue
                unit=(parameter.get("unit") or {}).get("n_name")
                unit={"s^(-1)":"s^-1","min^(-1)":"min^-1"}.get(unit,unit)
                rows.append({"record_id":identifier,"sequence":sequence,"smiles":next(iter(smiles[int(assigned)])),
                    "endpoint":endpoint,"value":parameter.get("n_start_value"),"unit":unit,"data_kind":"experimental",
                    "qualifier":"=","source":BASE+"/json/"+entry_id,"source_entry":entry_id,
                    "publication":json.dumps(entry.get("publication"),ensure_ascii=False),
                    "parent_id":defined.get("parent_id",entry_id),"measurement_group":"sabio_entry_"+entry_id,
                    "method":entry["general"].get("experiment_type"),"replicate_id":"not_reported",
                    "temperature_K":temperature,"pH":ph.get("start_value") if ph.get("end_value") is None else None,
                    "buffer":conditions.get("buffer"),"sequence_source":source,
                    "enzyme_type":enzyme.get("wildtype"),"mutant_description":enzyme.get("mutant_spec")})
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    columns=["record_id","sequence","smiles","endpoint","value","unit","data_kind","qualifier","source"]
    pd.DataFrame(rows,columns=None if rows else columns).to_csv(output/"mapped.csv",index=False)
    pd.DataFrame(excluded,columns=None if excluded else ["entry_id","record_id","reason"]).to_csv(output/"unmapped.csv",index=False)
    result={"rows":len(rows),"unmapped":len(excluded),"snapshot_sha256":file_hash(root/"manifest.json"),
            "mapping_sha256":file_hash(mapping),"mapped_sha256":file_hash(output/"mapped.csv")}
    (output/"mapping_report.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    return {"output":str(output),**result}

