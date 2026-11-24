"""Generate public reference-chain records; no experimental measurements are fabricated."""
import argparse,json
from pathlib import Path
from bur01_structmine.unified.families import profiles,reference,assess_structure
from bur01_structmine.unified.contracts import EnzymeRecord
from bur01_structmine.structure_search import load_ca_trace

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    for family in profiles():
        path=reference(family);trace=load_ca_trace(path,chain_id=family.reference_chain)
        record=EnzymeRecord(enzyme_id=family.family_id+"_reference_chain",family_id=family.family_id,
            sequence=trace.sequence,name=family.name+" public reference",
            source="public_reference",construct_id=family.family_id+"_resolved_coordinates",
            construct_description="Sequence resolved in the public structure; not a confirmed expression construct",
            chain=family.reference_chain,provenance={"reference_sha256":family.reference_sha256,"sources":family.sources})
        payload={"enzyme":[record.model_dump(mode="json")],"reaction":[],"assay":[]}
        (args.output/(family.family_id+"_records.json")).write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding="utf-8")
        (args.output/(family.family_id+"_reference.fasta")).write_text(">"+record.enzyme_id+"\n"+record.sequence+"\n",encoding="ascii")
        result=assess_structure(trace,family.family_id)
        (args.output/(family.family_id+"_structure_check.json")).write_text(json.dumps(result,indent=2),encoding="utf-8")
        print(family.family_id,len(trace.sequence),result["catalytic_identity_pass"])
if __name__=="__main__":main()
