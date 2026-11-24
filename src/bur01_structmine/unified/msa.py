"""Local monomer MSA reuse, with explicit unmapped terminal columns."""
from __future__ import annotations
import json,re
from pathlib import Path
from .datasets import file_hash


def read_monomer_a3m(path):
    records=[];name=None;parts=[]
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('#'):
            fields=line[1:].split()
            if len(fields)!=2 or ',' in fields[0] or fields[1]!='1':
                raise ValueError('Only a single-chain monomer MSA is supported')
            continue
        if line.startswith('>'):
            if name is not None:records.append((name,''.join(parts)))
            name=line[1:];parts=[]
        else:
            if name is None:raise ValueError('A3M sequence has no FASTA header')
            if not re.fullmatch('[A-Za-z-]+',line):raise ValueError('Invalid A3M residue or gap')
            parts.append(line)
    if name is not None:records.append((name,''.join(parts)))
    if not records or not records[0][1]:raise ValueError('MSA is empty')
    query=records[0][1]
    if not query.isupper() or '-' in query:raise ValueError('MSA query must be an ungapped uppercase sequence')
    if any(len(re.sub('[a-z]','',seq))!=len(query) for _,seq in records):
        raise ValueError('MSA aligned column count differs from its query')
    return records


def extend_msa(msa,reference_core,full_query,offset,output,query_name='query'):
    if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',full_query):
        raise ValueError('Full construct must contain canonical uppercase amino acids')
    if isinstance(offset,bool) or not isinstance(offset,int):raise ValueError('Core offset must be an integer')
    if not query_name or '\n' in query_name or '\r' in query_name:raise ValueError('Invalid query name')
    rows=read_monomer_a3m(msa)
    if rows[0][1]!=reference_core:raise ValueError('Cached MSA does not match the original core sequence')
    if offset<0 or offset+len(reference_core)>len(full_query):raise ValueError('Core mapping lies outside the full construct')
    suffix=len(full_query)-offset-len(reference_core)
    changes=[{'position':offset+i+1,'original':a,'query':b} for i,(a,b) in enumerate(zip(reference_core,full_query[offset:])) if a!=b]
    lines=[f'#{len(full_query)}\t1','>'+query_name,full_query]
    for name,seq in rows[1:]:lines.extend(['>'+name,'-'*offset+seq+'-'*suffix])
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    read_monomer_a3m(output)
    receipt={'version':'1.0.0','source_sha256':file_hash(msa),'output_sha256':file_hash(output),
             'alignment_rows':len(rows),'core_length':len(reference_core),'full_length':len(full_query),
             'uncovered_terminal_columns':{'n_terminal':offset,'c_terminal':suffix},'query_changes':changes,
             'homolog_rows_reused':True,'new_public_sequence_submission':False,
             'interpretation':'Terminal columns are gaps in homologs, not newly searched evolutionary evidence.'}
    output.with_suffix('.provenance.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    return receipt
