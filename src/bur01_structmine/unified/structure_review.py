"""Full-construct coordinate review; structural eligibility is separate from MD release."""
from __future__ import annotations
import hashlib,json
from itertools import product
from pathlib import Path
import numpy as np
from .datasets import file_hash
from .mechanism import atomic_json
from .families import assess_structure

SIDECHAINS={'A':'CB','C':'CB SG','D':'CB CG OD1 OD2','E':'CB CG CD OE1 OE2','F':'CB CG CD1 CD2 CE1 CE2 CZ',
 'G':'','H':'CB CG ND1 CD2 CE1 NE2','I':'CB CG1 CG2 CD1','K':'CB CG CD CE NZ','L':'CB CG CD1 CD2',
 'M':'CB CG SD CE','N':'CB CG OD1 ND2','P':'CB CG CD','Q':'CB CG CD OE1 NE2','R':'CB CG CD NE CZ NH1 NH2',
 'S':'CB OG','T':'CB OG1 CG2','V':'CB CG1 CG2','W':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2','Y':'CB CG CD1 CD2 CE1 CE2 CZ OH'}


def close_pairs(coordinates,cutoff):
    cells={};offsets=list(product([-1,0,1],repeat=3))
    for index,xyz in enumerate(coordinates):
        key=tuple(np.floor(xyz/cutoff).astype(int))
        for offset in offsets:
            neighbor=tuple(a+b for a,b in zip(key,offset))
            for other in cells.get(neighbor,[]):
                if np.linalg.norm(xyz-coordinates[other])<cutoff:yield other,index
        cells.setdefault(key,[]).append(index)


def review_structure(structure,sequence,output,family_id='bur01',chain='A',scores=None,uncovered_n_terminal=0,atom_rules=None):
    from ..structure_search import load_ca_trace
    trace=load_ca_trace(structure,chain_id=chain)
    issues=[];missing=[];peptide=[];points=[];labels=[];resindices=[]
    if trace.sequence!=sequence:issues.append('full_construct_sequence_mismatch')
    for index,residue in enumerate(trace.residues):
        expected=('N CA C O '+SIDECHAINS.get(residue.amino_acid,'')).split()
        absent=sorted(set(expected)-set(residue.atoms))
        if absent:missing.append({'position':index+1,'residue':residue.label,'missing_atoms':absent})
        for atom,xyz in residue.atoms.items():
            if atom.startswith('H'):continue
            points.append(xyz);labels.append(f'{chain}:{residue.label}:{atom}');resindices.append(index)
        if index and 'N' in residue.atoms and 'C' in trace.residues[index-1].atoms:
            d=float(np.linalg.norm(residue.atoms['N']-trace.residues[index-1].atoms['C']))
            if not 1.1<=d<=1.6:peptide.append({'positions':[index,index+1],'distance_A':d})
    if missing:issues.append('missing_heavy_atoms')
    if peptide:issues.append('peptide_connectivity')
    xyz=np.array(points)
    clashes=[]
    if not np.isfinite(xyz).all():issues.append('nonfinite_coordinates')
    else:
        for a,b in close_pairs(xyz,1.8):
            if abs(resindices[a]-resindices[b])<=1:continue
            clashes.append({'atoms':[labels[a],labels[b]],'distance_A':float(np.linalg.norm(xyz[a]-xyz[b]))})
    if clashes:issues.append('severe_nonlocal_heavy_atom_clashes')
    quality=assess_structure(trace,family_id)
    if not quality['fold_pass']:issues.append('family_fold')
    if not quality['catalytic_identity_pass']:issues.append('catalytic_identity')
    if any(row['status']!='compatible' for row in quality['geometry']):issues.append('family_catalytic_geometry')
    geometry=[]
    for rule in atom_rules or []:
        a,b=rule['positions'];names_a=rule['first_atoms'];names_b=rule['second_atoms']
        try:
            first=np.array([trace.residues[a-1].atoms[n] for n in names_a]);second=np.array([trace.residues[b-1].atoms[n] for n in names_b])
            distance=float(np.linalg.norm(first[:,None]-second[None],axis=-1).min())
            passed=rule['minimum_A']<=distance<=rule['maximum_A']
        except (KeyError,IndexError):distance=None;passed=False
        geometry.append({**rule,'distance_A':distance,'passed':passed})
    if any(not g['passed'] for g in geometry):issues.append('protocol_specific_catalytic_geometry')
    confidence={}
    if scores:
        values=np.array(json.loads(Path(scores).read_text())['plddt'],float)
        if len(values)!=len(sequence) or not np.isfinite(values).all():raise ValueError('Confidence does not match construct')
        core=values[uncovered_n_terminal:]
        catalytic=[s['index'] for s in quality['sites'] if s['index'] is not None]
        confidence={'mean':float(values.mean()),'covered_core_mean':float(core.mean()),
                    'uncovered_terminal_mean':float(values[:uncovered_n_terminal].mean()) if uncovered_n_terminal else None,
                    'catalytic_minimum':float(values[catalytic].min()) if catalytic else None,
                    'below_50_positions':(np.flatnonzero(values<50)+1).tolist()}
        if core.mean()<70 or not catalytic or values[catalytic].min()<80:issues.append('low_core_or_catalytic_confidence')
    result={'version':'1.0.0','status':'passed' if not issues else 'blocked','structure_review_pass':not issues,
            'production_ready':False,'remaining_md_checks':['assembly','protonation','solution_parameters','equilibration'],
            'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),'structure_sha256':file_hash(structure),
            'family':quality,'confidence':confidence,'missing_atoms':missing,'peptide_breaks':peptide,
            'severe_clashes':clashes,'atom_specific_geometry':geometry,'blocking_issues':issues,
            'uncovered_n_terminal':uncovered_n_terminal,'interpretation':'Computational structural eligibility, not experimental validation.'}
    output=Path(output);output.mkdir(parents=True,exist_ok=True);atomic_json(output/'structure_review.json',result)
    return result


def build_reference_assembly(monomer,reference,sequence,output,*,source_url,chains=None):
    """Place a modeled chain into a specified public biological assembly by CA fitting."""
    from copy import deepcopy
    from Bio.PDB import PDBParser,PDBIO,Structure,Model
    from ..structure_search import load_ca_traces,load_ca_trace
    from ..prediction_sequence import _alignment_index_map
    refs=load_ca_traces(reference)
    if chains:refs=[r for r in refs if r.chain_id in chains]
    if not refs or len(refs)>26:raise ValueError('Select a supported reference assembly')
    model=load_ca_trace(monomer,chain_id='A')
    if model.sequence!=sequence:raise ValueError('Monomer differs from the confirmed construct')
    structure=Structure.Structure('assembly');target=Model.Model(0);structure.add(target)
    original=PDBParser(QUIET=True).get_structure('monomer',str(monomer))[0]['A']
    alignments=[]
    for index,ref in enumerate(refs):
        mapping,identity,_=_alignment_index_map(ref.sequence,sequence)
        if identity<.9 or len(mapping)<.8*len(ref.sequence):raise ValueError('Insufficient sequence agreement for assembly transfer')
        x=np.array([model.coords[j] for i,j in mapping.items()]);y=np.array([ref.coords[i] for i,j in mapping.items()])
        xc=x.mean(0);yc=y.mean(0);u,_,vt=np.linalg.svd((x-xc).T@(y-yc))
        correction=np.diag([1.,1.,np.linalg.det(u@vt)]);rotation=u@correction@vt;translation=yc-xc@rotation
        rmsd=float(np.sqrt(np.mean(np.sum((x@rotation+translation-y)**2,axis=1))))
        if rmsd>2.:raise ValueError('Reference assembly fit exceeds 2 A CA RMSD')
        chain=deepcopy(original);chain.id=chr(65+index);chain.transform(rotation,translation);target.add(chain)
        alignments.append({'chain':chain.id,'reference_chain':ref.chain_id,'aligned_residues':len(mapping),
                           'sequence_identity':identity,'rmsd_A':rmsd,'rotation':rotation.tolist(),'translation':translation.tolist()})
    atoms=[a for a in structure.get_atoms() if a.element!='H'];xyz=np.array([a.coord for a in atoms]);clashes=[];contacts=0
    for i,j in close_pairs(xyz,4.5):
        if atoms[i].get_parent().get_parent().id==atoms[j].get_parent().get_parent().id:continue
        contacts+=1;distance=float(np.linalg.norm(xyz[i]-xyz[j]))
        if distance<1.8:
            labels=[]
            for atom in (atoms[i],atoms[j]):
                residue=atom.get_parent()
                labels.append(f'{residue.get_parent().id}:{residue.id[1]}:{atom.name}')
            clashes.append({'atoms':labels,'distance_A':distance})
    output=Path(output);output.mkdir(parents=True,exist_ok=True);pdb=output/'assembly.pdb'
    writer=PDBIO();writer.set_structure(structure);writer.save(str(pdb))
    report={'version':'1.0.0','status':'passed' if not clashes and (len(refs)==1 or contacts>0) else 'blocked',
            'source_url':source_url,'reference_sha256':file_hash(reference),'monomer_sha256':file_hash(monomer),
            'assembly_sha256':file_hash(pdb),'chains':len(refs),'alignments':alignments,'interchain_contacts_4_5A':contacts,
            'severe_interchain_clashes':clashes,'construct_sequences':{r['chain']:sequence for r in alignments},
            'interpretation':'Homology transfer of a deposited biological assembly; candidate oligomerization is not experimentally established.'}
    atomic_json(output/'assembly_review.json',report);return report
