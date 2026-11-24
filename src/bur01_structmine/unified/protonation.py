"""Local pKa estimates and explicit fixed-charge assumptions for protein preparation."""
from __future__ import annotations
import contextlib,json,math,os
from pathlib import Path
from importlib.metadata import version
from .datasets import file_hash
from .mechanism import atomic_json


def normalize_protein_pdb(source,destination):
    """Canonical PDB atom-name alignment prevents CD/CE/ND being read as metals."""
    lines=[]
    variants={'HID':'HIS','HIE':'HIS','HIP':'HIS','ASH':'ASP','GLH':'GLU','LYN':'LYS','CYM':'CYS'}
    for line in Path(source).read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            name=line[12:16].strip();element=line[76:78].strip().upper() or name.lstrip('0123456789')[:1]
            if element in {'H','D'}:continue
            if element not in {'C','N','O','S'}:raise ValueError('Unsupported protein element: '+element)
            atom_name=' '+name.ljust(3) if len(name)<4 else name
            line=line[:12]+atom_name+line[16:]
            line=line[:17]+variants.get(line[17:20],line[17:20])+line[20:]
        lines.append(line)
    Path(destination).write_text('\n'.join(lines)+'\n',encoding='ascii')


def estimate_pka(structure,output,pH=4.):
    from propka.run import single
    if version('propka')!='3.5.1':raise ValueError('Use PROPKA 3.5.1 for this frozen protocol')
    if not math.isfinite(pH) or not 0<=pH<=14:raise ValueError('pH must lie between zero and fourteen')
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    identity={'input_sha256':file_hash(structure),'pH':pH,'propka':'3.5.1'}
    old=output/'pka_report.json'
    if old.exists():
        previous=json.loads(old.read_text())
        if previous['identity']!=identity:raise ValueError('pKa estimate inputs changed')
        if file_hash(output/'normalized_input.pdb')!=previous['normalized_input_sha256']:raise ValueError('Normalized pKa input changed')
        return previous
    normalized=output/'normalized_input.pdb';normalize_protein_pdb(structure,normalized)
    prior=Path.cwd()
    try:
        os.chdir(output)
        with (output/'propka.log').open('w') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
            molecule=single(str(normalized),optargs=['--pH',str(pH)],write_pka=True)
    finally:os.chdir(prior)
    groups=[]
    for group in molecule.conformations['AVR'].groups:
        if not group.titratable:continue
        atom=group.atom
        groups.append({'label':group.label,'type':group.type,'residue':atom.res_name,'position':atom.res_num,
            'chain':atom.chain_id,'pka':group.pka_value,
            'coupled':[g.label for g in group.non_covalently_coupled_groups]})
    if not groups or any(not math.isfinite(r['pka']) for r in groups):raise ValueError('Invalid pKa estimates')
    result={'schema_version':'1.0.0','identity':identity,'status':'estimated','groups':groups,
        'normalized_input_sha256':file_hash(normalized),'production_ready':False,
        'interpretation':'Structure-based pKa estimates; fixed protonation and histidine tautomers require a recorded policy.'}
    atomic_json(old,result);return result


def assign_microstate(report, *, pH, overrides=(), reference=None):
    """Use symmetry-averaged estimates; conserved sites can inherit a parent policy.

    This is a fixed, conditional microstate, not a constant-pH population estimate.
    Overrides must identify the expected residue, selected Amber variant and reason.
    """
    from collections import defaultdict
    import statistics
    if not math.isfinite(pH) or not 0<=pH<=14:raise ValueError('Invalid pH')
    if reference and reference['pH']!=pH:raise ValueError('Parent reference pH differs')
    if len({int(r['position']) for r in overrides})!=len(overrides):raise ValueError('Duplicate protonation override')
    grouped=defaultdict(list)
    for group in report['groups']:
        # PROPKA exposes the C terminus as a COO group on the terminal amino acid.
        kind='C-' if group['label'].startswith('C-') else group['type']
        grouped[(group['position'],group['residue'],kind)].append(group)
    inherited={(r['position'],r['residue'],r['type']):r for r in (reference or {}).get('assignments',[])}
    overrides={int(r['position']):r for r in overrides}
    allowed={'ASP':{'ASP','ASH'},'GLU':{'GLU','GLH'},'HIS':{'HID','HIE','HIP'},
             'CYS':{'CYS','CYM'},'LYS':{'LYS','LYN'},'TYR':{'TYR'},'ARG':{'ARG'}}
    assignments=[];used=set()
    for (position,residue,kind),groups in sorted(grouped.items()):
        pka=statistics.mean(g['pka'] for g in groups)
        if not math.isfinite(pka):raise ValueError('Nonfinite pKa')
        protonated=pka>=pH
        if kind=='N+':
            if not protonated:raise ValueError('Neutral N terminus requires an explicit template')
            selected='charged_N_terminus'
        elif kind=='C-':
            if protonated:raise ValueError('Neutral C terminus requires an explicit template')
            selected='charged_C_terminus'
        elif residue=='ASP':selected='ASH' if protonated else 'ASP'
        elif residue=='GLU':selected='GLH' if protonated else 'GLU'
        elif residue=='HIS':selected='HIP' if protonated else 'HIE'
        elif residue=='CYS':selected='CYS' if protonated else 'CYM'
        elif residue=='LYS':selected='LYS' if protonated else 'LYN'
        elif residue in {'ARG','TYR'}:
            if not protonated:raise ValueError('Unsupported neutral arginine or deprotonated tyrosine')
            selected=residue
        else:raise ValueError('Unsupported titratable group '+str((residue,kind)))
        reason='Symmetry-averaged structure pKa relative to requested pH'
        parent=inherited.get((position,residue,kind))
        if parent:
            selected=parent['variant'];reason='Conserved-site state held at the parent reference policy'
        override=overrides.get(position) if kind not in {'N+','C-'} else None
        if override:
            if override['residue']!=residue or override['variant'] not in allowed.get(residue,set()) or not override.get('reason'):
                raise ValueError('Invalid explicit protonation override')
            selected=override['variant'];reason=override['reason'];used.add(position)
        assignments.append({'position':position,'residue':residue,'type':kind,'variant':selected,'reason':reason,
            'mean_pka':pka,'chain_pka':{g['chain']:g['pka'] for g in groups},
            'near_pH':abs(pka-pH)<=.5,'coupled_groups':sorted({x for g in groups for x in g['coupled']})})
    if used!=set(overrides):raise ValueError('Protonation override refers to an absent titratable residue')
    return {'schema_version':'1.0.0','pH':pH,'assignments':assignments,'production_ready':False,
            'assumptions':['Fixed protonation; structure pKa estimates do not measure microstate populations.',
                'Equivalent sites in identical chains use the mean pKa; conserved variant sites inherit the specified parent.',
                'Neutral non-overridden histidines use HIE; neutral acids use Amber ASH/GLH templates.',
                'Standard charged termini; unsupported terminal states stop preparation.']}


def apply_microstate(source,destination,policy):
    mapping={(r['position'],r['residue']):r['variant'] for r in policy['assignments'] if r['type'] not in {'N+','C-'}}
    lines=[]
    for line in Path(source).read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            name=mapping.get((int(line[22:26]),line[17:20]),line[17:20])
            line=line[:17]+name+line[20:]
        lines.append(line)
    Path(destination).write_text('\n'.join(lines)+'\n',encoding='ascii')
