"""Explicit ff19SB/OPC/GAFF2 free-enzyme systems with counted cosolvents."""
from __future__ import annotations
import json,math,os,re,subprocess,shutil,uuid
from pathlib import Path
import numpy as np
from .datasets import file_hash
from .mechanism import atomic_json
from .ligand_parameters import _run,normalize_rounded_charge
from .protonation import estimate_pka,assign_microstate,apply_microstate

CANONICAL={'HID':'HIS','HIE':'HIS','HIP':'HIS','ASH':'ASP','GLH':'GLU','LYN':'LYS','CYM':'CYS'}


def protein_lines(topology,xyz_A,sequences,*,hydrogens=False,canonical=True):
    """Retain explicit chain identities instead of Amber's merged residue numbering."""
    order=list(sequences);lengths=[len(sequences[c]) for c in order]
    ends=np.cumsum(lengths);lines=[];serial=1;previous=None
    for residue in topology.residues():
        if residue.index>=ends[-1]:break
        which=int(np.searchsorted(ends,residue.index,side='right'))
        chain=order[which];number=residue.index-(0 if which==0 else ends[which-1])+1
        if previous is not None and previous!=chain:lines.append('TER')
        previous=chain;name=CANONICAL.get(residue.name,residue.name) if canonical else residue.name
        for atom in residue.atoms():
            if atom.element is None or (not hydrogens and atom.element.symbol=='H'):continue
            x,y,z=xyz_A[atom.index];field=' '+atom.name.ljust(3) if len(atom.name)<4 else atom.name
            lines.append(f'ATOM  {serial:5d} {field} {name:3s} {chain}{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom.element.symbol:>2s}')
            serial+=1
    return lines+['TER','END']


def verify_construct(topology,sequences):
    from Bio.Data.PDBData import protein_letters_3to1
    expected=''.join(sequences.values());residues=list(topology.residues())[:len(expected)]
    if ''.join(protein_letters_3to1.get(CANONICAL.get(r.name,r.name),'?') for r in residues)!=expected:
        raise ValueError('Amber protein sequence differs from confirmed constructs')
    boundaries=set(np.cumsum([len(s) for s in sequences.values()])[:-1])
    for a,b in topology.bonds():
        if any(min(a.residue.index,b.residue.index)<i<=max(a.residue.index,b.residue.index) for i in boundaries):
            raise ValueError('Unexpected covalent bond between protein chains')


def component_counts(components,volume_nm3):
    if not math.isfinite(volume_nm3) or volume_nm3<=0:raise ValueError('Invalid box volume')
    names=[r['name'] for r in components]
    if len(names)!=len(set(names)):raise ValueError('Duplicate solution component')
    result={}
    for item in components:
        concentration=float(item['concentration_M'])
        if not math.isfinite(concentration) or not 0<concentration<=5:raise ValueError('Invalid solution concentration')
        number=int(round(concentration*.602214076*volume_nm3))
        if number<1:raise ValueError('Box cannot represent requested component; review integer composition')
        result[item['name']]=number
    return result


def prepare_protein_system(input,structure,parameter_dir,assembly_report,output,protonation_reference=None):
    from openmm import app,unit,NonbondedForce
    request=json.loads(Path(input).read_text());output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    review=json.loads(Path(assembly_report).read_text())
    if review['status']!='passed' or review.get('blocking_issues') or review['output_sha256']!=file_hash(structure):
        raise ValueError('Assembly has not passed coordinate checks for this file')
    if request.get('enzyme_state')!='free':raise ValueError('Only free enzyme is supported')
    for key in ['unvalidated_metals','covalent_ligands','nonstandard_residues']:
        if request.get(key):raise ValueError('Unsupported preparation item: '+key)
    sequences=request['construct_sequences']
    if len(set(sequences.values()))!=1:raise ValueError('This protonation symmetry policy supports identical chains only')
    amber=Path(os.environ['AQUA_AMBERHOME']);parameter_dir=Path(parameter_dir).resolve()
    env=os.environ.copy();env.update(AMBERHOME=str(amber),LD_LIBRARY_PATH=str(amber/'lib'),OMP_NUM_THREADS='2')
    env['PATH']=str(amber/'bin')+os.pathsep+env.get('PATH','')
    identity={'request':request,'structure_sha256':file_hash(structure),'assembly_sha256':file_hash(assembly_report),
              'parameter_reports':{},'reference_sha256':file_hash(protonation_reference) if protonation_reference else None}
    reports={}
    for item in request['components']:
        name=item['name'];folder=parameter_dir/name
        if not re.fullmatch('[A-Za-z0-9_]+',name) or folder.parent!=parameter_dir:raise ValueError('Invalid component name')
        report=json.loads((folder/'parameter_report.json').read_text())
        if report['status']!='numerical_checks_passed':raise ValueError('Parameters have not passed numerical checks')
        for name2,digest in report['sha256'].items():
            if file_hash(folder/name2)!=digest:raise ValueError('Parameter files changed')
        identity['parameter_reports'][name]=file_hash(folder/'parameter_report.json');reports[name]=report
    old=output/'identity.json'
    if old.exists() and json.loads(old.read_text())!=identity:raise ValueError('Protein system inputs changed')
    if (output/'system_build.json').is_file():
        result=json.loads((output/'system_build.json').read_text())
        for name,digest in result['sha256'].items():
            if file_hash(output/name)!=digest:raise ValueError('Built system changed')
        return {'output':str(output),**result}
    if old.exists():
        history=output/'history'/uuid.uuid4().hex;history.mkdir(parents=True)
        for p in output.iterdir():
            if p.is_file():shutil.copy2(p,history/p.name)
    atomic_json(old,identity)
    pka=estimate_pka(structure,output/'pka',request['pH'])
    reference=json.loads(Path(protonation_reference).read_text()) if protonation_reference else None
    policy=assign_microstate(pka,pH=request['pH'],overrides=request.get('protonation_overrides',[]),reference=reference)
    atomic_json(output/'protonation.json',policy)
    apply_microstate(output/'pka/normalized_input.pdb',output/'microstate.pdb',policy)
    commands=['source leaprc.protein.ff19SB','source leaprc.water.opc','source leaprc.gaff2']
    (output/'protein.leap').write_text('\n'.join(commands+['protein = loadpdb microstate.pdb','check protein',
        'saveamberparm protein protein.prmtop protein.inpcrd','quit'])+'\n')
    print('Building protonated protein with AmberTools',flush=True)
    _run([amber/'bin/tleap','-f','protein.leap'],output,env,'protein')
    protein=app.AmberPrmtopFile(str(output/'protein.prmtop'));coords=app.AmberInpcrdFile(str(output/'protein.inpcrd'))
    verify_construct(protein.topology,sequences)
    if protein.topology.getNumResidues()!=sum(map(len,sequences.values())):raise ValueError('Unexpected protein residues')
    protein_system=protein.createSystem(nonbondedMethod=app.NoCutoff)
    protein_mass=sum(protein_system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(protein_system.getNumParticles()))
    nb=next(f for f in protein_system.getForces() if isinstance(f,NonbondedForce))
    protein_charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(nb.getNumParticles()))
    if abs(protein_charge-round(protein_charge))>1e-4:raise ValueError('Nonintegral protein charge')
    xyz=np.array(coords.positions.value_in_unit(unit.angstrom))
    padding=float(request.get('padding_nm',1.2))
    if not 1.0<=padding<=2.0:raise ValueError('Use a reviewed box padding between 1 and 2 nm')
    edges=np.ceil((np.ptp(xyz,axis=0)/10+2*padding)*10)/10
    edges=np.maximum(edges,6.)
    if np.prod(edges)>1500:raise ValueError('Review systems larger than 1500 nm^3 separately')
    xyz+=edges*5-(xyz.max(axis=0)+xyz.min(axis=0))/2
    (output/'protein_centered.pdb').write_text('\n'.join(protein_lines(protein.topology,xyz,sequences,hydrogens=True,canonical=False))+'\n')
    # OpenMM canonicalizes Amber residue labels (HIP -> HIS, ASH -> ASP).
    # Restore the chosen templates before a second LEaP load, retaining their hydrogens.
    apply_microstate(output/'protein_centered.pdb',output/'protein_centered.pdb',policy)
    counts=component_counts(request['components'],float(np.prod(edges)));mass_solutes=0.;net_charge=int(round(protein_charge))
    templates={};rounding={};units=commands+['savepdb OPCBOX opc_box.pdb'];reload=list(commands)
    for index,(name,count) in enumerate(counts.items()):
        if index>=26:raise ValueError('Too many solute species')
        residue='QA'+chr(65+index);folder=parameter_dir/name;templates[residue]={'name':name,'count':count}
        mol2=re.sub(r'\bLIG\b',residue,(folder/'ligand.mol2').read_text());(output/(residue+'.mol2')).write_text(mol2)
        # Legacy individually accepted 1e-6 e residuals accumulate over many copies.
        # Normalize a derived template only; retain the verified source and raw copy.
        rounding[name]=normalize_rounded_charge(output/(residue+'.mol2'),reports[name]['formal_charge'])
        load=[f'loadamberparams "{folder / "ligand.frcmod"}"',f'{residue} = loadmol2 {residue}.mol2']
        units+=load+[f'savepdb {residue} {residue}.pdb'];reload+=load
        ligand=app.AmberPrmtopFile(str(folder/'ligand.prmtop')).createSystem(nonbondedMethod=app.NoCutoff)
        mass_solutes+=count*sum(ligand.getParticleMass(i).value_in_unit(unit.dalton) for i in range(ligand.getNumParticles()))
        net_charge+=count*reports[name]['formal_charge']
    (output/'units.leap').write_text('\n'.join(units+['quit'])+'\n');_run([amber/'bin/tleap','-f','units.leap'],output,env,'units')
    water=[];first=None
    for line in (output/'opc_box.pdb').read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')):continue
        if first is None:first=line[21:27]
        if line[21:27]!=first:break
        water.append(line)
    if len(water)!=4:raise ValueError('Expected four-site OPC')
    (output/'water.pdb').write_text('\n'.join(water+['TER','END'])+'\n')
    # The partial specific volume is only an initial packing estimate. NPT volume is measured later.
    protein_volume=protein_mass*.73/602.214076
    initial_density=float(request.get('initial_solvent_density_g_ml',.986))
    water_count=round((initial_density*602.214076*(np.prod(edges)-protein_volume)-mass_solutes-abs(net_charge)*35.45)/18.01528)
    if water_count<100:raise ValueError('Insufficient solvent')
    pack=['tolerance 2.1','filetype pdb','output packed.pdb','seed '+str(request['seed']),'nloop 1000',
          'structure protein_centered.pdb','number 1','fixed 0. 0. 0. 0. 0. 0.','end structure']
    for residue,item in templates.items():
        pack += [f'structure {residue}.pdb',f'number {item["count"]}',
                 'inside box 1.5 1.5 1.5 '+' '.join(str(v*10-1.5) for v in edges),'end structure']
    pack+=['structure water.pdb',f'number {water_count}','inside box 1.1 1.1 1.1 '+' '.join(str(v*10-1.1) for v in edges),'end structure']
    (output/'packmol.inp').write_text('\n'.join(pack)+'\n')
    print(f'Packing {water_count} waters, counted solutes, box {edges.tolist()} nm',flush=True)
    packing_inputs={name:file_hash(output/name) for name in ['packmol.inp','protein_centered.pdb','water.pdb',*[r+'.pdb' for r in templates]]}
    cache_path=output/'packing_cache.json';cache=json.loads(cache_path.read_text()) if cache_path.exists() else {}
    if cache.get('inputs')!=packing_inputs or not (output/'packed.pdb').is_file() or cache.get('packed_sha256')!=file_hash(output/'packed.pdb'):
        with (output/'packmol.inp').open('rb') as source,(output/'packmol.log').open('wb') as log:
            subprocess.run([str(amber/'bin/packmol')],stdin=source,stdout=log,stderr=subprocess.STDOUT,cwd=output,env=env,check=True,timeout=1800)
        atomic_json(cache_path,{'inputs':packing_inputs,'packed_sha256':file_hash(output/'packed.pdb')})
    else:print('Reusing checksum-verified completed packing',flush=True)
    raw=[r for r in (output/'packed.pdb').read_text().splitlines() if r.startswith(('ATOM  ','HETATM'))]
    protein_atoms=protein.topology.getNumAtoms();coordinates=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in raw[:protein_atoms]])
    if coordinates.shape!=xyz.shape or np.max(np.abs(coordinates-xyz))>.002:raise ValueError('Packing changed fixed protein coordinates')
    # Restore chain breaks from the original protein. Each solvent molecule has its own TER.
    # LEaP rebuilds protein hydrogens after packing. OpenMM also normalizes
    # terminal H1 -> H, which is not a valid input atom for an Amber N-terminal template.
    packed=[]
    for line in (output/'protein_centered.pdb').read_text().splitlines()[:-1]:
        if line.startswith(('ATOM  ','HETATM')) and line[76:78].strip()=='H':continue
        packed.append(line)
    last=None
    for line in raw[protein_atoms:]:
        key=line[17:27]
        if last is not None and key!=last:packed.append('TER')
        packed.append(line);last=key
    (output/'packed_terminated.pdb').write_text('\n'.join(packed+['TER','END'])+'\n')
    # Fresh templates: savepdb before loadpdb can mutate AmberTools 26 ligand units.
    build=reload+['system = loadpdb packed_terminated.pdb','set system box { '+' '.join(str(v*10) for v in edges)+' }']
    if net_charge:build+=['addionsrand system '+('Na+' if net_charge<0 else 'Cl-')+' 0']
    build+=['check system','saveamberparm system system.prmtop system.inpcrd','quit']
    (output/'system.leap').write_text('\n'.join(build)+'\n');_run([amber/'bin/tleap','-f','system.leap'],output,env,'system')
    topology=app.AmberPrmtopFile(str(output/'system.prmtop'));positions=app.AmberInpcrdFile(str(output/'system.inpcrd'))
    verify_construct(topology.topology,sequences)
    if topology.topology.getNumAtoms()!=len(positions.positions):raise ValueError('Topology coordinate count mismatch')
    system=topology.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1*unit.nanometer,constraints=app.HBonds,rigidWater=True)
    charge=sum(f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for f in system.getForces() if isinstance(f,NonbondedForce) for i in range(f.getNumParticles()))
    if abs(charge)>1e-4:raise ValueError('System is not neutral')
    from collections import Counter
    observed=Counter(r.name for r in topology.topology.residues())
    for residue,item in templates.items():
        if observed[residue]!=item['count']:raise ValueError('Solute count changed during preparation')
    with (output/'system_protein.pdb').open('w') as handle:
        handle.write('\n'.join(protein_lines(topology.topology,np.array(positions.positions.value_in_unit(unit.angstrom)),sequences))+'\n')
    report={'schema_version':'1.0.0','status':'built_awaiting_equilibration','production_ready':False,
        'atoms':system.getNumParticles(),'protein_residues':sum(map(len,sequences.values())),
        'box_nm':edges.tolist(),'protein_charge':protein_charge,'system_charge':charge,'component_counts':templates,
        'water_molecules':observed['WAT']+observed['HOH'],
        'counterions':{'Na+':observed['Na+']+observed['NA'],'Cl-':observed['Cl-']+observed['CL']},
        'initial_concentrations_M':{k:v/(.602214076*np.prod(edges)) for k,v in counts.items()},
        'packing_assumptions':{'protein_partial_specific_volume_ml_g':.73,'initial_solvent_density_g_ml':initial_density},
        'derived_template_charge_rounding':rounding,
        'sha256':{name:file_hash(output/name) for name in ['system.prmtop','system.inpcrd','system_protein.pdb','protonation.json']},
        'request':request,'forcefields':{'protein':'ff19SB','water':'OPC','ligand':'GAFF2'}}
    atomic_json(output/'system_build.json',report)
    return {'output':str(output),**report}
