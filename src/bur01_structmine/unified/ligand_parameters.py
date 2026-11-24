"""GAFF2/AM1-BCC preparation with separate numerical and solution-validation states."""
from __future__ import annotations
import json,math,os,re,subprocess,shutil,uuid,zipfile
from pathlib import Path
from .datasets import file_hash
from .mechanism import atomic_json


def acid_fractions(ph,pka):
    """Macroscopic charge-state fractions, not positional microstate populations."""
    if not pka or not all(math.isfinite(float(v)) for v in [ph,*pka]):raise ValueError('Finite pH and pKa values required')
    logs=[0.]
    for value in pka:logs.append(logs[-1]+float(ph)-float(value))
    weights=[10**(x-max(logs)) for x in logs];total=sum(weights)
    return [x/total for x in weights]


def molecule_from_spec(spec):
    from rdkit import Chem
    if isinstance(spec['formal_charge'],bool) or not isinstance(spec['formal_charge'],int):raise ValueError('Formal charge must be an integer')
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,63}',spec['name']):raise ValueError('Invalid molecule identifier')
    mol=Chem.MolFromSmiles(spec['smiles'])
    if mol is None or len(Chem.GetMolFrags(mol))!=1:raise ValueError('Provide one connected molecule')
    if not 1<=mol.GetNumHeavyAtoms()<=100:raise ValueError('Only small molecules up to 100 heavy atoms are supported')
    if any(a.GetSymbol() not in {'H','C','N','O'} for a in mol.GetAtoms()):
        raise ValueError('This adapter covers H/C/N/O molecules; other chemistry needs a reviewed route')
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):raise ValueError('Open-shell species need a separate parameter route')
    if Chem.GetFormalCharge(mol)!=spec['formal_charge']:raise ValueError('SMILES charge and requested charge differ')
    return Chem.AddHs(mol)


def _run(command,cwd,env,label):
    with (cwd/(label+'.log')).open('wb') as log:
        subprocess.run([str(x) for x in command],cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)


def normalize_rounded_charge(path,formal_charge):
    """Bounded uniform correction of small charge-rounding residuals; preserve raw MOL2."""
    path=Path(path);lines=path.read_text().splitlines();indices=[];active=False
    for index,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):active=line=='@<TRIPOS>ATOM';continue
        if active and line.strip():indices.append(index)
    charges=[float(lines[i].split()[8]) for i in indices]
    residual=formal_charge-sum(charges)
    if not indices or abs(residual)>.002 or abs(residual)/len(indices)>1e-4:
        raise ValueError('Charge discrepancy exceeds the permitted rounding correction')
    raw=path.with_name(path.stem+'.raw.mol2');shutil.copy2(path,raw)
    corrected=[round(q+residual/len(charges),8) for q in charges]
    corrected[-1]+=formal_charge-sum(corrected)
    for index,q in zip(indices,corrected):
        fields=lines[index].split();fields[8]=f'{q:.8f}';lines[index]=' '.join(fields)
    path.write_text('\n'.join(lines)+'\n')
    return {'raw_charge_sum':sum(charges),'target_charge':formal_charge,'uniform_correction_per_atom':residual/len(charges),
            'maximum_total_correction_e':.002,'maximum_per_atom_correction_e':1e-4,'raw_mol2_sha256':file_hash(raw)}


def prepare_molecule(spec,output,amber_home):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from openmm import app,unit,Platform,VerletIntegrator
    mol=molecule_from_spec(spec)
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    identity={'version':'1.0.0','molecule':spec,'charge_method':'AM1-BCC','forcefield':'GAFF2'}
    identity_path=output/'identity.json'
    if identity_path.exists() and json.loads(identity_path.read_text())!=identity:
        raise ValueError('Parameter directory belongs to a different molecule or method')
    atomic_json(identity_path,identity)
    report_path=output/'parameter_report.json'
    if report_path.exists():
        old=json.loads(report_path.read_text())
        if old.get('status')=='numerical_checks_passed' and all(file_hash(output/k)==v for k,v in old['sha256'].items()):return old
        history=output/'history'/uuid.uuid4().hex;history.mkdir(parents=True)
        for p in output.iterdir():
            if p.is_file():shutil.copy2(p,history/p.name)
    amber=Path(amber_home).resolve()
    env=os.environ.copy();env['AMBERHOME']=str(amber);env['PATH']=str(amber/'bin')+os.pathsep+env.get('PATH','')
    env['LD_LIBRARY_PATH']=str(amber/'lib')
    env['OMP_NUM_THREADS']='2'
    for executable in ['antechamber','parmchk2','tleap']:
        if not (amber/'bin'/executable).is_file():raise ValueError('Missing AmberTools executable: '+executable)
    params=AllChem.ETKDGv3();params.randomSeed=20261124
    conformers=list(AllChem.EmbedMultipleConfs(mol,numConfs=12,params=params))
    if not conformers:raise ValueError('Molecule embedding failed')
    energies=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000)
    valid=[(e,i) for i,(status,e) in enumerate(energies) if status==0 and math.isfinite(e)]
    if not valid:raise ValueError('No converged starting conformer')
    _,best=min(valid)
    with Chem.SDWriter(str(output/'input.sdf')) as writer:writer.write(mol,confId=conformers[best])
    _run([amber/'bin/antechamber','-i','input.sdf','-fi','sdf','-o','ligand.mol2','-fo','mol2',
          '-c','bcc','-at','gaff2','-nc',str(spec['formal_charge']),'-rn','LIG','-pf','y','-dr','yes'],output,env,'antechamber')
    rounding=normalize_rounded_charge(output/'ligand.mol2',spec['formal_charge'])
    _run([amber/'bin/parmchk2','-i','ligand.mol2','-f','mol2','-o','ligand.frcmod','-s','gaff2'],output,env,'parmchk2')
    frcmod=(output/'ligand.frcmod').read_text()
    unresolved=[line for line in frcmod.splitlines() if any(s in line.upper() for s in ['ATTN','NEED REVISION'])]
    leap='source leaprc.gaff2\nsource leaprc.water.opc\nloadamberparams ligand.frcmod\nLIG = loadmol2 ligand.mol2\ncheck LIG\nsaveamberparm LIG ligand.prmtop ligand.inpcrd\nsavepdb LIG ligand.pdb\nquit\n'
    (output/'leap.in').write_text(leap)
    _run([amber/'bin/tleap','-f','leap.in'],output,env,'tleap')
    text=(output/'tleap.log').read_text(errors='replace')
    if re.search(r'(?i)(FATAL|could not find|does not have a type|not saved)',text):unresolved.append('tleap reported unresolved parameters')
    prmtop=app.AmberPrmtopFile(str(output/'ligand.prmtop'));coords=app.AmberInpcrdFile(str(output/'ligand.inpcrd'))
    system=prmtop.createSystem(nonbondedMethod=app.NoCutoff,constraints=None)
    atoms=list(prmtop.topology.atoms())
    if len(atoms)!=mol.GetNumAtoms():unresolved.append('Atom count differs from input molecule')
    if sorted(a.element.symbol for a in atoms)!=sorted(a.GetSymbol() for a in mol.GetAtoms()):unresolved.append('Element counts differ from input')
    if len(list(prmtop.topology.bonds()))!=mol.GetNumBonds():unresolved.append('Bond count differs from input')
    nonbonded=next(f for f in system.getForces() if f.__class__.__name__=='NonbondedForce')
    charge=sum(nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
    if abs(charge-spec['formal_charge'])>1e-4:unresolved.append('Partial charges do not sum to formal charge')
    integrator=VerletIntegrator(.001*unit.picoseconds)
    simulation=app.Simulation(prmtop.topology,system,integrator,Platform.getPlatformByName('Reference'))
    simulation.context.setPositions(coords.positions)
    initial=simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    simulation.minimizeEnergy(maxIterations=2000)
    final=simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not all(math.isfinite(e) for e in [initial,final]):unresolved.append('Non-finite potential energy')
    names=['input.sdf','ligand.mol2','ligand.raw.mol2','ligand.frcmod','ligand.prmtop','ligand.inpcrd','ligand.pdb']
    ff=amber/'dat/leap/parm/gaff2.dat'
    metadata=next((amber/'conda-meta').glob('ambertools-*.json'),None)
    result={'version':'1.0.0','name':spec['name'],'status':'blocked' if unresolved else 'numerical_checks_passed',
            'release_status':'solution_validation_required','unresolved':unresolved,'formal_charge':spec['formal_charge'],
            'partial_charge_sum':charge,'atoms':len(atoms),'bonds':mol.GetNumBonds(),
            'initial_energy_kj_mol':initial,'minimized_energy_kj_mol':final,
            'sha256':{name:file_hash(output/name) for name in names},'gaff2_sha256':file_hash(ff),
            'ambertools':json.loads(metadata.read_text()).get('version') if metadata else 'unrecorded',
            'charge_rounding':rounding,
            'required_solution_checks':['aqueous_density','aggregation','realized_concentration','protonation_population_assumptions']}
    atomic_json(report_path,result)
    if unresolved:raise ValueError('Parameter checks failed: '+', '.join(unresolved))
    return result


def package_parameters(output):
    """Export the checked parameter files and their original numerical reports."""
    output=Path(output)
    report_path=output/'parameterization_report.json'
    report=json.loads(report_path.read_text())
    if report['status']!='numerical_checks_passed':raise ValueError('Parameter set is incomplete')
    paths=[report_path]
    for molecule in report['molecules']:
        folder=output/molecule['name'];paths.append(folder/'parameter_report.json')
        if folder.parent!=output or molecule['status']!='numerical_checks_passed':raise ValueError('Invalid parameter set')
        for name,digest in molecule['sha256'].items():
            path=folder/name
            if path.parent!=folder or file_hash(path)!=digest:raise ValueError('Parameter checksum mismatch')
            paths.append(path)
    temporary=output/'parameters.zip.tmp'
    with zipfile.ZipFile(temporary,'w',zipfile.ZIP_DEFLATED) as archive:
        for path in paths:archive.write(path,path.relative_to(output).as_posix())
    temporary.replace(output/'parameters.zip')
    return output/'parameters.zip'


def parameterize(input,output):
    request=json.loads(Path(input).read_text(encoding='utf-8'))
    if request.get('forcefield')!='GAFF2' or request.get('charge_method')!='AM1-BCC':raise ValueError('Unsupported parameter route')
    molecules=request.get('molecules',[])
    if not 1<=len(molecules)<=50 or len({s['name'] for s in molecules})!=len(molecules):raise ValueError('Provide 1–50 molecules with distinct names')
    home=os.environ.get('AQUA_AMBERHOME')
    if not home:raise ValueError('Configure AQUA_AMBERHOME on this compute server')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    results=[]
    for spec in request['molecules']:
        results.append(prepare_molecule(spec,output/spec['name'],home))
        atomic_json(output/'parameterization_report.json',{'status':'numerical_checks_passed' if len(results)==len(molecules) else 'incomplete',
                    'production_ready':False,'requested_molecules':len(molecules),'molecules':results})
    package_parameters(output)
    return {'output':str(output),'molecules':len(results),'status':'numerical_checks_passed','production_ready':False}
