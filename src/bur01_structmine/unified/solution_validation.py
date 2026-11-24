"""Bounded bulk-solution checks; measured density agreement is reported separately."""
from __future__ import annotations
import json,math,os,subprocess,shutil,uuid,re
from pathlib import Path
import numpy as np
from .datasets import file_hash
from .mechanism import atomic_json,OPENMM_VERSION
from .ligand_parameters import _run


def validate_solution(input,parameter_dir,output):
    from importlib.metadata import version
    import openmm
    from openmm import app,unit
    if version('OpenMM')!=OPENMM_VERSION:raise ValueError('Use the locked production OpenMM environment')
    request=json.loads(Path(input).read_text());output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    parameter_dir=Path(parameter_dir).resolve();amber=Path(os.environ['AQUA_AMBERHOME'])
    duration=float(request.get('npt_ns',1.));length=float(request.get('box_nm',7.))
    if not .02<=duration<=2 or not 3<=length<=10:raise ValueError('Bulk validation requires 0.02–2 ns and a 3–10 nm box')
    if abs(duration/.002-round(duration/.002))>1e-8:raise ValueError('Use a multiple of the 2 ps reporting interval')
    spec={'request':request,'parameter_files':{}}
    env=os.environ.copy();env['AMBERHOME']=str(amber);env['PATH']=str(amber/'bin')+os.pathsep+env.get('PATH','')
    env['LD_LIBRARY_PATH']=str(amber/'lib');env['OMP_NUM_THREADS']='2'
    commands=['source leaprc.gaff2','source leaprc.water.opc','savepdb OPCBOX opc_box.pdb']
    counts={};mass_solutes=0.;net_charge=0
    for index,item in enumerate(request['components']):
        name=item['name'];count=int(item['count']);folder=parameter_dir/name
        if not 0<count<=1000 or folder.parent!=parameter_dir:raise ValueError('Invalid component')
        report=json.loads((folder/'parameter_report.json').read_text())
        if report['status']!='numerical_checks_passed':raise ValueError('Component parameters have not passed numerical checks')
        for filename,digest in report['sha256'].items():
            if file_hash(folder/filename)!=digest:raise ValueError('Component parameters changed')
        spec['parameter_files'][name]=file_hash(folder/'parameter_report.json')
        if index>=26:raise ValueError('At most 26 solute species are supported')
        residue='QA'+chr(65+index);counts[residue]={'name':name,'count':count}
        prmtop=app.AmberPrmtopFile(str(folder/'ligand.prmtop'));system=prmtop.createSystem(nonbondedMethod=app.NoCutoff)
        mass_solutes+=count*sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
        net_charge+=count*report['formal_charge']
        # Rename the molecule and residue consistently in the intermediate template.
        mol2=(folder/'ligand.mol2').read_text()
        (output/(residue+'.mol2')).write_text(re.sub(r'\bLIG\b',residue,mol2))
        commands += [f'loadamberparams "{folder / "ligand.frcmod"}"',f'{residue} = loadmol2 {residue}.mol2',
                     f'savepdb {residue} {residue}.pdb']
    identity=output/'identity.json'
    if identity.exists() and json.loads(identity.read_text())!=spec:raise ValueError('Bulk-validation output belongs to another specification')
    if (output/'solution_validation.json').is_file():
        result=json.loads((output/'solution_validation.json').read_text())
        return {'output':str(output),**result}
    if identity.exists():
        history=output/'history'/uuid.uuid4().hex;history.mkdir(parents=True)
        for p in output.iterdir():
            if p.is_file():shutil.copy2(p,history/p.name)
    atomic_json(identity,spec)
    (output/'units.leap').write_text('\n'.join(commands+['quit'])+'\n');_run([amber/'bin/tleap','-f','units.leap'],output,env,'units')
    for residue in counts:
        names={line[17:20].strip() for line in (output/(residue+'.pdb')).read_text().splitlines() if line.startswith(('ATOM  ','HETATM'))}
        if names!={residue}:raise ValueError('Prepared residue name differs from its Amber template')
    water=[];first=None
    for line in (output/'opc_box.pdb').read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            key=line[21:27]
            if first is None:first=key
            if key!=first:break
            water.append(line)
    if len(water)!=4:raise ValueError('Expected a four-site OPC water template')
    (output/'water.pdb').write_text('\n'.join(water+['TER','END'])+'\n')
    target_density=float(request.get('initial_density_g_ml',.984))
    water_count=int(round((target_density*602.214076*length**3-mass_solutes-abs(net_charge)*22.99)/18.01528))
    if water_count<100:raise ValueError('Insufficient water for this mixture')
    pack=['tolerance 2.1','filetype pdb','output packed.pdb','seed '+str(request.get('seed',20261124)),'nloop 1000']
    edge=length*10
    for residue,item in counts.items():pack += [f'structure {residue}.pdb',f'number {item["count"]}',f'inside box 1.5 1.5 1.5 {edge-1.5} {edge-1.5} {edge-1.5}','end structure']
    pack += ['structure water.pdb',f'number {water_count}',f'inside box 1.1 1.1 1.1 {edge-1.1} {edge-1.1} {edge-1.1}','end structure']
    (output/'packmol.inp').write_text('\n'.join(pack)+'\n')
    with (output/'packmol.inp').open('rb') as source,(output/'packmol.log').open('wb') as log:
        subprocess.run([str(amber/'bin/packmol')],stdin=source,stdout=log,stderr=subprocess.STDOUT,cwd=output,env=env,check=True,timeout=900)
    # Every packed molecule is independent; no inferred inter-residue bonds.
    packed=[];last=None
    for line in (output/'packed.pdb').read_text().splitlines():
        if not line.startswith(('ATOM  ','HETATM')):continue
        key=line[17:27]
        if last is not None and last!=key:packed.append('TER')
        packed.append(line);last=key
    (output/'packed_terminated.pdb').write_text('\n'.join(packed+['TER','END'])+'\n')
    # AmberTools 26 savepdb mutates the loaded ligand template in this workflow.
    # Freshly reload templates and do not save individual units before loadpdb.
    build=[c for c in commands if not c.startswith('savepdb')]+["solution = loadpdb packed_terminated.pdb",f'set solution box {{ {edge} {edge} {edge} }}']
    if net_charge:build += ['addionsrand solution '+('Na+' if net_charge<0 else 'Cl-')+' 0']
    build += ['check solution','saveamberparm solution solution.prmtop solution.inpcrd','savepdb solution solution.pdb','quit']
    (output/'solution.leap').write_text('\n'.join(build)+'\n');_run([amber/'bin/tleap','-f','solution.leap'],output,env,'solution')
    topology=app.AmberPrmtopFile(str(output/'solution.prmtop'));coords=app.AmberInpcrdFile(str(output/'solution.inpcrd'))
    system=topology.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1*unit.nanometer,constraints=app.HBonds,rigidWater=True)
    temperature=float(request['temperature_K']);seed=int(request.get('seed',20261124))
    barostat=openmm.MonteCarloBarostat(1*unit.bar,temperature*unit.kelvin,25);barostat.setRandomNumberSeed(seed+1)
    system.addForce(barostat);barostat.setFrequency(0)
    integrator=openmm.LangevinMiddleIntegrator(temperature*unit.kelvin,1/unit.picosecond,.002*unit.picoseconds);integrator.setRandomNumberSeed(seed)
    sim=app.Simulation(topology.topology,system,integrator,openmm.Platform.getPlatformByName('CUDA'),{'Precision':'mixed','DeviceIndex':'0'})
    sim.context.setPositions(coords.positions);sim.context.setPeriodicBoxVectors(*coords.boxVectors)
    sim.minimizeEnergy(maxIterations=5000);sim.context.setVelocitiesToTemperature(temperature*unit.kelvin,seed)
    sim.step(20000)  # 40 ps NVT, excluded from NPT sampling.
    barostat.setFrequency(25);sim.context.reinitialize(preserveState=True)
    mass=sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
    count_steps=int(round(duration*500000));stride=1000;rows=[]
    sim.reporters=[app.DCDReporter(str(output/'bulk.dcd'),5000,enforcePeriodicBox=True)]
    for end in range(stride,count_steps+1,stride):
        sim.step(stride);state=sim.context.getState(getEnergy=True)
        volume=state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
        energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not all(math.isfinite(x) for x in [volume,energy]) or volume<=0:raise ValueError('Unstable bulk simulation')
        rows.append({'npt_ps':end*.002,'volume_nm3':volume,'density_g_ml':mass/(602.214076*volume),'energy_kj_mol':energy})
        if end%50000==0:
            atomic_json(output/'progress.json',{'npt_ns_completed':end*2e-6,'approved_npt_ns':duration,'atoms':system.getNumParticles()})
            print(f'Bulk NPT {end*2e-6:.1f}/{duration:g} ns; density {rows[-1]["density_g_ml"]:.5f} g/mL',flush=True)
    sim.saveState(str(output/'final_state.xml'));sim.saveCheckpoint(str(output/'final.chk'))
    density=np.array([r['density_g_ml'] for r in rows[len(rows)//2:]])
    blocks=[float(x.mean()) for x in np.array_split(density,4)]
    drift=abs(blocks[-1]-blocks[0])/float(density.mean());volume=float(np.mean([r['volume_nm3'] for r in rows[len(rows)//2:]]))
    concentrations={item['name']:item['count']/(.602214076*volume) for item in counts.values()}
    reference=request.get('density_reference');comparison=None
    if reference:
        if abs(reference['temperature_K']-temperature)>1e-8 or not reference.get('source_url'):raise ValueError('Density reference must match temperature and include its source')
        comparison={'reference':reference,'relative_error':abs(float(density.mean())/reference['density_g_ml']-1),
                    'within_tolerance':abs(float(density.mean())/reference['density_g_ml']-1)<=reference['relative_tolerance']}
    report={'version':'1.0.0','status':'numerical_checks_passed' if drift<.02 else 'density_drift_requires_longer_sampling',
            'production_ready':False,'npt_ns':duration,'nvt_ps':40,'temperature_K':temperature,'atoms':system.getNumParticles(),
            'component_counts':counts,'water_molecules':water_count,'counterion_charge':-net_charge,'mean_volume_nm3':volume,
            'density_mean_g_ml':float(density.mean()),'density_block_means':blocks,'relative_density_drift':drift,
            'realized_concentrations_M':concentrations,'density_comparison':comparison,
            'protonation_assumptions':request.get('protonation_assumptions',[]),
            'remaining_checks':['aggregation_and_bulk_distribution','microstate_population_sensitivity','protein_system_equilibration']}
    atomic_json(output/'timeseries.json',rows);atomic_json(output/'solution_validation.json',report)
    return {'output':str(output),**report}
