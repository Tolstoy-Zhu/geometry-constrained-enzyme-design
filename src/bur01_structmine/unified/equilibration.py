"""Bounded, restartable explicit-solvent equilibration. Failed checks never release MD."""
from __future__ import annotations
import json,math,time
from pathlib import Path
from importlib.metadata import version
import numpy as np
from .datasets import file_hash
from .mechanism import atomic_json,OPENMM_VERSION,runtime_fingerprint,checked_file
from .protein_system import protein_lines,verify_construct


def equilibration_schedule(protocol):
    dt=float(protocol['timestep_fs']);chunk=float(protocol['checkpoint_ps']);temperature=float(protocol['temperature_K'])
    if dt!=2. or chunk!=50.:raise ValueError('This protocol uses 2 fs steps and 50 ps checkpoints')
    if not math.isfinite(temperature) or not 250<=temperature<=400:raise ValueError('Unsupported equilibration temperature')
    stages=[('nvt',float(protocol['nvt_ps']),None)]
    ks=protocol['restraint_kj_mol_nm2']
    if not ks or any(not math.isfinite(k) or k<0 for k in ks):raise ValueError('Invalid restraint schedule')
    stages += [('restrained_npt',float(protocol['restrained_npt_ps'])/len(ks),float(k)) for k in ks]
    stages += [('unrestrained_npt',float(protocol['unrestrained_npt_ps']),0.)]
    if sum(ps for _,ps,_ in stages)>3000:raise ValueError('Equilibration exceeds the 3 ns per-system preparation limit')
    if float(protocol['unrestrained_npt_ps'])<500:raise ValueError('At least 500 ps unrestrained preparation is required')
    result=[];start=0
    for kind,ps,k in stages:
        if not math.isfinite(ps) or ps<=0 or abs(ps/chunk-round(ps/chunk))>1e-8:raise ValueError('Stages must be positive multiples of the checkpoint interval')
        for index in range(round(ps/chunk)):
            target=100+(temperature-100)*(index+1)/(ps/chunk) if kind=='nvt' else temperature
            steps=round(chunk*1000/dt)
            result.append({'stage':kind,'start_step':start,'end_step':start+steps,
                           'temperature_K':target,'restraint_k':float(ks[0]) if k is None else k})
            start+=steps
    return result


def aligned_rmsd(first,second):
    a=np.asarray(first)-np.mean(first,axis=0);b=np.asarray(second)-np.mean(second,axis=0)
    u,_,v=np.linalg.svd(a.T@b);correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@v)
    return float(np.sqrt(np.mean(np.sum((a@u@correction@v-b)**2,axis=1))))


def write_amber_restart(path,positions_A,box_A):
    xyz=np.asarray(positions_A);box=np.asarray(box_A)
    if xyz.ndim!=2 or xyz.shape[1]!=3 or not np.isfinite(xyz).all():raise ValueError('Invalid restart coordinates')
    lengths=np.linalg.norm(box,axis=1)
    if not np.allclose(box,np.diag(lengths),atol=1e-7):raise ValueError('This writer supports orthorhombic boxes')
    lines=['AquaEsterNet equilibrated coordinates',f'{len(xyz):6d}']
    values=xyz.reshape(-1)
    for start in range(0,len(values),6):lines.append(''.join(f'{x:12.7f}' for x in values[start:start+6]))
    lines.append(''.join(f'{x:12.7f}' for x in [*lengths,90.,90.,90.]))
    Path(path).write_text('\n'.join(lines)+'\n',encoding='ascii')


def assess_equilibration(rows,protocol):
    tail=[r for r in rows if r['stage']=='unrestrained_npt']
    if not tail:raise ValueError('No unrestrained sampling')
    tail=tail[len(tail)//2:];qc=protocol['qc'];issues=[]
    for row in tail:
        values=[row['temperature_K'],row['volume_nm3'],row['core_rmsd_A'],*row['geometry_A'].values()]
        if not all(math.isfinite(v) for v in values) or row['volume_nm3']<=0:
            return {'status':'blocked','blocking_issues':['nonfinite_or_invalid_tail_observables']}
    temperature=float(np.mean([r['temperature_K'] for r in tail]));vol=np.array([r['volume_nm3'] for r in tail])
    cv=float(vol.std()/vol.mean());rmsd=max(r['core_rmsd_A'] for r in tail)
    if abs(temperature-protocol['temperature_K'])>qc['temperature_mean_tolerance_K']:issues.append('tail_temperature')
    if cv>qc['volume_cv_max']:issues.append('tail_volume_variation')
    if rmsd>qc['core_rmsd_max_A']:issues.append('core_rmsd')
    geometry={}
    for name in tail[0]['geometry_A']:
        values=np.array([r['geometry_A'][name] for r in tail])
        valid=(values>=qc['geometry_min_A'])&(values<=qc['geometry_max_A'])
        geometry[name]={'final_A':float(values[-1]),'mean_A':float(values.mean()),'in_range_fraction':float(valid.mean())}
        if not valid[-1]:issues.append('final_geometry:'+name)
    return {'status':'passed' if not issues else 'blocked','blocking_issues':issues,
            'tail_temperature_K':temperature,'tail_volume_cv':cv,'tail_max_core_rmsd_A':rmsd,'geometry':geometry,
            'thresholds':qc,'tail_samples':len(tail),'tail_mean_volume_nm3':float(vol.mean())}


def equilibrate_system(system_dir,output,prerequisite_dir=None):
    import openmm
    from openmm import app,unit
    if version('OpenMM')!=OPENMM_VERSION:raise ValueError('Use the locked OpenMM environment')
    prerequisite_sha256=None
    if prerequisite_dir:
        prerequisite=Path(prerequisite_dir)/'equilibration_report.json'
        evidence=json.loads(prerequisite.read_text())
        if evidence.get('status')!='passed' or evidence.get('blocking_issues'):
            raise ValueError('Parent equilibration has not passed; candidate equilibration is held')
        prerequisite_sha256=file_hash(prerequisite)
    base=Path(system_dir).resolve();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    report=json.loads((base/'system_build.json').read_text())
    for name,digest in report['sha256'].items():
        if file_hash(base/name)!=digest:raise ValueError('Prepared system changed')
    request=report['request'];protocol=request['protocol'];schedule=equilibration_schedule(protocol)
    sequences=request['construct_sequences'];runtime=runtime_fingerprint('CUDA')
    identity={'system_build_sha256':file_hash(base/'system_build.json'),'protocol':protocol,'runtime':runtime,
              'seed':request['seed'],'equilibration_version':'1.0.0','prerequisite_sha256':prerequisite_sha256}
    progress_file=output/'progress.json';progress={'identity':identity,'completed_chunks':0,'chunks':[]}
    if progress_file.exists():
        progress=json.loads(progress_file.read_text())
        if progress['identity']!=identity:raise ValueError('Equilibration environment or inputs changed; binary checkpoint is not portable')
    topology=app.AmberPrmtopFile(str(base/'system.prmtop'));coords=app.AmberInpcrdFile(str(base/'system.inpcrd'))
    verify_construct(topology.topology,sequences);initial=np.array(coords.positions.value_in_unit(unit.nanometer))
    system=topology.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1*unit.nanometer,constraints=app.HBonds,rigidWater=True)
    barostat=openmm.MonteCarloBarostat(protocol['pressure_bar']*unit.bar,protocol['temperature_K']*unit.kelvin,25)
    barostat.setRandomNumberSeed(request['seed']+1);system.addForce(barostat)
    next_chunk=schedule[min(progress['completed_chunks'],len(schedule)-1)]
    barostat.setFrequency(0 if next_chunk['stage']=='nvt' else 25)
    restraint=openmm.CustomExternalForce('0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2');restraint.addGlobalParameter('k',1000.)
    for name in ['x0','y0','z0']:restraint.addPerParticleParameter(name)
    atoms=list(topology.topology.atoms());protein_count=report['protein_residues'];residues=list(topology.topology.residues())
    core=[];lookup={};offset=0
    for chain,sequence in sequences.items():
        for pos in range(1,len(sequence)+1):
            for atom in residues[offset+pos-1].atoms():lookup[(chain,pos,atom.name)]=atom.index
        offset+=len(sequence)
    for atom in atoms:
        if atom.residue.index>=protein_count:continue
        if atom.element and atom.element.symbol!='H':restraint.addParticle(atom.index,initial[atom.index])
    for (chain,pos,name),index in lookup.items():
        if pos>=protocol.get('core_first_position',1) and name=='CA':core.append(index)
    system.addForce(restraint)
    integrator=openmm.LangevinMiddleIntegrator(next_chunk['temperature_K']*unit.kelvin,1/unit.picosecond,.002*unit.picosecond)
    integrator.setRandomNumberSeed(request['seed'])
    simulation=app.Simulation(topology.topology,system,integrator,openmm.Platform.getPlatformByName('CUDA'),{'Precision':'mixed','DeviceIndex':'0'})
    simulation.context.setPositions(coords.positions);simulation.context.setPeriodicBoxVectors(*coords.boxVectors)
    if progress['completed_chunks']:
        simulation.loadCheckpoint(str(checked_file(progress['checkpoint'],output)))
    else:
        print('CUDA minimization started',flush=True)
        simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=10000)
        simulation.context.setVelocitiesToTemperature(100*unit.kelvin,request['seed'])
        print('CUDA minimization completed; starting bounded equilibration',flush=True)
    mass=sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
    dof=3*sum(system.getParticleMass(i)>0*unit.dalton for i in range(system.getNumParticles()))-system.getNumConstraints()-3
    report_steps=round(protocol['report_ps']*500)
    if report_steps<=0 or round(protocol['checkpoint_ps']*500)%report_steps:raise ValueError('Reporting interval must divide checkpoint interval')
    for index in range(progress['completed_chunks'],len(schedule)):
        chunk=schedule[index];name=f'chunk_{index:04d}';directory=output/name
        if directory.exists():directory=output/(name+'_'+__import__('uuid').uuid4().hex[:8])
        directory.mkdir();started=time.monotonic();frequency=0 if chunk['stage']=='nvt' else 25
        if frequency!=barostat.getFrequency():
            barostat.setFrequency(frequency);simulation.context.reinitialize(preserveState=True)
        simulation.context.setParameter('k',chunk['restraint_k']);integrator.setTemperature(chunk['temperature_K']*unit.kelvin)
        simulation.reporters=[app.DCDReporter(str(directory/'trajectory.dcd'),report_steps,enforcePeriodicBox=True)]
        rows=[]
        for end in range(chunk['start_step']+report_steps,chunk['end_step']+1,report_steps):
            simulation.step(report_steps)
            state=simulation.context.getState(getPositions=True,getEnergy=True,enforcePeriodicBox=False)
            xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
            box=np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer));lengths=np.diag(box)
            volume=state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
            energy=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            kinetic=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
            if not np.isfinite(xyz).all() or not all(math.isfinite(v) for v in [volume,energy,kinetic]) or volume<=0:
                raise ValueError('Nonfinite equilibration state')
            geometry={}
            for chain in sequences:
                for pair in protocol['atom_pairs']:
                    first=[lookup[(chain,pair['positions'][0],n)] for n in pair['first_atoms']]
                    second=[lookup[(chain,pair['positions'][1],n)] for n in pair['second_atoms']]
                    delta=xyz[first][:,None]-xyz[second][None];delta-=lengths*np.round(delta/lengths)
                    geometry[chain+':'+':'.join(map(str,pair['positions']))]=float(np.linalg.norm(delta,axis=-1).min()*10)
            rows.append({'step':end,'ps':end*.002,'stage':chunk['stage'],'energy_kj_mol':energy,
                'temperature_K':2*kinetic/(dof*.00831446261815324),'volume_nm3':volume,'density_g_ml':mass/(602.214076*volume),
                'core_rmsd_A':aligned_rmsd(xyz[core],initial[core])*10,'geometry_A':geometry})
        simulation.reporters=[];atomic_json(directory/'samples.json',rows)
        checkpoint=directory/'state.chk';simulation.saveCheckpoint(str(checkpoint))
        progress['chunks'].append({'directory':directory.name,'stage':chunk['stage'],'last_sample':rows[-1],
            'samples_sha256':file_hash(directory/'samples.json'),'wall_seconds':time.monotonic()-started})
        progress.update(completed_chunks=index+1,completed_ps=chunk['end_step']*.002,
                        checkpoint={'path':str(checkpoint.relative_to(output)),'sha256':file_hash(checkpoint)})
        atomic_json(progress_file,progress)
        print(f"{chunk['stage']} {progress['completed_ps']:g}/{schedule[-1]['end_step']*.002:g} ps; T={rows[-1]['temperature_K']:.2f} K; core RMSD={rows[-1]['core_rmsd_A']:.3f} A",flush=True)
    all_rows=[]
    for chunk in progress['chunks']:
        path=output/chunk['directory']/'samples.json'
        if file_hash(path)!=chunk['samples_sha256']:raise ValueError('Equilibration samples changed')
        all_rows+=json.loads(path.read_text())
    qc=assess_equilibration(all_rows,protocol)
    state=simulation.context.getState(getPositions=True,getEnergy=True,enforcePeriodicBox=False)
    xyz=np.asarray(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom));box=np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom))
    write_amber_restart(output/'equilibrated.inpcrd',xyz,box)
    (output/'equilibrated_protein.pdb').write_text('\n'.join(protein_lines(topology.topology,xyz,sequences))+'\n')
    simulation.saveState(str(output/'equilibrated_state.xml'))
    from .coordinates import validate_amber_coordinates
    coordinate_check=validate_amber_coordinates(base/'system.prmtop',output/'equilibrated.inpcrd',output/'equilibrated_protein.pdb',sequences)
    result={'schema_version':'1.0.0',**qc,'coordinate_validation':coordinate_check,'production_ready':False,
        'equilibration_ps':schedule[-1]['end_step']*.002,'production_ns':0.,'system_build_sha256':identity['system_build_sha256'],
        'protonation_mode':request['protonation_mode'],'runtime':runtime,
        'realized_concentrations_M':{r['name']:r['count']/(.602214076*qc['tail_mean_volume_nm3']) for r in report['component_counts'].values()},
        'remaining_checks':['independent_equilibration_review','production_batch_release'],
        'interpretation':'Finite fixed-microstate preparation check, not a microstate population or converged mechanistic result.'}
    atomic_json(output/'equilibration_report.json',result)
    return {'output':str(output),**result}
