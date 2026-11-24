"""Restrained assembly-coordinate repair, with no production or equilibrium claim."""
from __future__ import annotations
import json,os,math,shutil,uuid
from pathlib import Path
import numpy as np
from .mechanism import atomic_json
from .datasets import file_hash
from .ligand_parameters import _run
from .structure_review import close_pairs,review_structure


def preparation_pdb_lines(structure,neutral_catalytic_histidines):
    """Rebuild hydrogens after choosing tautomers; never mix old protonation states."""
    lines=[]
    for line in Path(structure).read_text().splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            element=line[76:78].strip() or line[12:16].strip().lstrip('0123456789')[:1]
            if element.upper() in {'H','D'}:continue
            if line[17:20] in {'HIS','HID','HIE','HIP'}:
                variant='HID' if int(line[22:26]) in neutral_catalytic_histidines else 'HIE'
                line=line[:17]+variant+line[20:]
        lines.append(line)
    return lines


def relax_assembly(structure,input,output):
    from openmm import app,unit,CustomExternalForce,LangevinMiddleIntegrator,Platform
    from ..pdb_utils import AA3_TO_AA1
    request=json.loads(Path(input).read_text());sequence=request['sequence'];chains=request['chains']
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    identity={'structure_sha256':file_hash(structure),'request':request}
    if (output/'identity.json').exists() and json.loads((output/'identity.json').read_text())!=identity:raise ValueError('Assembly repair inputs changed')
    if (output/'assembly_relaxation.json').exists():
        report=json.loads((output/'assembly_relaxation.json').read_text())
        if report['output_sha256']!=file_hash(output/'repaired_assembly.pdb'):raise ValueError('Repaired assembly changed')
        return {'output':str(output),'status':report['status'],'production_ready':False}
    if (output/'identity.json').exists():
        history=output/'history'/uuid.uuid4().hex;history.mkdir(parents=True)
        for path in output.iterdir():
            if path.is_file():shutil.copy2(path,history/path.name)
    atomic_json(output/'identity.json',identity)
    amber=Path(os.environ['AQUA_AMBERHOME']);env=os.environ.copy();env['AMBERHOME']=str(amber)
    env['LD_LIBRARY_PATH']=str(amber/'lib');env['PATH']=str(amber/'bin')+os.pathsep+env.get('PATH','')
    lines=preparation_pdb_lines(structure,request.get('neutral_catalytic_histidines',[]))
    (output/'input.pdb').write_text('\n'.join(lines)+'\n')
    leap='source leaprc.protein.ff19SB\nset default PBRadii mbondi2\nprotein = loadpdb input.pdb\ncheck protein\nsaveamberparm protein protein.prmtop protein.inpcrd\nquit\n'
    (output/'repair.leap').write_text(leap);_run([amber/'bin/tleap','-f','repair.leap'],output,env,'repair')
    topology=app.AmberPrmtopFile(str(output/'protein.prmtop'));coords=app.AmberInpcrdFile(str(output/'protein.inpcrd'))
    system=topology.createSystem(nonbondedMethod=app.NoCutoff,implicitSolvent=app.OBC2,constraints=None)
    residues=list(topology.topology.residues());atoms=list(topology.topology.atoms())
    if len(residues)!=len(chains)*len(sequence):raise ValueError('Amber preparation changed residue count')
    for index,residue in enumerate(residues):
        observed=AA3_TO_AA1.get({'HID':'HIS','HIE':'HIS','HIP':'HIS'}.get(residue.name,residue.name))
        if observed!=sequence[index%len(sequence)]:raise ValueError('Amber preparation changed construct identity')
    start=np.array(coords.positions.value_in_unit(unit.nanometer));force=CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    force.addGlobalParameter('k',5000.)
    for name in ['x0','y0','z0']:force.addPerParticleParameter(name)
    core=[]
    for index,atom in enumerate(atoms):
        position=atom.residue.index%len(sequence)+1
        if position>int(request.get('uncovered_n_terminal',0)) and atom.name in {'N','CA','C','O'}:
            force.addParticle(index,start[index]);core.append(index)
    system.addForce(force)
    integrator=LangevinMiddleIntegrator(313.15*unit.kelvin,1/unit.picosecond,.001*unit.picoseconds)
    sim=app.Simulation(topology.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'8'})
    sim.context.setPositions(coords.positions)
    before=sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=15000)
    state=sim.context.getState(getEnergy=True,getPositions=True);xyz=np.array(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    after=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not all(math.isfinite(e) for e in [before,after]) or not np.isfinite(xyz).all():raise ValueError('Assembly minimization is not finite')
    drift=float(np.sqrt(np.mean(np.sum((xyz[core]-start[core])**2,axis=1)))*10)
    output_lines=[];serial=1;last=None;heavy=[];heavy_chains=[];labels=[]
    for index,atom in enumerate(atoms):
        if atom.element.symbol=='H':continue
        chain=chains[atom.residue.index//len(sequence)];position=atom.residue.index%len(sequence)+1
        if last is not None and chain!=last:output_lines.append('TER')
        last=chain;resname={'HID':'HIS','HIE':'HIS','HIP':'HIS'}.get(atom.residue.name,atom.residue.name)
        x,y,z=xyz[index]*10;name=atom.name
        field=' '+name.ljust(3) if len(name)<4 else name
        output_lines.append(f'ATOM  {serial:5d} {field} {resname:3s} {chain}{position:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {atom.element.symbol:>2s}')
        serial+=1;heavy.append([x,y,z]);heavy_chains.append(chain);labels.append(f'{chain}:{position}:{name}')
    pdb=output/'repaired_assembly.pdb';pdb.write_text('\n'.join(output_lines+['TER','END'])+'\n')
    reviews=[review_structure(pdb,sequence,output/('review_'+chain),chain=chain,atom_rules=request.get('atom_rules')) for chain in chains]
    clashes=[{'atoms':[labels[a],labels[b]],'distance_A':float(np.linalg.norm(np.array(heavy[a])-heavy[b]))}
             for a,b in close_pairs(np.array(heavy),1.8) if heavy_chains[a]!=heavy_chains[b]]
    issues=[]
    if drift>1.:issues.append('core_backbone_moved_more_than_1_A')
    if clashes:issues.append('remaining_interchain_clashes')
    if any(not r['structure_review_pass'] for r in reviews):issues.append('chain_geometry_review')
    result={'version':'1.0.0','status':'passed' if not issues else 'blocked','production_ready':False,'blocking_issues':issues,
            'input_sha256':file_hash(structure),'output_sha256':file_hash(pdb),'core_backbone_rms_displacement_A':drift,
            'energy_before_kj_mol':before,'energy_after_kj_mol':after,'interchain_clashes':clashes,
            'method':'ff19SB, OBC2/mbondi2 implicit-solvent coordinate minimization; core backbone restraint 5000 kJ/mol/nm^2',
            'interpretation':'Coordinate repair only; production remains ff19SB/OPC explicit solvent after separate equilibration.',
            'protonation_assumptions':['Catalytic histidines HID; other histidines HIE for geometric repair only.','Production protonation must be recorded independently.'],
            'chain_reviews':reviews,'remaining_checks':['production_protonation','explicit_solution','equilibration']}
    atomic_json(output/'assembly_relaxation.json',result)
    return {'output':str(output),'status':result['status'],'production_ready':False}
