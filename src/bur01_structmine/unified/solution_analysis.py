"""Periodic bulk-solution descriptors, kept separate from protein-mechanism conclusions."""
from __future__ import annotations
import csv,json
from pathlib import Path
import numpy as np
from .datasets import file_hash
from .mechanism import atomic_json


def periodic_distances(first,second,lengths):
    lengths=np.asarray(lengths,float)
    if lengths.shape!=(3,) or not np.isfinite(lengths).all() or np.any(lengths<=0):
        raise ValueError('Three positive periodic box lengths are required')
    delta=np.asarray(first,float)[:,None,:]-np.asarray(second,float)[None,:,:]
    delta-=lengths*np.rint(delta/lengths)
    return np.linalg.norm(delta,axis=-1)


def largest_cluster(size,edges):
    parents=list(range(size))
    def root(index):
        while parents[index]!=index:
            parents[index]=parents[parents[index]];index=parents[index]
        return index
    for a,b in edges:
        parents[root(a)]=root(b)
    counts={}
    for index in range(size):
        key=root(index);counts[key]=counts.get(key,0)+1
    return max(counts.values(),default=0)


def solution_status(report):
    failures=[]
    if report['relative_density_drift']>=.02:failures.append('density_drift')
    comparison=report.get('density_comparison')
    if comparison and not comparison['within_tolerance']:failures.append('density_reference_mismatch')
    return {'status':'passed' if not failures else 'blocked','blocking_issues':failures}


def analyze_solution(run,output):
    import MDAnalysis as mda
    source=Path(run).resolve();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    report=json.loads((source/'solution_validation.json').read_text())
    universe=mda.Universe(str(source/'solution.prmtop'),str(source/'bulk.dcd'))
    components=report['component_counts']
    solutes={name:[] for name in (c['name'] for c in components.values())}
    for residue in universe.residues:
        if residue.resname in components:
            name=components[residue.resname]['name']
            solutes[name].append(residue.atoms.indices[residue.atoms.masses>1.5])
    for residue,component in components.items():
        if len(solutes[component['name']])!=component['count']:raise ValueError('Topology differs from recorded molecule counts')
    ethanol=[]
    for indices in solutes.get('ethanol',[]):
        oxygen=[i for i in indices if 15.5<universe.atoms[i].mass<16.5]
        if len(oxygen)!=1:raise ValueError('Expected one oxygen per ethanol molecule')
        ethanol.extend(oxygen)
    citrate=[atoms for name,groups in solutes.items() if name.startswith('citrate_') for atoms in groups]
    octanoic=[atoms for name,groups in solutes.items() if name in {'octanoic_acid','octanoate'} for atoms in groups]
    sodium=np.flatnonzero((universe.atoms.masses>22.5)&(universe.atoms.masses<23.5))
    citrate_oxygens=[i for group in citrate for i in group if 15.5<universe.atoms[i].mass<16.5]
    frames=[];rdf_counts=np.zeros(50);rdf_expected=np.zeros(50);bins=np.linspace(0,25,51)
    pair_contacts={};start=len(universe.trajectory)//2
    if len(universe.trajectory)-start<10:raise ValueError('At least ten saved analysis frames are required')
    for ts in universe.trajectory[start:]:
        if not np.isfinite(ts.positions).all():raise ValueError('Non-finite bulk trajectory coordinates')
        box=ts.dimensions[:3].astype(float)
        if not np.allclose(ts.dimensions[3:],90,atol=.01):raise ValueError('This bulk analysis requires an orthorhombic box')
        if min(box)<=2*bins[-1]:raise ValueError('RDF range exceeds half the shortest box length')
        row={'time_ps':float(ts.time),'citrate_largest_cluster':0,'octanoic_acid_contact':False}
        if len(ethanol)>1:
            xyz=ts.positions[ethanol]
            distances=periodic_distances(xyz,xyz,box)
            unique=distances[np.triu_indices(len(xyz),1)]
            rdf_counts+=np.histogram(unique,bins)[0]
            shell=4*np.pi/3*(bins[1:]**3-bins[:-1]**3)
            rdf_expected+=len(xyz)*(len(xyz)-1)/2*shell/np.prod(box)
            # Oxygen-based eight-cell occupancy is a finite-box descriptor.
            cells=np.floor((xyz%box)/(box/2)).astype(int)
            counts=np.bincount(cells[:,0]*4+cells[:,1]*2+cells[:,2],minlength=8)
            row['ethanol_octant_counts']=counts.tolist()
            row['ethanol_octant_cv']=float(counts.std()/counts.mean())
        edges=[]
        for i in range(len(citrate)):
            for j in range(i+1,len(citrate)):
                contact=bool(periodic_distances(ts.positions[citrate[i]],ts.positions[citrate[j]],box).min()<4.)
                pair_contacts.setdefault(f'{i}-{j}',[]).append(contact)
                if contact:edges.append((i,j))
        row['citrate_largest_cluster']=largest_cluster(len(citrate),edges)
        row['citrate_contact_pairs']=len(edges)
        if len(octanoic)>=2:
            row['octanoic_acid_contact']=any(periodic_distances(ts.positions[octanoic[i]],ts.positions[octanoic[j]],box).min()<4.
                                              for i in range(len(octanoic)) for j in range(i+1,len(octanoic)))
        if len(sodium) and citrate_oxygens:
            near=periodic_distances(ts.positions[sodium],ts.positions[citrate_oxygens],box).min(axis=1)<3.2
            row['sodium_near_citrate_fraction']=float(near.mean())
        frames.append(row)
    if len(ethanol)>1:
        with (output/'ethanol_oxygen_rdf.csv').open('w',newline='') as handle:
            writer=csv.writer(handle);writer.writerow(['r_A','g_r','pair_observations'])
            writer.writerows(zip((bins[1:]+bins[:-1])/2,rdf_counts/rdf_expected,rdf_counts.astype(int)))
    result={'schema_version':'1.0.0',**solution_status(report),'production_ready':False,
            'scope':'Bulk-solution numerical screening; not protein equilibration or converged aggregation thermodynamics.',
            'trajectory_sha256':file_hash(source/'bulk.dcd'),'topology_sha256':file_hash(source/'solution.prmtop'),
            'density_report_sha256':file_hash(source/'solution_validation.json'),
            'npt_ns':report['npt_ns'],'analyzed_frames':len(frames),'first_analyzed_ps':frames[0]['time_ps'],
            'last_analyzed_ps':frames[-1]['time_ps'],'periodic_boundary_method':'Minimum-image distances; orthorhombic box per frame.',
            'citrate_count':len(citrate),'citrate_heavy_atom_contact_cutoff_A':4.,
            'citrate_largest_cluster_mean':float(np.mean([r['citrate_largest_cluster'] for r in frames])),
            'citrate_largest_cluster_max':max(r['citrate_largest_cluster'] for r in frames),
            'citrate_pair_contact_fractions':{key:float(np.mean(value)) for key,value in pair_contacts.items()},
            'octanoic_acid_count':len(octanoic),'octanoic_acid_contact_frame_fraction':float(np.mean([r['octanoic_acid_contact'] for r in frames])) if len(octanoic)>=2 else None,
            'ethanol_octant_cv_mean':float(np.mean([r['ethanol_octant_cv'] for r in frames])) if ethanol else None,
            'sodium_near_citrate_fraction_mean':float(np.mean([r['sodium_near_citrate_fraction'] for r in frames])) if sodium.size and citrate_oxygens else None,
            'remaining_checks':['protein_protonation_and_equilibration'],
            'sampling_limit':f"One {report['npt_ns']:g} ns trajectory for this mixture. Contact persistence and positional microstate populations are not established."}
    atomic_json(output/'bulk_structure_timeseries.json',frames)
    atomic_json(output/'bulk_structure_analysis.json',result)
    return {'output':str(output),**result}
