"""Periodic-box occupancy and censored residence events, retaining each replica."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

def minimum_image(delta,box):
    box=np.asarray(box,float)
    if box.shape!=(3,3) or np.linalg.det(box)<=0:raise ValueError("Invalid periodic box")
    fraction=np.asarray(delta)@np.linalg.inv(box)
    return (fraction-np.round(fraction))@box

def residence_events(present,times):
    present=np.asarray(present,bool);times=np.asarray(times,float)
    if len(present)!=len(times) or len(times)<2 or np.any(np.diff(times)<=0):
        raise ValueError("Ordered frames with positive time increments are required")
    events=[];start=None
    for i,occupied in enumerate(present):
        if occupied and start is None:start=i
        if start is not None and (not occupied or i==len(present)-1):
            end=i if not occupied else i+1
            # Lower-bound observed duration, not an invented exit time.
            duration=float(times[end-1]-times[start])
            events.append({"start_ps":float(times[start]),"last_observed_ps":float(times[end-1]),
                           "observed_duration_ps":duration,"left_censored":start==0,
                           "right_censored":end==len(present)})
            start=None
    return events

def channel_bottleneck(points,protein_xyz,radii):
    points=np.asarray(points,float);protein_xyz=np.asarray(protein_xyz,float);radii=np.asarray(radii,float)
    if points.ndim!=2 or points.shape[1]!=3 or len(points)<2:raise ValueError("Channel centerline is required")
    clearance=np.linalg.norm(points[:,None]-protein_xyz[None],axis=-1)-radii[None]
    profile=clearance.min(1)
    return {"bottleneck_A":float(profile.min()),"clearance_A":profile.tolist()}

def analyze_trajectory(topology,trajectories,output,site_selection,solvent_selections,
                       radius_A=5.,replica_id="replica",channel_manifest=None,channel_stride=100):
    import MDAnalysis as mda
    from MDAnalysis.transformations import unwrap
    from MDAnalysis.lib.distances import distance_array
    from .datasets import file_hash
    universe=mda.Universe(str(topology),[str(p) for p in trajectories])
    protein=universe.select_atoms("protein")
    site=universe.select_atoms(site_selection)
    if not len(protein) or not len(site) or not len(protein.bonds):
        raise ValueError("Protein, catalytic region and bonded topology are required for PBC repair")
    universe.trajectory.add_transformations(unwrap(protein))
    groups={}
    for name,selection in solvent_selections.items():
        atoms=universe.select_atoms(selection)
        if not len(atoms):raise ValueError("No atoms match solvent selection: "+name)
        # One marker atom per molecule; prevent atom counts masquerading as molecule occupancy.
        if len(set(atoms.resindices))!=len(atoms):raise ValueError("Select exactly one marker per solvent molecule")
        groups[name]=atoms
    times=[];series={k:[] for k in groups};volumes=[];channels=[]
    if channel_stride<1:raise ValueError('Positive channel sampling stride is required')
    from .channels import widest_path
    try: elements=list(protein.elements)
    except mda.exceptions.NoDataError:
        elements=[a.name.lstrip('0123456789')[0] for a in protein]
    for ts in universe.trajectory:
        if ts.dimensions is None or not np.isfinite(ts.positions).all():raise ValueError("Missing box or invalid trajectory coordinates")
        if times and ts.time<=times[-1]:raise ValueError("Trajectory times overlap; supply ordered non-overlapping segments")
        times.append(float(ts.time));volumes.append(float(ts.volume))
        # Place all bonded protein components in the same periodic image as the first component.
        fragments=protein.fragments;anchor=fragments[0].center_of_geometry()
        for fragment in fragments[1:]:
            shift=minimum_image(fragment.center_of_geometry()-anchor,ts.triclinic_dimensions)-(fragment.center_of_geometry()-anchor)
            fragment.positions+=shift
        if (len(times)-1)%channel_stride==0 and not channel_manifest:
            channels.append({'frame':len(times)-1,'time_ps':float(ts.time),
                             **widest_path(protein.positions,elements,site.center_of_geometry())})
        for key,atoms in groups.items():
            distances=distance_array(atoms.positions,site.positions,box=ts.dimensions)
            series[key].append((distances.min(1)<=radius_A).tolist())
    result={"schema_version":"1.0.0","replica_id":replica_id,"frames":len(times),"radius_A":radius_A,
            "pbc":"bonded_protein_unwrap_and_minimum_image_distances","observables":{},
            "msm_tpt":{"status":"insufficient_sampling_validation","rates":None,
                       "required":["independent_replicas_with_transitions","ITS_plateau","Chapman_Kolmogorov"]}}
    if len(times)<2:raise ValueError("At least two trajectory frames are required")
    for key,rows in series.items():
        occupied=np.asarray(rows,bool);events=[]
        for molecule in range(occupied.shape[1]):
            for event in residence_events(occupied[:,molecule],times):
                events.append({"molecule_index":int(groups[key].resindices[molecule]),**event})
        concentrations=len(groups[key])/(np.asarray(volumes)*1e-27*6.02214076e23)
        result["observables"][key]={"mean_molecules_in_region":float(occupied.sum(1).mean()),
                                   "fraction_frames_occupied":float(occupied.any(1).mean()),
                                   "entry_events":sum(not e["left_censored"] for e in events),
                                   "residence_events":events,"molecules_in_box":len(groups[key]),
                                   "realized_concentration_M_mean":float(concentrations.mean())}
    result['channel']={'status':'recomputed_grid_paths','method':'aqua_grid_widest_path_v1','frame_results':channels,
                       'resolution_A':.8,'probe_radius_A':1.2,'stride_frames':channel_stride,
                       'scope':'geometric_descriptor_requires_grid_and_probe_sensitivity_checks'}
    if channel_manifest:
        channel=json.loads(Path(channel_manifest).read_text())
        expected=[file_hash(p) for p in trajectories]
        if channel.get("trajectory_sha256")!=expected or channel.get("topology_sha256")!=file_hash(topology):
            raise ValueError("Channel result was not generated from these trajectories and assembly")
        if not channel.get("tool_version") or not channel.get("frame_results"):
            raise ValueError("Channel result requires tool version and per-frame measurements")
        for row in channel["frame_results"]:
            if not np.isfinite(row["bottleneck_A"]) or not 0<=row["frame"]<len(times):raise ValueError("Invalid channel result")
        result["channel"]={"status":"recomputed","provenance":channel,
                           "mean_bottleneck_A":float(np.mean([r["bottleneck_A"] for r in channel["frame_results"]]))}
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/"occupancy.npz",time_ps=np.asarray(times),**{k:np.asarray(v) for k,v in series.items()})
    (output/"analysis.json").write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    return result
