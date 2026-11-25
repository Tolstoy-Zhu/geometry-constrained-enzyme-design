"""Recompute a widest solvent path on an atomic clearance grid, without literature centerlines."""
from __future__ import annotations
import heapq
import numpy as np
VDW={"C":1.70,"N":1.55,"O":1.52,"S":1.80,"P":1.80}

def clearance_grid(xyz,elements,spacing_A=.8,padding_A=5.,max_voxels=2_000_000,origin_fraction=0.):
    xyz=np.asarray(xyz,float);elements=np.asarray(elements)
    if xyz.ndim!=2 or xyz.shape[1]!=3 or len(elements)!=len(xyz) or not np.isfinite(xyz).all():
        raise ValueError('Invalid protein coordinates')
    if not .3<=spacing_A<=2 or not 0<=origin_fraction<1 or padding_A<=0:
        raise ValueError("Invalid channel grid settings")
    keep=np.isin(elements,list(VDW));xyz=xyz[keep];elements=elements[keep]
    if not len(xyz):raise ValueError("No recognized protein heavy atoms")
    lower=xyz.min(0)-padding_A-origin_fraction*spacing_A;upper=xyz.max(0)+padding_A
    shape=np.ceil((upper-lower)/spacing_A).astype(int)+1
    size=int(np.prod(shape))
    if size>max_voxels:raise ValueError("Channel grid exceeds memory budget")
    from scipy.spatial import cKDTree
    ijk=np.stack(np.unravel_index(np.arange(size),shape),axis=1)
    points=lower+ijk*spacing_A
    clearance=np.full(size,np.inf)
    # Separate element trees give exact minimum surface clearance for the chosen radii.
    for element,radius in VDW.items():
        atoms=xyz[elements==element]
        if len(atoms):clearance=np.minimum(clearance,cKDTree(atoms).query(points,workers=1)[0]-radius)
    return {'xyz':xyz,'elements':elements,'points':points,'ijk':ijk,'shape':shape,
            'clearance':clearance,'spacing_A':spacing_A,'origin_fraction':origin_fraction,
            'atom_tree':cKDTree(xyz),'radii':np.asarray([VDW[str(e)] for e in elements]),'edge_cache':{}}


def grid_path(grid,seed,probe_radius_A=1.2,continuous_edges=False):
    seed=np.asarray(seed,float)
    if seed.shape!=(3,) or not np.isfinite(seed).all() or not np.isfinite(probe_radius_A) or probe_radius_A<=0:
        raise ValueError('Invalid channel seed or probe')
    points=grid['points'];ijk=grid['ijk'];shape=grid['shape'];clearance=grid['clearance']
    size=len(points);spacing_A=grid['spacing_A']
    near=np.linalg.norm(points-seed,axis=1)
    starts=np.flatnonzero((clearance>=probe_radius_A)&(near<=3.0))
    if not len(starts):return {"status":"no_accessible_seed","bottleneck_A":None}
    start=int(starts[np.argmin(near[starts])])
    best=np.full(size,-np.inf);distance=np.full(size,np.inf);previous=np.full(size,-1,dtype=np.int64)
    best[start]=clearance[start];distance[start]=0
    queue=[(-best[start],0.,start)];end=None
    strides=np.array([shape[1]*shape[2],shape[2],1])
    while queue:
        neg,width_distance,node=heapq.heappop(queue);width=-neg
        if width<best[node]-1e-10 or width_distance>distance[node]+1e-10:continue
        loc=ijk[node]
        if np.any(loc==0) or np.any(loc==shape-1):end=node;break
        for axis in range(3):
            for sign in (-1,1):
                if not 0<=loc[axis]+sign<shape[axis]:continue
                nxt=int(node+sign*strides[axis])
                if clearance[nxt]<probe_radius_A:continue
                value=min(width,clearance[nxt]);length=width_distance+spacing_A
                if continuous_edges:
                    # Distance-to-surface is 1-Lipschitz. The endpoint lower bound
                    # certifies wide edges cheaply; only narrow edges need atoms.
                    if min(clearance[node],clearance[nxt])-spacing_A/2<value:
                        edge=tuple(sorted((node,nxt)))
                        if edge not in grid['edge_cache']:
                            a=points[node];b=points[nxt];upper=min(clearance[node],clearance[nxt])
                            nearby=grid['atom_tree'].query_ball_point((a+b)/2,upper+max(VDW.values())+spacing_A/2+1e-8)
                            atoms=grid['xyz'][nearby];delta=b-a
                            fraction=np.clip((atoms-a)@delta/(delta@delta),0.,1.)
                            edge_value=float((np.linalg.norm(atoms-a-fraction[:,None]*delta,axis=1)-grid['radii'][nearby]).min())
                            grid['edge_cache'][edge]=edge_value
                        value=min(value,grid['edge_cache'][edge])
                    if value<probe_radius_A:continue
                if value>best[nxt]+1e-10 or (abs(value-best[nxt])<=1e-10 and length<distance[nxt]):
                    best[nxt]=value;distance[nxt]=length;previous[nxt]=node
                    heapq.heappush(queue,(-value,length,nxt))
    if end is None:return {"status":"no_open_path_at_probe_radius","bottleneck_A":None}
    path=[];node=end
    while node>=0:path.append(node);node=int(previous[node])
    path=path[::-1]
    return {"status":"path_found","bottleneck_A":float(best[end]),"length_A":float(distance[end]),
            "centerline_A":points[path].tolist(),"seed_offset_A":float(near[start]),
            "spacing_A":spacing_A,"probe_radius_A":probe_radius_A,"method":"aqua_grid_continuous_edge_v2" if continuous_edges else "aqua_grid_widest_path_v1",
            "scope":"geometric_sampling_descriptor_not_transport_rate"}


def widest_path(xyz,elements,seed,spacing_A=.8,probe_radius_A=1.2,padding_A=5.,max_voxels=2_000_000):
    return grid_path(clearance_grid(xyz,elements,spacing_A,padding_A,max_voxels),seed,probe_radius_A,continuous_edges=True)


def continuous_path_clearance(points,xyz,elements):
    """Exact segment-to-sphere clearance, including locations between grid nodes."""
    points=np.asarray(points,float);xyz=np.asarray(xyz,float)
    if points.ndim!=2 or points.shape[1]!=3 or len(points)<2 or not np.isfinite(points).all():
        raise ValueError('A finite path with at least two nodes is required')
    radii=np.asarray([VDW[str(e)] for e in elements])
    values=[]
    for first,last in zip(points[:-1],points[1:]):
        delta=last-first;length2=float(delta@delta)
        fraction=np.clip((xyz-first)@delta/length2,0.,1.) if length2 else np.zeros(len(xyz))
        values.append(float((np.linalg.norm(xyz-first-fraction[:,None]*delta,axis=1)-radii).min()))
    return min(values)
