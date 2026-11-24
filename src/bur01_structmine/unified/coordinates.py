"""Verify that checked protein coordinates are the actual Amber simulation input."""
from __future__ import annotations
import numpy as np

def validate_amber_coordinates(prmtop_path,inpcrd_path,pdb_path,construct_sequences,chain_order=None):
    from openmm import app,unit
    topology=app.AmberPrmtopFile(str(prmtop_path)).topology
    coordinates=app.AmberInpcrdFile(str(inpcrd_path))
    protein=app.PDBFile(str(pdb_path))
    if coordinates.boxVectors is None:raise ValueError("Equilibrated input has no periodic box")
    xyz=np.asarray(coordinates.positions.value_in_unit(unit.angstrom))
    if len(xyz)!=topology.getNumAtoms() or not np.isfinite(xyz).all():
        raise ValueError("Amber topology and coordinates have different atom counts or invalid coordinates")
    from Bio.Data.PDBData import protein_letters_3to1
    names=dict(protein_letters_3to1)
    names.update({"HID":"H","HIE":"H","HIP":"H","ASH":"D","GLH":"E","CYM":"C","CYX":"C","LYN":"K"})
    amber=[r for r in topology.residues() if r.name in names]
    order=chain_order or list(construct_sequences)
    if set(order)!=set(construct_sequences) or len(order)!=len(construct_sequences):
        raise ValueError("Topology construct order must list each confirmed chain once")
    expected="".join(construct_sequences[chain] for chain in order)
    if "".join(names[r.name] for r in amber)!=expected:
        raise ValueError("Amber topology does not contain the confirmed complete protein constructs")
    chains={c.id:list(c.residues()) for c in protein.topology.chains()}
    pdb_res=[r for chain in order for r in chains.get(chain,[]) if r.name in names]
    if len(pdb_res)!=len(amber):raise ValueError("Protein PDB does not match the Amber protein residue count")
    pdb_xyz=np.asarray(protein.positions.value_in_unit(unit.angstrom))
    errors=[]
    for source,target in zip(amber,pdb_res):
        def heavy(residue):
            return {a.name:a.index for a in residue.atoms() if a.element is not None and a.element.symbol!="H"}
        a,b=heavy(source),heavy(target)
        if set(a)!=set(b):raise ValueError("Protein PDB and topology have different heavy atoms at residue "+source.id)
        errors.extend(np.linalg.norm(xyz[a[name]]-pdb_xyz[b[name]]) for name in a)
    maximum=float(max(errors,default=float("inf")))
    if maximum>.01:
        raise ValueError(f"Protein validation coordinates differ from the equilibrated system by {maximum:.4f} A")
    return {"status":"passed","protein_residues":len(amber),"atoms":len(xyz),"max_heavy_atom_difference_A":maximum}
