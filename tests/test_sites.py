from __future__ import annotations

import io

import pytest
from Bio.PDB import MMCIFIO

from bur01_structmine.sites import ligand_contacts, nonprotein_components, parse_structure_bytes, protein_sites, site_selections, viewer_selection


def atom(serial, atom_name, residue, chain, number, insertion, x, record="ATOM", element="C"):
    return f"{record:<6}{serial:5d} {atom_name:^4} {residue:>3} {chain}{number:4d}{insertion:1}   {x:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{20:6.2f}          {element:>2}\n"


def fixture():
    return ("MODEL        1\n" + atom(1, "CA", "ALA", "A", 10, " ", 0) + atom(2, "CA", "SER", "A", 10, "A", 1)
            + atom(3, "C1", "LIG", "B", 5, " ", 3, "HETATM") + atom(4, "O", "HOH", "B", 6, " ", 5, "HETATM", "O")
            + atom(5, "NA", "NA", "B", 7, " ", 8, "HETATM", "NA") + "ENDMDL\nMODEL        2\n"
            + atom(6, "CA", "ALA", "Z", 1, " ", 10) + "ENDMDL\nEND\n").encode()


def test_model_chain_insertion_code_and_distances():
    doc = parse_structure_bytes(fixture(), "test.pdb")
    assert doc.models == ["0", "1"]
    assert doc.chains("1") == ["Z"]
    ligands = nonprotein_components(doc, "0")
    assert ligands.default_selected.tolist() == [True, False, False]
    ligand_ids = [ligands.iloc[0].residue_id]
    contacts = ligand_contacts(doc, "0", ["A"], ligand_ids, 4)
    assert contacts.minimum_distance_angstrom.tolist() == [2, 3]
    assert contacts.insertion_code.tolist() == ["A", ""]
    key = contacts.iloc[0].residue_id
    assert viewer_selection(key)["icode"] == "A"
    marks = site_selections(doc, "0", {key: "focus"}, ligand_ids, 4, "example", "CCO")
    assert marks[0].residue_number == 10 and marks[0].insertion_code == "A"
    assert marks[0].minimum_distance_angstrom == 2
    assert len(ligand_contacts(doc, "0", ["A"], ligand_ids, 2.5)) == 1


def test_mmcif_and_long_chain_id_preserved():
    doc = parse_structure_bytes(fixture(), "test.pdb")
    doc.structure[0]["A"].id = "LONG"
    writer = MMCIFIO()
    writer.set_structure(doc.structure)
    buffer = io.StringIO()
    writer.save(buffer)
    cif = parse_structure_bytes(buffer.getvalue().encode(), "test.mmcif")
    assert "LONG" in cif.chains("0")
    assert protein_sites(cif, "0", ["LONG"]).insertion_code.tolist() == ["", "A"]


def test_no_ligand_no_fabricated_contact_and_invalid_cutoff():
    doc = parse_structure_bytes(fixture(), "test.pdb")
    assert ligand_contacts(doc, "0", ["A"], []).empty
    with pytest.raises(ValueError):
        ligand_contacts(doc, "0", ["A"], [], 0)
