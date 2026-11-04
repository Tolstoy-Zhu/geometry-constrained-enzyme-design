"""Bur01-family, catalytic-site-constrained local sequence design."""
from __future__ import annotations

import hashlib
import itertools
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .config import project_root
from .design_models import ALPHABET, conditional_scores
from .structure_search import compare_ca_traces, load_ca_trace

TRIAD_ATOMS = (("OG",), ("OD1", "OD2"), ("ND1", "NE2"))
# Element indices are atomic numbers used by LigandMPNN.
ELEMENTS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
            "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm "
            "Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No "
            "Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").upper().split()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def ligand_groups(path: Path) -> dict:
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    model = next(parser.get_structure("parent", str(path)).get_models())
    groups = {}
    for chain in model:
        for r in chain:
            hetero, number, ins = r.id
            if not hetero.startswith("H_") or r.resname.upper() in {"HOH", "DOD", "WAT"}:
                continue
            key = f"{chain.id}:{r.resname}:{number}{ins.strip()}"
            groups[key] = [(a.element.upper(), a.coord) for a in r if a.element.upper() not in {"H", "D"}]
    return groups


def ligand_context(path: Path, selected):
    groups = ligand_groups(path) if selected else {}
    if set(selected) - groups.keys():
        raise ValueError(f"Unknown ligand. Available IDs: {', '.join(groups)}")
    atoms = [atom for key in dict.fromkeys(selected) for atom in groups[key]]
    if selected and not atoms:
        raise ValueError("Selected ligand has no heavy atoms.")
    if any(e not in ELEMENTS or e == "X" for e, _ in atoms):
        raise ValueError("Selected ligand contains an unidentified element.")
    if len(atoms) > 1000:
        raise ValueError("Select at most 1000 context heavy atoms.")
    xyz = np.asarray([c for _, c in atoms], dtype=np.float32).reshape(-1, 3)
    if not np.isfinite(xyz).all():
        raise ValueError("Ligand has non-finite coordinates.")
    return xyz, np.asarray([ELEMENTS.index(e) for e, _ in atoms], dtype=int)


def alignment_support(path: Path, parent: str):
    """Unique aligned homologs, excluding the parent; report both identity and coverage."""
    records = list(SeqIO.parse(str(path), "fasta"))
    rows = [re.sub(r"[a-z.]", "", str(r.seq)) for r in records]
    if not rows or rows[0].replace("-", "") != parent:
        raise ValueError("The first MSA record must match the parent coordinate sequence.")
    cols = [index for index, aa in enumerate(rows[0]) if aa != "-"]
    if any(len(s) != len(rows[0]) or set(s) - set(ALPHABET + "-XBZJUO") for s in rows):
        raise ValueError("Invalid aligned FASTA/A3M rows after removing insertions.")
    homologs = {tuple(s[i] for i in cols) for s in rows[1:]} - {tuple(parent)}
    if len(homologs) < 3:
        raise ValueError("Provide at least three distinct non-parent homologs for conservation screening.")
    letters = np.asarray(sorted(homologs))
    valid = np.isin(letters, list(ALPHABET)).sum(axis=0)
    same = (letters == np.array(list(parent))).sum(axis=0)
    support = np.divide(same, valid, out=np.zeros(len(parent)), where=valid > 0)
    return len(homologs), valid, support


def protein_trace(path, chain):
    trace = load_ca_trace(path, chain_id=chain)
    if not 20 <= len(trace.sequence) <= 1000 or set(trace.sequence) - set(ALPHABET):
        raise ValueError("Use one 20–1000 residue protein chain with standard amino acids.")
    for left, right in zip(trace.residues, trace.residues[1:]):
        step = right.residue_number - left.residue_number
        if step not in (0, 1) or np.linalg.norm(left.ca - right.ca) > 4.5:
            raise ValueError(f"Discontinuous coordinate chain at {left.label}/{right.label}; supply a complete parent construct.")
    for r in trace.residues:
        if not all(a in r.atoms for a in ("N", "CA", "C", "O")) or not np.isfinite(list(r.atoms.values())).all():
            raise ValueError(f"Missing backbone atom or non-finite coordinates at {r.label}.")
    return trace


def site_policy(trace, ligand_xyz, msa_path, fixed, site_limit, radius):
    ref = load_ca_trace(project_root() / "data/query/Bur01m.pdb", chain_id="A")
    comparison = compare_ca_traces(ref, trace)
    if min(comparison["qtmscore"], comparison["ttmscore"]) < 0.5 or comparison["mapped_catalytic_identity_count"] != 3:
        raise ValueError("The parent does not pass Bur01 fold and Ser/Asp/His identity checks.")
    labels = [r.label for r in trace.residues]
    triad = [labels.index(label) for label in comparison["mapped_catalytic_residues"].split(";")]
    core = []
    for index, names in zip(triad, TRIAD_ATOMS):
        atoms = trace.residues[index].atoms
        if any(a not in atoms for a in names):
            raise ValueError(f"Missing catalytic side-chain coordinates at {labels[index]}.")
        core.extend(atoms[a] for a in names)
    core = np.asarray(core)
    if ligand_xyz.size and np.linalg.norm(ligand_xyz[:, None] - core[None], axis=-1).min() > radius:
        raise ValueError("The selected ligand context is outside the catalytic region.")
    if set(fixed) - set(labels):
        raise ValueError("Unknown fixed parent residue: " + ", ".join(sorted(set(fixed) - set(labels))))
    motif = set()
    ser = triad[0]
    if ser >= 2 and re.fullmatch("G.S.G", trace.sequence[ser-2:ser+3]):
        motif.update(range(ser - 2, ser + 3))
    count, valid, support = alignment_support(msa_path, trace.sequence) if msa_path else (0, None, None)
    rows = []
    for i, r in enumerate(trace.residues):
        atoms = np.asarray([v for k, v in r.atoms.items() if not k.lstrip("0123456789").startswith(("H", "D"))])
        core_d = float(np.linalg.norm(atoms[:, None] - core[None], axis=-1).min())
        ligand_d = float(np.linalg.norm(atoms[:, None] - ligand_xyz[None], axis=-1).min()) if ligand_xyz.size else None
        anchor_d = ligand_d if ligand_d is not None else core_d
        reasons = []
        if i in triad:
            reasons.append("catalytic_triad")
        if i in motif:
            reasons.append("serine_motif")
        if core_d <= 4.0:
            reasons.append("catalytic_contact")
        if r.label in fixed:
            reasons.append("user_fixed")
        if count and valid[i] >= 3 and valid[i] / count >= 0.70 and support[i] >= 0.85:
            reasons.append("msa_conserved")
        if anchor_d > radius:
            reasons.append("outside_shell")
        rows.append(dict(position=i + 1, residue=r.label, core_distance_A=core_d, ligand_distance_A=ligand_d,
                         anchor_distance_A=anchor_d, homolog_count=count, msa_observed=int(valid[i]) if count else None,
                         msa_parent_fraction=float(support[i]) if count and valid[i] else None,
                         status=";".join(reasons) or "designable"))
    eligible = sorted((i for i, r in enumerate(rows) if r["status"] == "designable"),
                      key=lambda i: (rows[i]["anchor_distance_A"], i))
    for i in eligible[site_limit:]:
        rows[i]["status"] = "outside_site_budget"
    if not eligible:
        raise ValueError("No designable positions remain under the current constraints.")
    return rows, eligible[:site_limit], triad


def propose(trace, sites, scores, control, max_edits, count, priority):
    singles = []
    landscape = []
    for row, (index, values) in enumerate(zip(sites, scores)):
        wt = trace.sequence[index]
        wt_index = ALPHABET.index(wt)
        odds = values - values[wt_index]
        delta = odds - (control[row] - control[row, wt_index]) if control is not None else None
        choices = sorted((j for j in range(20) if j != wt_index), key=lambda j: -odds[j])
        for j in choices:
            landscape.append(dict(position=index + 1, residue=trace.residues[index].label, alternative=ALPHABET[j],
                                  log_odds=float(odds[j]), ligand_context_delta=float(delta[j]) if delta is not None else None))
        best = max(choices[:3], key=lambda j: delta[j]) if priority == "ligand_response" else choices[0]
        singles.append((index, ALPHABET[best], float(odds[best]), float(delta[best]) if delta is not None else None))
    # Rank within each mutation count. Interleaving prevents sums from automatically favouring one count.
    ranked = []
    for size in range(1, max_edits + 1):
        options = []
        for edits in itertools.combinations(singles, size):
            if any(np.linalg.norm(trace.coords[a[0]] - trace.coords[b[0]]) < 8.0
                   for a, b in itertools.combinations(edits, 2)):
                continue
            score = sum(e[2] for e in edits)
            response = sum(e[3] for e in edits) if control is not None else None
            options.append((score, response, sorted(edits)))
        key_index = 1 if priority == "ligand_response" else 0
        ranked.append(sorted(options, key=lambda item: (-item[key_index], -item[0], item[2]))[:count])
    candidates = []
    for group in itertools.zip_longest(*ranked):
        for option in group:
            if option is None:
                continue
            score, response, edits = option
            seq = list(trace.sequence)
            for index, aa, _, _ in edits:
                seq[index] = aa
            seq = "".join(seq)
            analysis = ProteinAnalysis(seq)
            digest = sha(seq.encode("ascii"))
            candidates.append(dict(candidate_id=f"AD_{digest[:16]}", sequence_sha256=digest, sequence=seq,
                                   mutations_pdb=",".join(trace.residues[i].label + aa for i, aa, _, _ in edits),
                                   mutations_sequence=",".join(f"{trace.sequence[i]}{i+1}{aa}" for i, aa, _, _ in edits),
                                   substitutions=len(edits), additive_log_odds=score, ligand_context_delta=response, gravy=analysis.gravy(),
                                   estimated_charge_pH7=analysis.charge_at_pH(7.0), structure_status="awaiting_coordinates"))
            if len(candidates) == count:
                return candidates, landscape
    return candidates, landscape


def run_design(parent, output, *, chain="A", ligands=(), fixed=(), msa=None, site_limit=16,
               radius=10.0, max_edits=2, count=12, seed=17, device="cpu", priority="compatibility",
               parent_id=None, expected_sequence=None):
    if not (1 <= site_limit <= 48 and 4 < radius <= 20 and 1 <= max_edits <= 3 and 1 <= count <= 100
            and 0 <= seed < 2**32):
        raise ValueError("Invalid design settings (sites 1–48, radius (4,20], edits 1–3, count 1–100, seed uint32).")
    if priority not in {"compatibility", "ligand_response"} or (priority == "ligand_response" and not ligands):
        raise ValueError("Ligand-response ranking requires explicitly selected ligand coordinates.")
    parent, output = Path(parent), Path(output)
    suffix = parent.suffix.lower()
    if suffix not in {".pdb", ".cif", ".mmcif", ".ent"}:
        raise ValueError("Expected PDB, CIF or mmCIF coordinates.")
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        run = Path(tmp) / "run"
        run.mkdir()
        structure = run / ("parent" + suffix)
        shutil.copyfile(parent, structure)
        msa_path = run / "homologs.a3m" if msa else None
        if msa_path:
            shutil.copyfile(msa, msa_path)
        trace = protein_trace(structure, chain)
        if expected_sequence and trace.sequence != expected_sequence.strip().upper():
            raise ValueError("The selected coordinate chain does not match the selected candidate sequence.")
        xyz, elements = ligand_context(structure, ligands)
        rows, sites, triad = site_policy(trace, xyz, msa_path, fixed, site_limit, radius)
        scores, control, model = conditional_scores(trace, sites, xyz, elements, seed=seed, device=device)
        candidates, landscape = propose(trace, sites, scores, control, max_edits, count, priority)
        for name, table in (("sites", rows), ("substitutions", landscape), ("candidates", candidates)):
            pd.DataFrame(table).to_csv(run / f"{name}.tsv", sep="\t", index=False)
        (run / "candidates.fasta").write_text("".join(f">{r['candidate_id']}\n{r['sequence']}\n" for r in candidates), encoding="ascii")
        manifest = dict(method="aqua_design_v1", parent_id=parent_id or parent.name, parent_file=structure.name,
                        parent_sha256=sha(structure.read_bytes()),
                        parent_sequence_sha256=sha(trace.sequence.encode()), chain=chain, ligands=list(ligands),
                        fixed=list(fixed), triad_positions=[i + 1 for i in triad], triad_labels=[trace.residues[i].label for i in triad],
                        msa_sha256=sha(msa_path.read_bytes()) if msa_path else None, radius_A=radius, site_limit=site_limit,
                        evaluated_sites=len(sites), max_edits=max_edits, requested_candidates=count, candidates=len(candidates),
                        priority=priority, ligand_control="masked_same_checkpoint" if ligands else None,
                        model=model, coordinate_scope="first_model_selected_chain", score_scope="parent_conditioned_inverse_folding")
        (run / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        run.rename(output)
    return output


def review_design(run, candidate_id, structure, chain="A", origin="unspecified", prediction_manifest=None):
    """Attribute coordinates to an exact candidate sequence, then compare with its parent."""
    run, structure = Path(run), Path(structure)
    if origin not in {"unspecified", "experimental", "predicted", "template_model"}:
        raise ValueError("Unknown coordinate origin.")
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    entries = pd.read_csv(run / "candidates.tsv", sep="\t", dtype=str)
    chosen = entries[entries.candidate_id == candidate_id]
    if len(chosen) != 1:
        raise ValueError("Unknown or ambiguous design candidate.")
    record = chosen.iloc[0]
    query_path = run / manifest["parent_file"]
    if sha(query_path.read_bytes()) != manifest["parent_sha256"]:
        raise ValueError("Archived parent structure checksum mismatch.")
    query = protein_trace(query_path, manifest["chain"])
    target = protein_trace(structure, chain)
    if target.sequence != record.sequence or sha(target.sequence.encode()) != record.sequence_sha256:
        raise ValueError("Coordinates do not match this design sequence.")
    pair = compare_ca_traces(query, target, catalytic_query_residues=())
    # Designs have substitutions only. Direct sequence indices also support insertion-coded parent residues.
    core = []
    for pos, names in zip(manifest["triad_positions"], TRIAD_ATOMS):
        residue = target.residues[pos - 1]
        core.append(np.array([residue.atoms[a] for a in names]) if all(a in residue.atoms for a in names) else None)
    def distance(i, j):
        if core[i] is None or core[j] is None:
            return None
        return float(np.linalg.norm(core[i][:, None] - core[j][None], axis=-1).min())
    report = dict(candidate_id=candidate_id, sequence_sha256=record.sequence_sha256, structure_sha256=sha(structure.read_bytes()),
                  coordinate_origin=origin, chain=chain, qTM=pair["qtmscore"], tTM=pair["ttmscore"], RMSD_A=pair["rmsd"],
                  aligned_residues=pair["alnlen"], Ser_His_A=distance(0, 2), Asp_His_A=distance(1, 2))
    if prediction_manifest is not None:
        pm = Path(prediction_manifest)
        pred = json.loads(pm.read_text(encoding="utf-8"))
        if (origin != "predicted" or pred.get("status") != "completed" or pred.get("sequence_sha256") != record.sequence_sha256
                or sha((pm.parent / "rank_001.pdb").read_bytes()) != report["structure_sha256"]):
            raise ValueError("AF2 manifest or rank-one coordinates do not match this design.")
        report.update(af2_prediction_id=pred["prediction_id"], mean_plddt=pred.get("mean_plddt"), pTM=pred.get("ptm"))
        confidence = pm.parent / "rank_001_plddt.json"
        if confidence.is_file():
            values = np.asarray(json.loads(confidence.read_text(encoding="utf-8")).get("plddt", []), dtype=float)
            if values.shape == (len(target.residues),) and np.all((values >= 0) & (values <= 100)):
                report["triad_plddt"] = [float(values[p - 1]) for p in manifest["triad_positions"]]
    data = json.dumps(report, indent=2).encode("utf-8")
    directory = run / "reviews"
    directory.mkdir(exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=directory, delete=False) as handle:
        tmp = Path(handle.name)
        handle.write(data)
    try:
        tmp.replace(directory / f"{sha(data)}.json")
    finally:
        tmp.unlink(missing_ok=True)
    return report
