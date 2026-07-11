# AquaEsterNet-StructMine Model Guide

## 1. What It Is

AquaEsterNet-StructMine is not a single neural network. It is a reproducible structure-mining and evidence-integration workflow constrained by the Bur01 structure. It answers:

> Which proteins are sufficiently similar to Bur01 in fold, coverage, catalytic core, and annotation evidence to justify experiments?

The names have distinct scopes:

- **AquaEsterNet**: the broader framework for water-phase ester mechanism and candidate discovery.
- **StructMine**: the module for structure search, chain-level alignment, annotation, triage, and export.
- **WASE-screen**: the current transparent screen, not the future formal WASE model.

## 2. Inputs And Outputs

### Input

The core query is PDB 8HFW chain A:

```text
File: data/query/Bur01m.pdb
Residues: 12-294
Length: 283 aa
Catalytic triad: Ser111 / Asp241 / His274
```

Search targets can be Foldseek databases, small PDB/mmCIF folders, custom predicted-structure databases, or sequence sets represented through ProstT5.

### Output

Each candidate preserves structural metrics, query and target chains, aligned sequences, catalytic-residue mappings, database identity, reviewed status, AlphaFold confidence, final sequence, WASE-screen evidence, configuration, input hashes, software versions, database fingerprints, warnings, and failures.

## 3. Core Principle

```text
8HFW chain A
      |
      v
Foldseek database search / TM-align small-set alignment
      |
      v
TM-score + RMSD + bidirectional coverage filters
      |
      v
Ser111/Asp241/His274 mapping and geometry checks
      |
      v
RCSB / UniProt / AlphaFold / NCBI annotation
      |
      v
Transparent WASE-screen triage
      |
      v
Run manifest + Excel/TSV/FASTA + web visualization
```

### 3.1 Foldseek Search

Foldseek encodes protein structures with searchable 3Di representations and uses efficient prefiltering and structural alignment for database-scale discovery. Production results come from real Foldseek and are not replaced by the Python preview.

The workflow retains Foldseek E-value, bit score, probability, TM-score, RMSD, coverage, and identity fields and records the exact Foldseek and database versions in each run manifest.

### 3.2 TM-align Preview

Small local sets use `tmtools==0.3.0`, which calls the real TM-align algorithm. It reports query-normalized TM-score, target-normalized TM-score, RMSD, and aligned sequences. Local rows do not fabricate Foldseek E-value, bit score, or probability.

### 3.3 Chain-Level Processing

Every protein chain is handled independently:

- The query chain must be explicit.
- Each target protein chain is aligned separately.
- Chains are never concatenated into a synthetic protein.
- `query_chain` and `target_chain` are exported.

### 3.4 Catalytic-Core Mapping

Ser111, Asp241, and His274 are projected through the structural alignment. The workflow records mapped identities, residue numbers, Ser-His and Asp-His side-chain distances, and geometry contacts.

This is catalytic-core evidence only. It does not establish identical substrate tunnels, dynamics, kinetics, or selectivity.

### 3.5 Annotation Integration

Target identifiers are parsed into PDB, UniProt, or NCBI accessions. Configured clients retrieve experimental-structure metadata, reviewed status, function, cross-references, AlphaFold model versions and pLDDT, organisms, and sequences.

Requests use batching, retries, exponential backoff, disk caching, and explicit error columns. The selected sequence provenance is stored in `final_sequence_source`.

### 3.6 AlphaFold2-PTM Prediction Layer

When only a sequence is available, v0.3 uses LocalColabFold 1.6.1 and AlphaFold2-PTM. Prediction is not folded into WASE-screen as an undocumented feature. It creates a separately versioned query structure with its own evidence:

```text
monomer FASTA
  -> format, length, and X-content validation
  -> public MMseqs2 MSA cached by sequence hash
  -> Evoformer sequence/MSA processing and recycling
  -> Structure Module coordinates
  -> pLDDT, pTM, and PAE confidence outputs
  -> rank_001 Amber relaxation
  -> dynamic Foldseek/TM-align query
```

- `pLDDT` is local per-residue confidence, not activity.
- `pTM` reports global topology confidence, not catalytic correctness.
- `PAE` exposes uncertainty in relative residue and domain placement.
- The standard five models use distinct AF2-PTM parameters and are ranked by pLDDT.

Global sequence alignment maps 8HFW `Ser111/Asp241/His274` to the new sequence. Manual corrections and residue identities are recorded. Identity mismatches never receive fabricated catalytic credit.

Public MSA retrieval and local GPU inference are separate stages. MSA requests are serialized and retried up to three times with default 15- and 30-second exponential delays; original service errors are retained in events and manifests. A process lock allows one GPU job, cancellation terminates the corresponding WSL process, stale locks recover, and CUDA memory failures never fall back to CPU. JAX explicitly disables ColabFold unified memory for the 16 GB RTX 4080; this changes allocation behavior, not model algorithms. CPU OpenMM relaxes rank_001 after ranking.

## 4. WASE-screen v0.1

### Structural Fit, 45 Points

```text
25 * min(qTM, tTM)
+ 10 * sqrt(qcov * tcov)
+ 10 * clamp(1 - RMSD / 5, 0, 1)
```

Bidirectional TM-score limits length-normalization bias, coverage rejects local-fragment matches, and the RMSD term rewards tighter superpositions.

### Catalytic Core, 30 Points

- GXSXG-like motif: 5 points.
- Mapped Ser/Asp/His identities: 15 points.
- Catalytic geometry contacts: 10 points.

### Annotation Confidence, 15 Points

- Esterase, lipase, or hydrolase-related annotation.
- Swiss-Prot reviewed or experimental-structure evidence.
- AlphaFold confidence.

### Experimental Readiness, 10 Points

- Available sequence without unknown X residues.
- Plausible protein length.
- Adequate target coverage.
- Traceable accession.

### Tier Gate

High priority requires:

```text
wase_screen_score >= 70
wase_evidence_coverage >= 75
min(qTM, tTM) >= 0.50
catalytic_core_score >= 20
```

High partial scores cannot override inadequate evidence coverage.

### Missing Formal WASE Evidence

Formal WASE still requires substrate-tunnel specialization, alcohol retention, water exclusion, product escape, substrate compatibility, fermentation-source priors, adaptive-site enrichment, and ancestral gains. These remain explicit in `wase_missing_evidence`.

## 5. Innovation

1. **Structure-constrained acceptance instead of a zero-threshold demo**: fixed positives and negatives require P00520 to be rejected.
2. **Bidirectional TM-score and chain-level alignment**: reduces length and multi-chain artifacts.
3. **Catalytic mapping in the candidate record**: integrates overall fold with Ser/Asp/His identity and geometry.
4. **Traceable cross-database evidence**: cache state, errors, input hashes, software, and database versions are recorded.
5. **Shared CLI and web service layer**: both interfaces execute the same scientific logic.
6. **Explicit WASE-screen boundary**: uncertainty and missing evidence are reported without claiming S/H-ratio prediction.
7. **Bidirectional prediction-mining lineage**: rank_001, prediction ID, run ID, input hashes, and parameter checksums remain connected.
8. **Dynamic catalytic references**: new sequences map the 8HFW triad instead of assuming every query uses residue numbers 111/241/274.
9. **Privacy and GPU state are first-class evidence**: consent, cache state, serialized MSA access, locks, cancellation, and memory failures are recorded.

## 6. Usage

Web workbench:

```text
Double-click start_ui.bat
```

The workbench uses the stable address `http://127.0.0.1:8501`. Double-click `install_ui_autostart.bat` once to install a per-user Windows logon task and background watchdog. Streamlit is restarted after an unexpected exit, so no console window needs to remain open. Use `ui_status.bat`, `stop_ui.bat`, and `uninstall_ui_autostart.bat` for status, stop, and scheduled-task removal.

The AF Prediction view supports FASTA paste/upload, existing candidates, presets, templates, per-job privacy consent, catalytic-site preview and correction, predict-only, predict-and-mine, cancellation, pLDDT, PAE, and a confidence-colored structure.

For an expert demonstration, double-click `install_online_ui.bat`. Models and databases remain local while remote browsers use the same service layer through an authenticated HTTPS tunnel. See `ONLINE_DEMO_GUIDE_EN.md` for operation and security limits.

AlphaFold2 prediction:

```powershell
python -m bur01_structmine.cli predict --fasta query.fasta --preset standard --confirm-public-msa

python -m bur01_structmine.cli predict --candidate-run <run_id> --accession Q9FVW3 --preset standard --confirm-public-msa

python -m bur01_structmine.cli predict-and-mine --fasta query.fasta --preset standard --confirm-public-msa --mining-config configs/alphafold_swissprot.yaml
```

The standard preset is fixed at five AF2-PTM models, three recycles, one seed, and rank_001 Amber relaxation. Templates are enabled only with `--templates`.

Environment diagnosis:

```powershell
python -m bur01_structmine.cli doctor --config configs/alphafold_swissprot.yaml
```

Formal benchmark:

```powershell
python -m bur01_structmine.cli benchmark --engine foldseek
```

Small-set preview:

```powershell
python -m bur01_structmine.cli benchmark --engine tmtools
```

Production search:

```powershell
python prepare_databases.py --config configs/alphafold_swissprot.yaml --only alphafold_swissprot
python -m bur01_structmine.cli production-run --config configs/alphafold_swissprot.yaml
```

Validated reference: AlphaFold/Swiss-Prot `v6 2025_03` contains 590,183 structures. Standard 8HFW prediction produced five models; rank_001 has mean pLDDT `98.07`, pTM `0.950`, and `min(qTM,tTM)=0.9925` against experiment. Standard Q9FVW3 prediction has mean pLDDT `83.54`, pTM `0.740`, and `min(qTM,tTM)=0.7815` against AlphaFold DB v6. The predicted-8HFW Foldseek benchmark passed every condition, and the full search produced 999 raw hits, 12 filtered candidates, and 12 valid annotations with zero warnings.

## 7. Interpreting Candidates

1. Inspect `min(qtmscore,ttmscore)`, RMSD, and bidirectional coverage.
2. Confirm complete catalytic-triad mapping and plausible geometry.
3. Check reviewed or experimental evidence and AlphaFold confidence.
4. Read `wase_evidence_coverage` and `wase_missing_evidence`, not only the total score.
5. Confirm sequence completeness, plausible length, and accession provenance.
6. Before experiments, inspect the overlay, PAE, tunnels, and substrate poses manually.

## 8. Boundaries

- Structural similarity does not directly imply esterification activity.
- GXSXG-like motifs can occur in unrelated proteins.
- Low-confidence AlphaFold regions can distort coverage and local geometry.
- Annotation networks can fail; inspect `annotation_status` and `annotation_errors`.
- WASE-screen is a prioritization tool, not an experimental conclusion.
