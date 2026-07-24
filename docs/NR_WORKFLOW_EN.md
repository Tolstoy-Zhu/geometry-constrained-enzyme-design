# AquaEsterNet-StructMine v0.4 NR Workflow

## Purpose And Evidence Boundary

Many NCBI NR proteins have sequence records but no experimental or predicted structure. v0.4 therefore adds an independent sequence-expansion layer after structure mining. It asks which natural sequences are related to Bur01 or structure-supported seeds, which representatives conserve the mapped catalytic identities, and which sequences merit AF2 validation.

Sequence evidence is not structural evidence. Every NR-only candidate remains:

```text
wase_tier = Not evaluated - structure required
wase_scope = not_evaluated_sequence_only
```

Only AF2 followed by Foldseek/TM-align can add structural evidence. WASE-screen remains transparent triage and cannot establish enzyme activity, S/H ratio, or a water-phase esterification mechanism.

## Fixed Reproduction Preset

The `hong2025` preset follows the public VenusMine MMseqs2 search and clustering steps:

```text
iterations               2
max-seqs                 2500
search-type              1
cluster min sequence id  0.50
```

v0.4 does not add an untrained PLM composite score. Ranking is limited to E-value, bits, identity, query/target coverage, length, GXSXG-like motif, mapped Ser/Asp/His identity, and annotation state.

References:

- VenusMine search: https://github.com/ai4protein/VenusMine/tree/main/Step3%20MMseqs%20Search
- VenusMine clustering: https://github.com/ai4protein/VenusMine/tree/main/Step4%20MMseqs%20Cluster
- VenusMine paper: https://doi.org/10.1038/s41467-025-61599-z
- CPDiffusion: https://doi.org/10.1038/s41421-024-00728-2
- MMseqs2: https://github.com/soedinglab/MMseqs2
- NCBI BLAST API guidance: https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html

## Local Demonstration Snapshot

The workstation does not download complete NR. A single confirmed NCBI remote `blastp nr` request against Bur01 retains at most 2,500 hits and retrieves their sequences:

```powershell
$env:AQUA_NCBI_EMAIL = "name@example.org"
python -m bur01_structmine.cli nr-demo-build `
  --query-fasta data/query/Bur01m_seed.fasta `
  --max-hits 2500 `
  --confirm-ncbi
```

A valid contact email and explicit per-request consent are mandatory. Requests are serial, polling respects a 60-second minimum interval, and network failures use bounded retries. The snapshot manifest records the query hash, RID, parameters, timestamps, file hashes, and the explicit scope `nr_remote_snapshot`.

Run expansion from a completed structure run:

```powershell
python -m bur01_structmine.cli nr-expand `
  --source-run <structure_run_id> `
  --config configs/nr_local_demo.yaml
```

## Complete-NR Server Job

Export a credential-free server package:

```powershell
python -m bur01_structmine.cli nr-export-job `
  --source-run <structure_run_id> `
  --config configs/nr_full_server.yaml
```

The package contains generic Bash, Slurm, and PBS launchers, fixed seeds, configuration and metadata snapshots, resource checks, post-processing, checksums, and bilingual instructions. It contains no SSH host, password, token, or private key.

Database download is an explicit management action:

```bash
bash prepare_nr_database.sh
```

This calls `mmseqs databases NR`. Normal searches never download or update NR. The server job requires `NR_SNAPSHOT_DATE`, writes `server_result_manifest.json` and `RESULT_SHA256SUMS`, and creates a return ZIP.

Import:

```powershell
python -m bur01_structmine.cli nr-import-job `
  --bundle <nr_run_id>_server_result.zip
```

The importer rejects path traversal, missing files, checksum mismatch, a non-`full_nr` scope, incompatible project version, and an inconsistent run ID. A formal v0.4.0 Release is blocked until a complete-NR result has been imported and accepted.

## Processing Logic

1. Export Bur01 plus all available Foldseek target sequences from a structure run.
2. Select the best hit per accession, then deduplicate complete sequences by SHA-256.
3. Cluster at 50% sequence identity and retain membership provenance.
4. Directly realign every representative to 8HFW chain A.
5. Record mapped position, observed identity, and missing reason separately for Ser111, Asp241, and His274.
6. Build the AF2-ready queue with transparent thresholds.

P00520 must remain on hold under the default sequence screen because its significance, coverage, or catalytic evidence is insufficient.

## AF2 And Mining Lineage

```powershell
python -m bur01_structmine.cli predict `
  --nr-run <nr_run_id> `
  --accession <accession> `
  --preset standard `
  --confirm-public-msa
```

The prediction manifest records `source_nr_run_id` and accession; the NR manifest records `child_prediction_ids`. `predict-and-mine` extends the lineage from NR to AF2 rank_001 and then to Foldseek/TM-align and WASE-screen.

## Novelty Audit

```powershell
python -m bur01_structmine.cli nr-novelty `
  --fasta designed_sequence.fasta `
  --backend ncbi_remote `
  --email name@example.org `
  --confirm-ncbi
```

The output reports nearest accession, identity, coverage, E-value, title/organism, and snapshot provenance. A natural NR sequence matching itself at 100% is expected and is never counted as novelty evidence.

## Outputs And Status

Each run writes `results/nr_runs/<nr_run_id>/` with seeds, raw and deduplicated hits, clusters, direct Bur01 alignments, sequence-screen TSV/XLSX/FASTA, AF2 queue, events, manifest, warnings, and checksums.

As of 2026-07-24:

- MMseqs2 `18.8cc5c` and NCBI BLAST+ `2.17.0` pass runtime diagnostics.
- 65 isolated tests pass, including fake-MMseqs2 expansion/novelty execution, hardened archive import, and verified server lineage.
- A real NCBI snapshot still requires the user's explicit transfer confirmation and contact email.
- Complete NR still requires an external server/HPC result; the development branch must not be merged or tagged before that gate passes.

`configs/ncbi_prostt5.yaml` remains a legacy small-custom-FASTA route and is not complete NR.
