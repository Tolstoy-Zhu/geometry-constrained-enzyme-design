# Third-Party Software And Data Notices

This file is an inventory, not a replacement for the complete license text shipped by each dependency. Installed package metadata and upstream repositories remain authoritative.

## Core Software

| Component | Pinned/observed version | License | Upstream |
|---|---:|---|---|
| Foldseek | `7f4b94644d6b210af4530e04a9fdc9e4bb64ea66` | GPL-3.0 | https://github.com/steineggerlab/foldseek |
| ColabFold | 1.6.1 | MIT | https://github.com/sokrypton/ColabFold |
| AlphaFold / alphafold-colabfold | 2.3.13 | Apache-2.0 code; parameters separately licensed | https://github.com/google-deepmind/alphafold |
| JAX | 0.5.3 | Apache-2.0 | https://github.com/jax-ml/jax |
| OpenMM | 8.5.2 | See upstream bundled notices | https://github.com/openmm/openmm |
| MMseqs2 | 18.8cc5c | GPL-3.0 | https://github.com/soedinglab/MMseqs2 |
| NCBI BLAST+ | 2.17.0 | Public domain / NCBI terms and bundled notices | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ |
| HHsuite | 3.3.0 | GPL-3.0 | https://github.com/soedinglab/hh-suite |
| tmtools | 0.3.0 | GPL-3.0 | https://github.com/jvkersch/tmtools |
| Streamlit | 1.59.1 | Apache-2.0 | https://github.com/streamlit/streamlit |
| Plotly.py | 6.9.0 | MIT | https://github.com/plotly/plotly.py |
| py3Dmol | 2.5.5 | MIT | https://github.com/3dmol/3Dmol.js |
| Pydantic | 2.13.4 | MIT | https://github.com/pydantic/pydantic |
| pandas | 3.0.3 | BSD-3-Clause | https://github.com/pandas-dev/pandas |
| NumPy | 2.4.6 | BSD-3-Clause and bundled component licenses | https://github.com/numpy/numpy |
| Biopython | 1.87 | Biopython License Agreement | https://github.com/biopython/biopython |
| requests | 2.34.2 | Apache-2.0 | https://github.com/psf/requests |
| PyYAML | 6.0.3 | MIT | https://github.com/yaml/pyyaml |
| openpyxl | 3.1.5 | MIT | https://foss.heptapod.net/openpyxl/openpyxl |
| tqdm | 4.68.4 | MPL-2.0 and MIT | https://github.com/tqdm/tqdm |
| pytest | 9.1.1 | MIT | https://github.com/pytest-dev/pytest |
| Playwright for Python | 1.61.0 | Apache-2.0 | https://github.com/microsoft/playwright-python |
| cloudflared | 2026.7.1 | Apache-2.0 | https://github.com/cloudflare/cloudflared |

Exact transitive dependency versions are recorded in `requirements.lock.txt`.
The isolated WSL prediction environment is recorded in `configs/prediction_runtime.lock.yaml`, `/opt/colabfold/pip-freeze.txt`, and `/opt/colabfold/conda-lock.json`.
The NR runtime source, version, official MD5, and locally verified SHA-256 are recorded in `configs/nr_runtime.lock.yaml` and `/opt/ncbi-blast/2.17.0/INSTALL_PROVENANCE.json`.
UI acceptance uses the locally installed Microsoft Edge browser under Microsoft's product terms; no browser binary is redistributed.

## Cloudflare Tunnel

- The optional online demonstration route uses the separately downloaded `cloudflared` executable. The binary is not committed to this repository or included in release archives.
- Version, download URL, byte length, SHA-256, and license metadata are recorded in `configs/cloudflared_runtime.lock.json`.
- A Quick Tunnel is a Cloudflare network service governed by Cloudflare's applicable online-services terms and has no uptime guarantee. It is used here only for temporary expert demonstrations.
- A durable public hostname requires an owner-managed Cloudflare account, named Tunnel, DNS hostname, and recommended Access policy. Tunnel credentials remain local and are not distributed with the project.

## AlphaFold Parameters And Public MSA

- AlphaFold source code is Apache-2.0. ColabFold 1.6.1 downloads the pTM parameters from `https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar`; those parameter files are licensed under CC BY-NC-SA 4.0 as stated in the bundled `params/LICENSE`. Parameter files and that 3.47 GB archive are not included in this repository or its release archive.
- ColabFold's public MMseqs2 MSA service is a third-party network service. A cache miss can transmit the submitted protein sequence. The application requires per-job confirmation and serializes requests from this machine.
- Predicted structures and confidence values remain computational outputs. Their use does not transfer or replace the licenses and citation requirements of AlphaFold, ColabFold, sequence databases, or templates.
- The installer applies a documented one-line ColabFold 1.6.1 network patch that changes only the result-archive read timeout from 6.02 to 180 seconds. `COLABFOLD_PATCH_PROVENANCE.json` records the original and patched source hashes; model algorithms and MSA content are unchanged.

## WSL And Ubuntu

- WSL is distributed by Microsoft under its applicable product terms: https://github.com/microsoft/WSL
- The imported Ubuntu 22.04 root filesystem contains many packages under their individual licenses. Package-level notices remain inside the distribution.

## Biological Databases And Structures

The project retrieves metadata or structures from RCSB PDB, UniProt, AlphaFold DB, and NCBI. Their data-use and attribution policies apply independently:

- RCSB PDB: https://www.rcsb.org/pages/policies
- UniProt: https://www.uniprot.org/help/license
- AlphaFold DB: https://alphafold.ebi.ac.uk/
- NCBI: https://www.ncbi.nlm.nih.gov/home/about/policies/
- NCBI BLAST API guidance: https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html

The NCBI remote NR route can transmit the submitted protein sequence. It requires a contact email and explicit confirmation for every cache miss, serializes requests, and applies NCBI-compatible polling. Query-derived local snapshots, complete NR databases, real candidates, caches, logs, and returned server archives are excluded from Git and public Release assets.

Downloaded databases and structures are not relicensed by this project. Source URLs, accessions, versions, and checksums are preserved where available.
