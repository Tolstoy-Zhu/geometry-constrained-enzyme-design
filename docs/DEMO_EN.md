# Read-only Bur01 Demonstration

## Start and Stop

Double-click `start_demo.bat`. It uses the existing project Python environment, selects an available loopback port and opens a browser. After presentation, double-click `stop_demo.bat`; closing the browser alone does not stop the server.

For real experimental imports or computation, use `start_ui.bat` and stop it separately with `stop_ui.bat`. This update installs no logon task, watchdog recovery or public tunnel. The demo listens only on 127.0.0.1; that URL is not accessible from another computer. Present locally on a configured computer or share the screen.

## Five Steps

1. **Input structure**: public experimental PDB 8HFW chain A, residues 12–294 and Ser111/Asp241/His274, with rotatable coordinates.
2. **Candidate filtering**: verified real Foldseek benchmark records, positive/negative controls and exclusion reasons. This is a recorded computation, not a new live database search.
3. **Structure comparison**: local TM-align overlay, qTM/tTM/RMSD and residue alignment. Structure similarity does not establish enzyme function.
4. **Score explanation**: four rubric components, gaps, geometry status and version. Predictive error remains unestimated.
5. **Experimental assessment**: explicitly simulated measurements demonstrate replicate summaries, error metrics, ranking and Top-k. They are not Bur01 assay results.

Filtering, rotation and inspection are allowed; task submission, reannotation and production writes have no demo entry point. `data/demo/demo_manifest.json` verifies SHA-256 for every included computational asset. Missing or changed assets raise an error.

## Offline Scope

3Dmol 2.5.5, structures, computation summaries and simulated experiments are local. Starting or navigating the demo does not invoke WSL, GPU inference, external prediction or download services. Browser WebGL rendering remains necessary and is distinct from AlphaFold GPU inference. A local Python/Streamlit process is required; this is not standalone HTML transferable as a localhost URL.

## Sites and Ligands

In the full workbench, expand **Substrates and sites** under **Structure**. Upload PDB/CIF/mmCIF, choose a model, protein chains and existing non-protein components. Water and common ions are not selected by default. Other components are ligand/cofactor options, not automatically confirmed reaction substrates.

The default 4-Angstrom heavy-atom contact table supports row selection and 3D highlights. Residues can also be selected directly. Exports preserve model, chain, residue number, insertion code, structure hash, cutoff, ligand identifiers and keep/focus/exclude marks. SMILES is metadata, without chemical-validity or binding prediction claims. No ligand coordinates means metadata-only status, never an invented binding pose.

```powershell
python -m bur01_structmine.cli demo
python -m bur01_structmine.cli demo --no-browser
python -m bur01_structmine.cli stop-demo
```

Service state/logs are confined to `.cache/demo_service/`. Stop verifies PID, random token and project script before terminating its own process tree. Omit `--port` to select an available port automatically.
