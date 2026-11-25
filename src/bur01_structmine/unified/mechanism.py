"""Prepared Amber/OpenMM adapter with explicit scientific and runtime gates."""
from __future__ import annotations
import hashlib,json,os,platform,subprocess,sys
from importlib.metadata import version as distribution_version
from pathlib import Path
from .datasets import file_hash
from .mechanism_plan import MechanismBatch
from .mechanism_analysis import analyze_trajectory

OPENMM_VERSION="8.2.0"

def atomic_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")
    tmp.replace(path)

def checked_file(entry,base):
    from .transfer import resolve_reference
    path=resolve_reference(entry["path"],base)
    if not path.is_file():raise ValueError("Missing prepared file: "+str(path))
    if file_hash(path)!=entry["sha256"]:raise ValueError("Prepared file checksum mismatch: "+str(path))
    return path.resolve()

def inspect_prepared(prepared_manifest):
    """No legacy 'PASS' flag alone is sufficient to release a production run."""
    path=Path(prepared_manifest)
    data=json.loads(path.read_text(encoding="utf-8"))
    missing=[];files={}
    if data.get("schema_version")!="1.0.0":missing.append("unsupported_prepared_schema")
    if data.get("enzyme_state")!="free":missing.append("only_free_enzyme_supported")
    if data.get("forcefields")!={"protein":"ff19SB","water":"OPC","ligand":"GAFF2"}:
        missing.append("locked_forcefield_route_not_matched")
    if not data.get("protonation_assumptions"):missing.append("protonation_assumptions")
    missing.extend("unresolved:"+str(x) for x in data.get("unresolved",[]))
    for name in ("prmtop","inpcrd","protein_pdb","reaction"):
        try:files[name]=str(checked_file(data[name],path.parent))
        except (KeyError,ValueError):missing.append(name+"_missing_or_checksum_mismatch")
    gates=("construct","assembly","parameters","geometry","equilibration")
    for name in gates:
        gate=data.get("validation",{}).get(name,{})
        if gate.get("status")!="passed":missing.append(name+"_validation")
        try:
            evidence=json.loads(checked_file(gate["evidence"],path.parent).read_text(encoding="utf-8"))
            if evidence.get("status")!="passed":missing.append(name+"_evidence_not_passed")
            if evidence.get("unresolved") or evidence.get("blocking_issues"):missing.append(name+"_unresolved_issues")
        except (KeyError,ValueError,OSError):missing.append(name+"_evidence_missing")
    sequences=data.get("construct_sequences",{})
    if not sequences:missing.append("confirmed_full_construct_sequences")
    if "protein_pdb" in files and sequences:
        from ..structure_search import load_ca_trace
        from .families import assess_structure
        for chain,sequence in sequences.items():
            try:
                trace=load_ca_trace(Path(files["protein_pdb"]),chain_id=chain)
                if trace.sequence!=sequence:missing.append("construct_sequence_mismatch:"+chain)
                quality=assess_structure(trace,data.get("family_id","bur01"))
                if not quality["catalytic_identity_pass"] or any(s["missing_atoms"] for s in quality["sites"]):
                    missing.append("catalytic_atoms:"+chain)
                if any(g["status"]!="compatible" for g in quality["geometry"]):
                    missing.append("catalytic_geometry:"+chain)
            except (ValueError,OSError):missing.append("construct_coordinates:"+chain)
    coordinate_check=None
    if all(k in files for k in ("prmtop","inpcrd","protein_pdb")) and sequences:
        try:
            from .coordinates import validate_amber_coordinates
            coordinate_check=validate_amber_coordinates(files["prmtop"],files["inpcrd"],files["protein_pdb"],
                                sequences,data.get("topology_construct_order"))
        except (ImportError,ValueError,OSError,KeyError) as exc:
            missing.append("actual_simulation_coordinates:"+str(exc))
    return {"status":"passed" if not missing else "blocked","missing":sorted(set(missing)),
            "coordinate_validation":coordinate_check,
            "manifest_sha256":file_hash(path),"files":files,"contract":data}

def prepare_system(input,output):
    """Call an administrator-configured, versioned preparation adapter."""
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    request=json.loads(Path(input).read_text(encoding="utf-8"))
    if request.get("enzyme_state","free")!="free":raise ValueError("Acylated-state preparation is not supported")
    unsupported=[k for k in ("unvalidated_metals","covalent_ligands","nonstandard_residues") if request.get(k)]
    if unsupported:raise ValueError("Preparation stopped: "+", ".join(unsupported))
    adapter_path=os.environ.get("AQUA_PREPARATION_ADAPTER")
    if not adapter_path:raise ValueError("Configure a local preparation adapter with validated molecule parameters")
    adapter=json.loads(Path(adapter_path).read_text(encoding="utf-8"))
    script=checked_file(adapter["script"],Path(adapter_path).parent)
    # The adapter path is deployment configuration, never executable code from a task archive.
    command=[adapter.get("python",sys.executable),str(script),"--request",str(Path(input).resolve()),"--output",str(output.resolve())]
    with (output/"preparation.log").open("ab") as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    result=inspect_prepared(output/"prepared.json")
    atomic_json(output/"preparation_check.json",result)
    if result["status"]!="passed":raise ValueError("Preparation did not pass: "+", ".join(result["missing"]))
    return {"output":str(output),"prepared_manifest":str(output/"prepared.json")}

def runtime_fingerprint(platform_name):
    import openmm
    gpu=None
    if platform_name=="CUDA":
        result=subprocess.run(["nvidia-smi","--query-gpu=uuid,name,driver_version","--format=csv,noheader"],capture_output=True,check=True)
        gpu=result.stdout.decode("utf-8").strip()
    return {"openmm":openmm.__version__,"python":platform.python_version(),"system":platform.platform(),
            "platform":platform_name,"gpu":gpu}

def run_prepared(prepared_manifest,output,production_ns=10.,segment_ns=10.,seed=20261124,
                 platform_name="CUDA",approved_budget_ns=10.,replica_id="replica_1"):
    import openmm
    from openmm import app,unit
    import numpy as np
    # OpenMM 8.2.0 exposes __version__ as "8.2"; validate the installed distribution.
    if distribution_version("OpenMM")!=OPENMM_VERSION:raise ValueError("This protocol requires locked OpenMM "+OPENMM_VERSION)
    if not 0<production_ns<=approved_budget_ns or not 0<segment_ns<=10:
        raise ValueError("Production duration exceeds the explicit budget or segment limit")
    review=inspect_prepared(prepared_manifest)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    atomic_json(output/"preflight.json",review)
    if review["status"]!="passed":raise ValueError("Prepared system is blocked: "+", ".join(review["missing"]))
    cfg=review["contract"];files=review["files"];protocol=cfg["protocol"]
    # Equilibrated Amber coordinates are required; production never silently performs unreviewed repair.
    prmtop=app.AmberPrmtopFile(files["prmtop"]);inpcrd=app.AmberInpcrdFile(files["inpcrd"])
    if inpcrd.boxVectors is None:raise ValueError("A periodic equilibrated system is required")
    dt_fs=float(protocol["timestep_fs"])
    if not 0<dt_fs<=2:raise ValueError("This non-HMR protocol supports timesteps up to 2 fs")
    total=int(round(production_ns*1e6/dt_fs));segment=int(round(segment_ns*1e6/dt_fs))
    interval=int(protocol.get("trajectory_interval_steps",5000))
    if interval<=0 or interval>total:raise ValueError("Trajectory interval exceeds this run")
    system=prmtop.createSystem(nonbondedMethod=app.PME,nonbondedCutoff=1.0*unit.nanometer,
                              constraints=app.HBonds,rigidWater=True)
    barostat=openmm.MonteCarloBarostat(float(protocol["pressure_bar"])*unit.bar,
                                      float(protocol["temperature_K"])*unit.kelvin)
    barostat.setRandomNumberSeed(int(seed));system.addForce(barostat)
    integrator=openmm.LangevinMiddleIntegrator(float(protocol["temperature_K"])*unit.kelvin,
                    float(protocol.get("friction_ps",1))/unit.picosecond,dt_fs*unit.femtosecond)
    integrator.setRandomNumberSeed(int(seed))
    runtime=runtime_fingerprint(platform_name)
    platform_obj=openmm.Platform.getPlatformByName(platform_name)
    properties={"Precision":"mixed","DeviceIndex":"0"} if platform_name=="CUDA" else {}
    simulation=app.Simulation(prmtop.topology,system,integrator,platform_obj,properties)
    simulation.context.setPositions(inpcrd.positions)
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    status_path=output/"progress.json"
    identity={"manifest_sha256":review["manifest_sha256"],"seed":int(seed),"replica_id":replica_id,
              "runtime":runtime,"production_ns":production_ns,"dt_fs":dt_fs}
    progress={"identity":identity,"completed_steps":0,"segments":[]}
    if status_path.exists():
        progress=json.loads(status_path.read_text(encoding="utf-8"))
        if progress["identity"]!=identity:
            raise ValueError("Checkpoint contract or environment changed; use a separately reviewed portable-state continuation")
        checkpoint=checked_file(progress["checkpoint"],output)
        simulation.loadCheckpoint(str(checkpoint))
    else:simulation.context.setVelocitiesToTemperature(float(protocol["temperature_K"])*unit.kelvin,int(seed))
    start=progress["completed_steps"]
    while start<total:
        end=min(total,start+segment);name=f"segment_{start:012d}_{end:012d}"
        dcd=output/(name+".dcd")
        # A segment without an atomically committed progress entry is abandoned after an interruption.
        if dcd.exists():dcd=dcd.with_name(name+"_"+__import__("uuid").uuid4().hex[:8]+".dcd")
        simulation.reporters=[app.DCDReporter(str(dcd),interval,enforcePeriodicBox=True),
                 app.StateDataReporter(str(dcd.with_suffix(".csv")),max(1,min(5000,end-start)),
                        step=True,time=True,potentialEnergy=True,temperature=True,volume=True,speed=True)]
        simulation.step(end-start)
        state=simulation.context.getState(getPositions=True,getEnergy=True)
        if not np.isfinite(state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)).all():
            raise ValueError("Trajectory quality failed: non-finite coordinates")
        if not np.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)):
            raise ValueError("Trajectory quality failed: non-finite energy")
        checkpoint=output/(name+".chk");simulation.saveCheckpoint(str(checkpoint))
        simulation.saveState(str(output/(name+".xml")))
        simulation.reporters=[]
        progress["segments"].append({"path":str(dcd),"sha256":file_hash(dcd),"start_step":start,"end_step":end})
        progress.update(completed_steps=end,checkpoint={"path":checkpoint.name,"sha256":file_hash(checkpoint)})
        atomic_json(status_path,progress);start=end
    summary={"status":"production_completed_awaiting_analysis","replica_id":replica_id,
             "production_ns":production_ns,"seed":seed,"segments":progress["segments"],
             "runtime":runtime,"mechanistic_conclusion":"not_assessed"}
    atomic_json(output/"trajectory_report.json",summary)
    return {"output":str(output),**summary}

def execute_mechanism(output,action="check",**inputs):
    output=Path(output)
    if action=='assess':
        from .mechanism_assessment import assess_batch
        return assess_batch(output=output,**inputs)
    if action=='build-system':
        from .protein_system import prepare_protein_system
        return prepare_protein_system(output=output,**inputs)
    if action=='equilibrate':
        from .equilibration import equilibrate_system
        return equilibrate_system(output=output,**inputs)
    if action=='parameters':
        from .ligand_parameters import parameterize
        return parameterize(output=output,**inputs)
    if action=='solution':
        from .solution_validation import validate_solution
        return validate_solution(output=output,**inputs)
    if action=='solution-analyze':
        from .solution_analysis import analyze_solution
        return analyze_solution(output=output,**inputs)
    if action=='relax-assembly':
        from .assembly_relaxation import relax_assembly
        return relax_assembly(output=output,**inputs)
    if action=="batch":
        from .batch import run_batch
        return run_batch(output=output,**inputs)
    if action=="check":
        result=inspect_prepared(inputs["prepared_manifest"]);output.mkdir(parents=True,exist_ok=True)
        atomic_json(output/"preflight.json",result)
        if result["status"]!="passed":raise ValueError("Preparation blocked: "+", ".join(result["missing"]))
        return {"output":str(output),"status":"prepared"}
    if action=="prepare":return prepare_system(output=output,**inputs)
    if action=="run":return run_prepared(output=output,**inputs)
    if action=="analyze":
        result=analyze_trajectory(output=output,**inputs)
        return {"output":str(output),"status":"analyzed","observables":result["observables"]}
    raise ValueError("Unsupported mechanism action")
