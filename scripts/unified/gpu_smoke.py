"""Numerical installation/restart check; this is not an enzyme simulation."""
import argparse,json,platform,time
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument("--output",type=Path,required=True)
args=p.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
result={"scope":"numerical_gpu_and_checkpoint_smoke_test","python":platform.python_version()}
try:
    import torch
    x=torch.ones((512,512),device="cuda")
    value=float((x@x).mean())
    assert value==512
    result["torch"]={"status":"passed","version":torch.__version__,"cuda":torch.version.cuda,
                     "gpu":torch.cuda.get_device_name(0),"capability":torch.cuda.get_device_capability(0)}
except Exception as exc:result["torch"]={"status":"failed","error":str(exc)}
try:
    import numpy as np
    import openmm
    from openmm import unit
    system=openmm.System();system.addParticle(12);system.addParticle(12)
    force=openmm.HarmonicBondForce();force.addBond(0,1,0.15,1000);system.addForce(force)
    integrator=openmm.LangevinMiddleIntegrator(300,1,0.002);integrator.setRandomNumberSeed(2026)
    context=openmm.Context(system,integrator,openmm.Platform.getPlatformByName("CUDA"),{"Precision":"mixed","DeviceIndex":"0"})
    context.setPositions([[0,0,0],[.16,0,0]]);context.setVelocitiesToTemperature(300,2026)
    integrator.step(20);checkpoint=context.createCheckpoint()
    integrator.step(30);first=context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    context.loadCheckpoint(checkpoint);integrator.step(30)
    second=context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    error=float(np.max(np.abs(first-second)))
    assert error<1e-5
    result["openmm"]={"status":"passed","version":openmm.__version__,"platform":context.getPlatform().getName(),
                      "restart_max_position_difference_nm":error,"steps":80}
except Exception as exc:result["openmm"]={"status":"failed","error":str(exc)}
result["status"]="passed" if all(result[k]["status"]=="passed" for k in ["torch","openmm"]) else "failed"
args.output.write_text(json.dumps(result,indent=2),encoding="utf-8")
print(json.dumps(result,indent=2))
raise SystemExit(0 if result["status"]=="passed" else 1)
