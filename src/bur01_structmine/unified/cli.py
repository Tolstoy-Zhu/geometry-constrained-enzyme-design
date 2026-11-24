"""Project, training and mechanism commands sharing the workbench service layer."""
from __future__ import annotations
import json
from pathlib import Path

COMMANDS={"projects","tasks","dataset","train","mechanism","records","server","families","kinetics"}

def add_commands(subparsers):
    families=subparsers.add_parser("families",help="List configured enzyme families")
    projects=subparsers.add_parser("projects",help="Manage enzyme-engineering projects")
    projects.add_argument("action",choices=["list","create","index"])
    projects.add_argument("--name");projects.add_argument("--family",default="bur01")
    projects.add_argument("--executor",default="local");projects.add_argument("--root",default=".")
    projects.add_argument("--project",default="default")
    tasks=subparsers.add_parser("tasks",help="Inspect durable local and SSH tasks")
    tasks.add_argument("action",choices=["list","show","sync","cancel","retry"])
    tasks.add_argument("--task");tasks.add_argument("--project",default="default")
    server=subparsers.add_parser("server",help="Local-only SSH server profiles")
    server.add_argument("action",choices=["list","save","check"])
    server.add_argument("--profile",type=Path);server.add_argument("--name")
    records=subparsers.add_parser("records",help="Import validated experiments and export kinetic data")
    records.add_argument("action",choices=["import","export-assays"])
    records.add_argument("--input",type=Path);records.add_argument("--output",type=Path)
    records.add_argument("--project",default="default")
    for name in ["dataset","train","mechanism","kinetics"]:
        cmd=subparsers.add_parser(name,help="Submit "+name+" operations")
        cmd.add_argument("action",choices=({"dataset":["fetch-catapro","fetch-sabio","map-sabio","build","features"],
                                          "train":["submit","register"],
                                          "mechanism":["plan","parameters","solution","solution-analyze","relax-assembly","check","prepare","run","analyze","batch"],
                                          "kinetics":["predict"]})[name])
        cmd.add_argument("--spec",type=Path,help="JSON input fields; paths are validated and staged")
        cmd.add_argument("--input",type=Path);cmd.add_argument("--output",type=Path)
        cmd.add_argument("--query");cmd.add_argument("--project",default="default")
        cmd.add_argument("--family",default="bur01");cmd.add_argument("--executor",default="local")
        cmd.add_argument("--reaction-config",type=Path)

def dispatch(args):
    from .store import Store,index_history
    from .execution import submit,server_profiles,save_server,ServerProfile,RemoteExecutor
    from .services import import_records,export_assays,reconcile_task,cancel_task,retry_task,register_model
    from .families import profiles,get_profile
    result=None;store=Store()
    if args.command=="families":result=[p.model_dump() for p in profiles()]
    elif args.command=="projects":
        if args.action=="list":result=store.projects()
        elif args.action=="create":
            get_profile(args.family)
            result={"project_id":store.create_project(args.name or "",args.family,args.executor)}
        else:result={"indexed":index_history(args.root,store,args.project)}
    elif args.command=="server":
        if args.action=="list":result=[p.model_dump() for p in server_profiles().values()]
        elif args.action=="save":
            profile=ServerProfile.model_validate_json(args.profile.read_text(encoding="utf-8"))
            save_server(profile);result={"saved":profile.name}
        else:result=RemoteExecutor(server_profiles()[args.name]).check()
    elif args.command=="tasks":
        if args.action=="list":result=store.tasks(args.project)
        elif args.action=="show":result=store.task(args.task)
        elif args.action=="sync":result=reconcile_task(args.task,store)
        elif args.action=="cancel":
            cancel_task(args.task,args.project,store);result={"cancel_requested":args.task}
        else:result={"task_id":retry_task(args.task,store)}
    elif args.command=="records":
        result=import_records(args.input,args.project,store) if args.action=="import" else export_assays(args.output,args.project,store)
    elif args.command=="dataset" and args.action.startswith("fetch-"):
        from .datasets import fetch_catapro,fetch_sabio
        result=fetch_catapro(args.output) if args.action=="fetch-catapro" else fetch_sabio(args.query,args.output)
    elif args.command=="train" and args.action=="register":
        result=register_model(args.input,args.project,store)
    elif args.command=="dataset" and args.action=="map-sabio":
        from .sabio import map_assays
        result=map_assays(args.input,args.spec,args.output)
    else:
        inputs=json.loads(args.spec.read_text(encoding="utf-8")) if args.spec else {}
        if args.command=="mechanism" and args.action=="plan":
            from .mechanism_plan import freeze_batch
            result=freeze_batch(inputs,args.output)
        else:
            kind="features" if args.action=="features" else args.command
            if kind=="mechanism":inputs["action"]=args.action
            if args.input:inputs["input"]=str(args.input)
            result={"task_id":submit(kind,inputs,project_id=args.project,family_id=args.family,executor=args.executor,store=store,
                                     reaction_config=args.reaction_config)}
    print(json.dumps(result,ensure_ascii=False,indent=2,default=str))
    return 0
