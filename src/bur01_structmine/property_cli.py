"""Local property model setup, batch prediction and interactive browser access."""

import json
import subprocess
import sys
import urllib.request
import webbrowser
from pathlib import Path

from .config import project_root


def add_property_commands(parsers):
    setup = parsers.add_parser(
        "properties-setup", help="Fetch and verify local SWI, TmProt and EnzRank model assets"
    )
    setup.add_argument(
        "--install-runtime",
        action="store_true",
        help="Also create and install the isolated Python environment",
    )
    profile = parsers.add_parser("properties", help="Profile protein sequences or AquaDesign batches")
    inputs = profile.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--run",
        dest="runs",
        type=Path,
        action="append",
        help="Repeat for multiple design runs",
    )
    inputs.add_argument("--fasta", type=Path, help="One or more protein sequences in FASTA format")
    profile.add_argument(
        "--first-as-parent",
        action="store_true",
        help="Compare aligned-length FASTA variants with the first record",
    )
    profile.add_argument("--substrate", action="append", help="Name=SMILES, repeat for each substrate")
    profile.add_argument("--output", type=Path, help="New local report directory")
    profile.add_argument(
        "--device", choices=["cpu", "cuda"], default="cpu", help="TmProt device; EnzRank uses CPU"
    )
    for command, port, app, description in (
        ("web", 8501, "app.py", "Open StructureMine with mining, sequence design and property prediction tabs"),
        ("properties-web", 8765, "properties_app.py", "Open the standalone property pages"),
    ):
        web = parsers.add_parser(command, help=description)
        web.add_argument("--port", type=int, default=port)
        web.add_argument("--headless", action="store_true", help="Serve the app without opening a browser window")
        web.set_defaults(web_app=app)


def existing_workbench(url):
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        with opener.open(url + "/app/static/aquaesternet-web.json", timeout=1) as response:
            return json.loads(response.read(1024)) == {"app": "aquaesternet-model-web", "version": 1}
    except (ValueError, OSError):
        return False


def dispatch_property(args):
    from .property_assets import install_assets, install_runtime
    from .property_screen import DEFAULT_SUBSTRATES, parse_substrates, profile_fasta, profile_runs

    try:
        if args.command == "properties-setup":
            if args.install_runtime:
                print(install_runtime())
            print(install_assets())
        elif args.command in {"web", "properties-web"}:
            if not 1 <= args.port <= 65535:
                raise ValueError("Invalid local port.")
            url = f"http://127.0.0.1:{args.port}"
            if existing_workbench(url):
                print(url)
                if not args.headless:
                    webbrowser.open(url)
                return 0
            return subprocess.call(
                [
                    sys.executable,
                    "-m",
                    "streamlit",
                    "run",
                    str(project_root() / "ui" / args.web_app),
                    "--server.address",
                    "127.0.0.1",
                    "--server.port",
                    str(args.port),
                    "--server.headless",
                    str(args.headless).lower(),
                    "--server.enableStaticServing",
                    "true",
                    "--browser.gatherUsageStats",
                    "false",
                ],
                cwd=project_root(),
            )
        else:
            pairs = parse_substrates(args.substrate) if args.substrate else DEFAULT_SUBSTRATES
            if args.fasta:
                print(
                    profile_fasta(
                        args.fasta.read_text(encoding="utf-8-sig"),
                        pairs,
                        args.output,
                        first_as_parent=args.first_as_parent,
                        device=args.device,
                        progress=print,
                    )
                )
            elif args.first_as_parent:
                raise ValueError("--first-as-parent requires --fasta.")
            else:
                print(profile_runs(args.runs, pairs, args.output, device=args.device, progress=print))
        return 0
    except (ValueError, OSError, KeyError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
