"""Local property model setup, batch prediction and interactive browser access."""

import subprocess
import sys
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
    profile = parsers.add_parser("properties", help="Profile AquaDesign candidates and their parents")
    profile.add_argument(
        "--run",
        dest="runs",
        type=Path,
        action="append",
        required=True,
        help="Repeat for multiple design runs",
    )
    profile.add_argument("--substrate", action="append", help="Name=SMILES, repeat for each substrate")
    profile.add_argument("--output", type=Path, help="New local report directory")
    profile.add_argument(
        "--device", choices=["cpu", "cuda"], default="cpu", help="TmProt device; EnzRank uses CPU"
    )
    web = parsers.add_parser(
        "properties-web", help="Open the local property prediction and comparison workbench"
    )
    web.add_argument("--port", type=int, default=8765)


def dispatch_property(args):
    from .property_assets import install_assets, install_runtime
    from .property_screen import DEFAULT_SUBSTRATES, parse_substrates, profile_runs

    try:
        if args.command == "properties-setup":
            if args.install_runtime:
                print(install_runtime())
            print(install_assets())
        elif args.command == "properties-web":
            if not 1 <= args.port <= 65535:
                raise ValueError("Invalid local port.")
            return subprocess.call(
                [
                    sys.executable,
                    "-m",
                    "streamlit",
                    "run",
                    str(project_root() / "ui/properties_app.py"),
                    "--server.address",
                    "127.0.0.1",
                    "--server.port",
                    str(args.port),
                    "--server.headless",
                    "false",
                    "--browser.gatherUsageStats",
                    "false",
                ]
            )
        else:
            pairs = parse_substrates(args.substrate) if args.substrate else DEFAULT_SUBSTRATES
            print(profile_runs(args.runs, pairs, args.output, device=args.device, progress=print))
        return 0
    except (ValueError, OSError, KeyError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
