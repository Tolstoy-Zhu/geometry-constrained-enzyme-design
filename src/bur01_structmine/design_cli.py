"""Command-line entry points for AquaEsterNet local sequence design."""
from __future__ import annotations

import json
from pathlib import Path

from .config import project_root


def add_design_commands(subparsers):
    subparsers.add_parser("design-setup", help="Download and verify the pinned local inverse-folding assets")
    design = subparsers.add_parser("design", help="Design limited substitutions around a Bur01-like catalytic core")
    design.add_argument("--parent", type=Path, default=project_root() / "data/query/Bur01m.pdb")
    design.add_argument("--chain", default="A")
    design.add_argument("--ligand", dest="ligands", action="append", default=[], help="Chain:resname:number, repeat for multiple ligands")
    design.add_argument("--fixed", action="append", default=[], help="Extra protected parent residue such as A240")
    design.add_argument("--msa", type=Path, help="Aligned FASTA/A3M starting with the parent")
    design.add_argument("--sites", dest="site_limit", type=int, default=16)
    design.add_argument("--radius", type=float, default=10.0, help="Design shell radius in angstroms")
    design.add_argument("--max-edits", type=int, default=2)
    design.add_argument("--count", type=int, default=12)
    design.add_argument("--seed", type=int, default=17)
    design.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    design.add_argument("--priority", choices=["compatibility", "ligand_response"], default="compatibility")
    design.add_argument("--output", type=Path, help="A new output directory (default: results/designs/<run_id>)")
    review = subparsers.add_parser("design-review", help="Attribute coordinates to an exact design and compare to its parent")
    review.add_argument("--run", type=Path, required=True)
    review.add_argument("--candidate", required=True)
    review.add_argument("--structure", type=Path, required=True)
    review.add_argument("--chain", default="A")
    review.add_argument("--origin", choices=["unspecified", "experimental", "predicted", "template_model"], default="unspecified")


def dispatch_design(args):
    from .design import review_design, run_design
    from .design_models import install_models
    from .run_context import make_run_id
    try:
        if args.command == "design-setup":
            print(install_models())
        elif args.command == "design-review":
            print(json.dumps(review_design(args.run, args.candidate, args.structure, args.chain, args.origin), indent=2))
        else:
            options = vars(args).copy()
            options.pop("command")
            options["output"] = options["output"] or project_root() / "results/designs" / make_run_id("design")
            print(run_design(**options))
        return 0
    except (ValueError, OSError) as exc:
        import sys
        print(str(exc), file=sys.stderr)
        return 2
