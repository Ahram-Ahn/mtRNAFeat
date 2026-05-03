"""mtrnafeat CLI: single entry point, subcommand dispatch."""
from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

from mtrnafeat import __version__
from mtrnafeat.config import Config, load_config
from mtrnafeat.core.manifest import write_run_manifest

SUBCOMMANDS = {
    "doctor": "mtrnafeat.commands.doctor",
    "validate-inputs": "mtrnafeat.commands.validate_inputs",
    "stats": "mtrnafeat.commands.stats",
    "landscape": "mtrnafeat.commands.landscape",
    "kinetic": "mtrnafeat.commands.kinetic",
    "window": "mtrnafeat.commands.window",
    "local-probability": "mtrnafeat.commands.local_probability",
    "features": "mtrnafeat.commands.features",
    "significance": "mtrnafeat.commands.significance",
    "structure-deviation": "mtrnafeat.commands.structure_deviation",
    "compare": "mtrnafeat.commands.compare",
    "tis": "mtrnafeat.commands.tis",
    "substitution": "mtrnafeat.commands.substitution",
    "cofold": "mtrnafeat.commands.cofold",
    "gene-panel": "mtrnafeat.commands.gene_panel",
    "plot": "mtrnafeat.commands.plot",
    "run-all": "mtrnafeat.commands.pipeline",
}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mtrnafeat",
        description=f"mtrnafeat {__version__} — mitochondrial mRNA structural-feature analysis",
    )
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)
    for name in SUBCOMMANDS:
        sp = sub.add_parser(name, help=f"Run the {name} subcommand")
        sp.add_argument("--config", type=Path, default=None, help="Path to YAML config")
        sp.add_argument("--outdir", type=Path, default=None, help="Override outdir")
        sp.add_argument("--seed", type=int, default=None, help="Override seed")
        sp.add_argument("rest", nargs=argparse.REMAINDER, help="Subcommand-specific args")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    overrides: dict = {}
    if args.outdir is not None:
        overrides["outdir"] = str(args.outdir)
    if args.seed is not None:
        overrides["seed"] = args.seed
    cfg: Config = load_config(args.config, overrides if overrides else None)
    Path(cfg.outdir).mkdir(parents=True, exist_ok=True)
    if args.config is not None:
        os.environ["MTRNAFEAT_CONFIG_PATH"] = str(args.config)

    module = importlib.import_module(SUBCOMMANDS[args.command])
    rest = args.rest[1:] if args.rest and args.rest[0] == "--" else args.rest

    # Reproducibility manifest: package + tool versions, git commit,
    # command and argv. Written before dispatch so a crashing stage
    # still leaves an audit trail of what was attempted.
    write_run_manifest(
        Path(cfg.outdir),
        command=args.command,
        argv=list(argv) if argv is not None else list(sys.argv[1:]),
        config_path=args.config,
        seed=cfg.seed,
    )

    return int(module.run(cfg, rest) or 0)


if __name__ == "__main__":
    sys.exit(main())
