"""`mtrnafeat run-all` — orchestrates the full analysis pipeline.

Every stage is independent; there are no cross-stage dependencies after
the round-3 cleanup that retired nascent / signals / concordance.

`--parallel` fires every command concurrently as a subprocess.

DrTransformer kinetic folding is opt-in only — never auto-runs in
`run-all`. Invoke `mtrnafeat kinetic` explicitly if you have the binary.

Usage:
    mtrnafeat run-all --config configs/all.yaml --outdir runs/all
    mtrnafeat run-all --parallel --config configs/all.yaml --outdir runs/all
    mtrnafeat run-all --parallel --skip significance,cofold --config ...
"""
from __future__ import annotations

import importlib
import multiprocessing as mp
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from mtrnafeat.config import Config


INDEPENDENT = (
    "stats", "landscape", "features",
    "window", "significance", "tis", "compare",
    "substitution", "cofold", "gene_panel",
)

# Per-stage extra args appended after `--` when invoking the subcommand.
# significance: --per-window unlocks the refined Panel-B-style overlays;
# without it, the stage emits only z_per_gene.csv.
STAGE_EXTRAS: dict[str, list[str]] = {
    "significance": ["--per-window"],
}


def _parse(args: list[str] | None) -> dict:
    parsed = {"parallel": False, "skip": set()}
    if not args:
        return parsed
    it = iter(args)
    for tok in it:
        if tok == "--parallel":
            parsed["parallel"] = True
        elif tok == "--skip":
            parsed["skip"].update(next(it).split(","))
    return parsed


def _run_subcommand_in_process(name: str, config_path: str | None, outdir: str, seed: int | None) -> tuple[str, int, int]:
    """Subprocess entry — invoked by the parallel pool."""
    cmd = [sys.executable, "-m", "mtrnafeat.cli", name.replace("_", "-"),
           "--outdir", outdir]
    if config_path:
        cmd += ["--config", config_path]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    extras = STAGE_EXTRAS.get(name)
    if extras:
        cmd += ["--", *extras]
    started = time.time()
    proc = subprocess.run(cmd, capture_output=False)
    return name, int(time.time() - started), int(proc.returncode)


def _sequential(cfg: Config, skip: set[str]) -> None:
    print(f"[mtrnafeat] running pipeline (sequential) → {cfg.outdir}")
    failures: list[tuple[str, str]] = []
    for name in INDEPENDENT:
        if name in skip:
            print(f"  - {name}: SKIPPED")
            continue
        print(f"  → {name}")
        mod = importlib.import_module(f"mtrnafeat.commands.{name}")
        try:
            rc = int(mod.run(cfg, list(STAGE_EXTRAS.get(name, []))) or 0)
        except Exception as exc:
            failures.append((name, f"exception: {exc!r}"))
            print(f"  ✗ {name} RAISED: {exc!r}")
            continue
        if rc != 0:
            failures.append((name, f"exit {rc}"))
            print(f"  ✗ {name} FAILED (exit {rc})")
    if failures:
        print(f"[mtrnafeat] pipeline complete WITH FAILURES: " + ", ".join(f"{n} ({why})" for n, why in failures))
    else:
        print("[mtrnafeat] pipeline complete.")


def _parallel(cfg: Config, skip: set[str], config_path: str | None) -> None:
    print(f"[mtrnafeat] running pipeline (parallel) → {cfg.outdir}")
    n_indep = max(2, min(len(INDEPENDENT) + 1, mp.cpu_count()))
    outdir = str(cfg.outdir)
    seed = cfg.seed

    independent_to_run = [n for n in INDEPENDENT if n not in skip]
    finished: dict[str, int] = {}
    failures: list[tuple[str, int]] = []
    started = time.time()

    with ProcessPoolExecutor(max_workers=n_indep) as ex:
        futures = {
            ex.submit(_run_subcommand_in_process, name, config_path, outdir, seed): name
            for name in independent_to_run
        }
        for fut in as_completed(futures):
            name, elapsed, rc = fut.result()
            finished[name] = elapsed
            wall = int(time.time() - started)
            if rc == 0:
                print(f"[mtrnafeat]   ✓ {name} done in {elapsed}s (wall {wall}s)")
            else:
                failures.append((name, rc))
                print(f"[mtrnafeat]   ✗ {name} FAILED (exit {rc}) after {elapsed}s (wall {wall}s)")

    total = int(time.time() - started)
    print(f"[mtrnafeat] pipeline complete — wall {total}s")
    print(f"[mtrnafeat] per-stage time: " + ", ".join(f"{k}={v}s" for k, v in finished.items()))
    if failures:
        print(f"[mtrnafeat] FAILURES: " + ", ".join(f"{n} (exit {rc})" for n, rc in failures))


def run(cfg: Config, args: list[str] | None = None) -> int:
    parsed = _parse(args)
    skip = set(parsed["skip"])
    if parsed["parallel"]:
        config_path = os.environ.get("MTRNAFEAT_CONFIG_PATH")
        _parallel(cfg, skip, config_path)
    else:
        _sequential(cfg, skip)
    return 0
