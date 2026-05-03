"""`mtrnafeat window` — whole-transcript fold-and-compare.

Per (species, gene), folds the transcript three ways and emits a per-
position paired-fraction track plus a summary CSV:

* DMS (from the .db file, energy recalculated under Vienna)
* Vienna full (no max_bp_span, pure thermodynamic baseline)
* Engine span (`--engine vienna|rnastructure` × `--span N`)

Args (after `--`):
    --span INT             override cfg.max_bp_span (default: 300)
    --rolling-window INT   override cfg.rolling_window (default: 25)
    --engine NAME          "rnastructure" (default) | "vienna"
"""
from __future__ import annotations

from dataclasses import replace

import pandas as pd

from mtrnafeat.analysis import window
from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene, file_safe_gene
from mtrnafeat.io.db_parser import list_genes
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import window_plot
from mtrnafeat.viz.style import plot_path

_VALID_ENGINES = ("rnastructure", "vienna")


def _parse(args: list[str] | None) -> dict:
    out: dict = {}
    if not args:
        return out
    it = iter(args)
    for tok in it:
        if tok == "--span":
            out["max_bp_span"] = int(next(it))
        elif tok == "--rolling-window":
            out["rolling_window"] = int(next(it))
        elif tok == "--engine":
            out["fold_engine"] = next(it)
        else:
            raise SystemExit(
                f"window: unknown flag {tok!r}. "
                "Supported: --span, --rolling-window, --engine."
            )
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    overrides = _parse(args)
    if overrides:
        cfg = replace(cfg, **overrides)
    if cfg.fold_engine not in _VALID_ENGINES:
        raise ValueError(
            f"Unknown fold_engine: {cfg.fold_engine!r} "
            f"(expected one of {_VALID_ENGINES})"
        )
    out = cfg.outdir / "window"
    out.mkdir(parents=True, exist_ok=True)

    pos_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    failed: list[dict] = []
    coverage: dict[str, dict[str, list[str]]] = {}

    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        species_genes = set(list_genes(path))
        configured = [canonical_gene(g) for g in cfg.target_genes]
        coverage[species] = {
            "found": [g for g in configured if g in species_genes],
            "missing": [g for g in configured if g not in species_genes],
        }
        for gene in cfg.target_genes:
            target = canonical_gene(gene)
            if target not in species_genes:
                print(f"[mtrnafeat] window: skip {species} {target} (not in {fname})", flush=True)
                continue
            window.step_log(species, gene)
            try:
                res = window.fold_transcript(species, path, gene, cfg)
            except Exception as exc:  # noqa: BLE001 — keep batch going
                print(f"[mtrnafeat] window: FAIL {species} {target}: {exc}", flush=True)
                failed.append({"Species": species, "Gene": target, "Error": str(exc)})
                continue
            pos_df = window.per_position_table(res, cfg.rolling_window)
            summary_df = window.summarize_transcript(res)
            pos_frames.append(pos_df)
            summary_frames.append(summary_df)
            window_plot.plot_transcript_pairing(
                res, pos_df,
                out_path=plot_path(out, f"window_{species}_{file_safe_gene(gene)}", cfg.plot_format),
                rolling_window=cfg.rolling_window,
                dpi=cfg.dpi,
            )

    # Summarize per-species gene coverage so silent skips become visible.
    for species, cov in coverage.items():
        n_found = len(cov["found"])
        n_missing = len(cov["missing"])
        msg = f"[mtrnafeat] window: {species} coverage — {n_found} found"
        if n_missing:
            msg += f", {n_missing} missing ({', '.join(cov['missing'])})"
        print(msg, flush=True)

    if not pos_frames and not failed:
        # Every (species, gene) pair was silently skipped — no .db gene
        # matched the target list. Refuse to write empty output so users
        # don't mistake an empty run for a successful one.
        raise SystemExit(
            "window: no transcripts processed — every configured gene was "
            "missing from the .db files. Check `target_genes:` and `db_files:` "
            "in your config."
        )

    if pos_frames:
        canonical_csv(pd.concat(pos_frames, ignore_index=True),
                      out / "window_per_position.csv")
    if summary_frames:
        canonical_csv(pd.concat(summary_frames, ignore_index=True),
                      out / "window_summary.csv")
    if failed:
        canonical_csv(pd.DataFrame(failed), out / "window_failed.csv")
    return 0
