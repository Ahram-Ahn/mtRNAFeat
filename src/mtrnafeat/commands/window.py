"""`mtrnafeat window` — sliding-window scan at a single max_bp_span.

Args (after `--`):
    --window INT       override cfg.window_nt (default: 120)
    --step INT         override cfg.step_nt (default: 10)
    --span INT         override cfg.max_bp_span (default: 300)
"""
from __future__ import annotations

from dataclasses import replace

import pandas as pd

from mtrnafeat.analysis import window
from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import list_genes
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import window_plot
from mtrnafeat.viz.style import plot_path


def _parse(args: list[str] | None) -> dict:
    out: dict = {}
    if not args:
        return out
    it = iter(args)
    for tok in it:
        if tok == "--window":
            out["window_nt"] = int(next(it))
        elif tok == "--step":
            out["step_nt"] = int(next(it))
        elif tok == "--span":
            out["max_bp_span"] = int(next(it))
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    overrides = _parse(args)
    if overrides:
        cfg = replace(cfg, **overrides)
    out = cfg.outdir / "window"
    out.mkdir(parents=True, exist_ok=True)

    all_frames: list[pd.DataFrame] = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        species_genes = set(list_genes(path))
        for gene in cfg.target_genes:
            target = canonical_gene(gene)
            if target not in species_genes:
                continue
            df = window.scan_for_gene(species, path, gene, cfg)
            all_frames.append(df)
            try:
                annot = annotation_for(species, gene)
            except KeyError:
                continue
            window_plot.plot_full_transcript(
                df, annot["l_utr5"], annot["l_cds"], annot["l_tr"],
                species, target, cfg.max_bp_span,
                plot_path(out, f"window_{species}_{target}", cfg.plot_format),
                dpi=cfg.dpi,
            )

    if all_frames:
        merged = pd.concat(all_frames, ignore_index=True)
        canonical_csv(merged, out / "window_scan_metrics.csv")
        canonical_csv(window.summarize_windows(merged), out / "window_scan_summary.csv")
    return 0
