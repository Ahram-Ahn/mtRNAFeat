"""`mtrnafeat significance` — per-gene structural-significance analysis.

Two layers, both restricted to ``cfg.target_genes``:

1. **Sequence-level z-score** (always emitted, fast). Folds each whole
   transcript and a dinucleotide-shuffled null pool, reports
   ``Z_MFE = (observed − null_mean) / null_std``. Output:
   ``z_per_gene.csv``.

2. **Per-gene cotranscriptional / sliding scan** (with ``--scan``).
   For each gene, walks the transcript in the chosen mode and tracks
   Vienna's MFE, ensemble diversity, and paired fraction. Z-scores the
   first-difference signals so structural rearrangement events stand out.
   Outputs one PNG per (species, gene) plus a long-form CSV.

Args (after ``--``):
    --scan                       enable the per-gene cotranscriptional scan
    --mode prefix|sliding        scan mode (default sliding; prefix grows
                                 the prefix length 0..N to mimic nascent RNA)
    --window N                   sliding-mode window size (default 120)
    --step N                     prefix-mode growth step or sliding stride
                                 (default 30)
    --z-threshold T              peak threshold for the smoothed Δ signals
                                 (default 2.0)
    --per-window                 deprecated alias for ``--scan`` (preserved
                                 so ``mtrnafeat run-all`` keeps working)
"""
from __future__ import annotations

from mtrnafeat.analysis import cotrans, significance
from mtrnafeat.config import Config
from mtrnafeat.constants import file_safe_gene
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import cotrans_plot
from mtrnafeat.viz.style import plot_path


def _parse(args: list[str] | None) -> dict:
    out = {
        "scan": False,
        "mode": "sliding",
        "window": 120,
        "step": 30,
        "z_threshold": 2.0,
    }
    if not args:
        return out
    it = iter(args)
    for tok in it:
        if tok in ("--scan", "--per-window"):
            out["scan"] = True
        elif tok == "--mode":
            out["mode"] = next(it)
        elif tok == "--window":
            out["window"] = int(next(it))
        elif tok == "--step":
            out["step"] = int(next(it))
        elif tok == "--z-threshold":
            out["z_threshold"] = float(next(it))
    return out


def run(cfg: Config, args: list[str] | None = None) -> int:
    parsed = _parse(args)
    out = cfg.outdir / "significance"
    out.mkdir(parents=True, exist_ok=True)

    canonical_csv(significance.per_gene_significance(cfg), out / "z_per_gene.csv")

    if not parsed["scan"]:
        return 0

    params = cotrans.ScanParams(
        mode=parsed["mode"],
        window=int(parsed["window"]),
        step=int(parsed["step"]),
        z_threshold=float(parsed["z_threshold"]),
    )
    df = cotrans.per_gene_cotrans_scan(cfg, params)
    if df.empty:
        return 0
    canonical_csv(df, out / "cotrans_per_window.csv")

    # One PNG per (species, gene)
    for (species, gene), gene_df in df.groupby(["Species", "Gene"]):
        cotrans_plot.plot_one_gene(
            gene_df.reset_index(drop=True),
            out_path=plot_path(
                out, f"cotrans_{species}_{file_safe_gene(gene)}", cfg.plot_format
            ),
            z_threshold=params.z_threshold,
            dpi=cfg.dpi,
        )
    return 0
