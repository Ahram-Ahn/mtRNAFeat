"""`mtrnafeat significance` — per-gene structural-significance analysis.

Two layers, both restricted to ``cfg.target_genes``. Only **layer 1** is
a statistical test. Layer 2 is an exploratory within-transcript scan and
its Z-scores are not p-values — see the note below.

1. **Sequence-level z-score / null-model p-value** (always emitted, fast).
   Folds each whole transcript and a dinucleotide-shuffled null pool of
   size ``cfg.n_shuffles`` (Workman & Krogh 1999). Reports
   ``Z_MFE = (observed − null_mean) / null_std`` and the empirical
   ``P_Empirical = mean(null_MFE ≤ observed_MFE)``. Output:
   ``z_per_gene.csv``. **This is the only null-model-backed output of
   this command.**

2. **Per-gene local structural-change scan** (with ``--scan``).
   For each gene, walks the transcript in the chosen mode and tracks
   Vienna's MFE, ensemble diversity, and paired fraction. The
   ``Z_Delta_*`` columns standardize each gene's smoothed first-difference
   signals against *that same gene's* distribution to surface candidate
   transition windows. They are within-gene outlier scores, **not**
   statistical p-values — windows are not independent samples and there
   is no null model. The CSV tags every row with
   ``Z_score_type="within_gene_window_standardized_delta"`` and
   ``Is_statistical_pvalue=False``. For a real null-model test, use
   layer 1's ``P_Empirical``. Outputs one PNG per (species, gene) plus
   ``cotrans_per_window.csv``.

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
        else:
            raise SystemExit(
                f"significance: unknown flag {tok!r}. See module docstring "
                "for the supported flag list."
            )
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
