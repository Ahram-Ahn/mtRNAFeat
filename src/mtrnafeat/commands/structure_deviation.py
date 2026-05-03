"""`mtrnafeat structure-deviation` — region-discovery on the
RNAplfold-vs-DMS deviation track.

For every (species, gene) target, computes the smoothed signed
deviation ``P_model(i) − P_DMS(i)``, calls intervals where the
deviation passes a threshold, classifies each interval into one of
{model_high_dms_low, model_low_dms_high, concordant_paired,
concordant_open, mixed_deviation, ambiguous}, and emits:

    structure_deviation_per_position.csv
    structure_deviation_regions.csv
    structure_deviation_gene_summary.csv
    structure_deviation_gene_region_matrix.csv
    structure_deviation_{species}_{gene}.{svg|png}    (one per gene)
    structure_deviation_lollipop_{species}.{svg|png}  (one per species)
    structure_deviation_heatmap.{svg|png}             (cross-gene)

This stage is the interpretation / region-discovery layer on top of
``local-probability``. The two are complementary:

    local-probability  — "what does local pairing probability look like?"
    structure-deviation — "where does DMS diverge from the local model,
                          and what biological class is each region?"

Args (after ``--``):
    --threshold F           deviation magnitude threshold for region
                            calling (default cfg.structure_deviation_threshold)
    --min-region-length N   minimum region length in nt
    --merge-gap N           merge gap between same-sign regions
    --high-threshold F      paired-fraction high cutoff for class labels
    --low-threshold F       paired-fraction low cutoff for class labels
    --top-labels N          how many top-|Δ| regions to label per gene plot
    --null {none,dinuc}     null model for per-region empirical p-values
                            (default cfg.structure_deviation_null_model = "none";
                            "dinuc" runs an Altschul-Erikson dinucleotide-shuffle
                            null + Westfall-Young max-statistic correction within
                            gene + Benjamini-Hochberg q-values across genes)
    --n-null N              null draws per gene when --null dinuc
                            (default cfg.structure_deviation_n_null;
                            100 is a reasonable starting point)
    --no-plots              skip plot generation (CSV only)
    --only-species NAME     restrict to one species
    --only-gene NAME        restrict to one gene
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.analysis import deviation, local_probability as lp_analysis
from mtrnafeat.config import Config
from mtrnafeat.constants import file_safe_gene
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import structure_deviation as viz
from mtrnafeat.viz.style import plot_path


def _parse(args: list[str] | None) -> dict:
    out: dict = {}
    if not args:
        return out
    it = iter(args)
    for tok in it:
        if tok == "--threshold":
            out["threshold"] = float(next(it))
        elif tok == "--min-region-length":
            out["min_region_length"] = int(next(it))
        elif tok == "--merge-gap":
            out["merge_gap"] = int(next(it))
        elif tok == "--high-threshold":
            out["high_threshold"] = float(next(it))
        elif tok == "--low-threshold":
            out["low_threshold"] = float(next(it))
        elif tok == "--top-labels":
            out["top_labels"] = int(next(it))
        elif tok == "--null":
            out["null"] = str(next(it))
        elif tok == "--n-null":
            out["n_null"] = int(next(it))
        elif tok == "--no-plots":
            out["no_plots"] = True
        elif tok == "--only-species":
            out["only_species"] = str(next(it))
        elif tok == "--only-gene":
            out["only_gene"] = str(next(it))
        else:
            raise SystemExit(
                f"structure-deviation: unknown flag {tok!r}. "
                "See module docstring or `mtrnafeat structure-deviation --help`."
            )
    return out


def _override_cfg(cfg: Config, parsed: dict) -> Config:
    """Apply CLI overrides without mutating the input config."""
    fields_map = {
        "threshold": "structure_deviation_threshold",
        "min_region_length": "structure_deviation_min_region_length",
        "merge_gap": "structure_deviation_merge_gap",
        "high_threshold": "structure_deviation_high_threshold",
        "low_threshold": "structure_deviation_low_threshold",
        "top_labels": "structure_deviation_top_labels",
        "null": "structure_deviation_null_model",
        "n_null": "structure_deviation_n_null",
    }
    overrides: dict = {}
    for cli_key, cfg_key in fields_map.items():
        if cli_key in parsed:
            overrides[cfg_key] = parsed[cli_key]
    if not overrides:
        return cfg
    from dataclasses import replace
    return replace(cfg, **overrides)


def _filter_targets(cfg: Config, parsed: dict) -> Config:
    """Apply ``--only-species`` / ``--only-gene`` filters."""
    from dataclasses import replace
    db = dict(cfg.db_files)
    if "only_species" in parsed:
        wanted = parsed["only_species"]
        db = {k: v for k, v in db.items() if k == wanted}
        if not db:
            raise SystemExit(f"--only-species {wanted!r} not in cfg.db_files")
    targets = cfg.target_genes
    if "only_gene" in parsed:
        targets = (parsed["only_gene"],)
    return replace(cfg, db_files=db, target_genes=targets)


def run(cfg: Config, args: list[str] | None = None) -> int:
    parsed = _parse(args)
    cfg = _override_cfg(cfg, parsed)
    cfg = _filter_targets(cfg, parsed)
    out = cfg.outdir / "structure_deviation"
    out.mkdir(parents=True, exist_ok=True)

    bundle = deviation.scan_all(cfg)
    if bundle["per_position"].empty:
        return 0

    canonical_csv(bundle["per_position"],
                  out / "structure_deviation_per_position.csv")
    canonical_csv(bundle["regions"],
                  out / "structure_deviation_regions.csv")
    canonical_csv(bundle["gene_summary"],
                  out / "structure_deviation_gene_summary.csv")
    canonical_csv(bundle["gene_region_matrix"],
                  out / "structure_deviation_gene_region_matrix.csv")

    if parsed.get("no_plots"):
        return 0

    regions_df: pd.DataFrame = bundle["regions"]

    # Build the same per-window aggregation that ``local-probability``
    # plots use, so both stages display the deviation track at the same
    # scale. Region calling itself is unchanged — it still runs on the
    # per-position 25-nt rolling track via ``deviation.scan_all``. The
    # per-window frame here is purely a visualization aid.
    win_frames: list[pd.DataFrame] = []
    scan_window = int(cfg.local_probability_scan_window_nt)
    scan_step = int(cfg.local_probability_scan_step_nt)
    for result in bundle["results"]:
        # Build a LocalProbResult-shaped object via the existing
        # ``compute_one_gene`` returned data; the per-window helper
        # only reads ``.species``, ``.gene``, ``.sequence``,
        # ``.p_paired``, ``.dms_paired_binary``, ``.dms_structure``.
        lp_like = lp_analysis.LocalProbResult(
            species=result.species,
            gene=result.gene,
            sequence=result.sequence,
            p_paired=result.p_model_raw,
            window=result.rnaplfold_window,
            max_bp_span=result.rnaplfold_max_bp_span,
            cutoff=result.rnaplfold_cutoff,
            dms_structure=result.dms_structure,
            dms_paired_binary=result.p_dms_raw.astype("int8"),
        )
        try:
            annot = annotation_for(result.species, result.gene)
        except KeyError:
            annot = None
        win_df = lp_analysis.per_window_agreement_table(
            lp_like, annot, win=scan_window, step_nt=scan_step,
        )
        if not win_df.empty:
            win_frames.append(win_df)
    per_window_all = (pd.concat(win_frames, ignore_index=True)
                      if win_frames else pd.DataFrame())

    # Per-gene 4-panel plots
    for result in bundle["results"]:
        gene_regions = (
            regions_df[(regions_df["Species"] == result.species)
                       & (regions_df["Gene"] == result.gene)]
            if not regions_df.empty else pd.DataFrame()
        )
        viz.plot_one_gene(
            result, gene_regions,
            out_path=plot_path(
                out,
                f"structure_deviation_{result.species}_{file_safe_gene(result.gene)}",
                cfg.plot_format,
            ),
            cfg=cfg, dpi=cfg.dpi,
            per_window_df=per_window_all,
            scan_window=scan_window,
        )

    # Per-species lollipop summary
    if not regions_df.empty:
        for species in sorted(regions_df["Species"].unique()):
            viz.plot_lollipop(
                regions_df,
                out_path=plot_path(
                    out,
                    f"structure_deviation_lollipop_{species}",
                    cfg.plot_format,
                ),
                species=species, cfg=cfg, dpi=cfg.dpi,
            )

    # Cross-gene heatmap
    if not bundle["gene_region_matrix"].empty:
        viz.plot_heatmap(
            bundle["gene_region_matrix"],
            out_path=plot_path(
                out, "structure_deviation_heatmap", cfg.plot_format,
            ),
            value="Mean_Deviation",
            cfg=cfg, dpi=cfg.dpi,
        )
    return 0
