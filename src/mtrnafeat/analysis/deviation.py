"""Structure-deviation analysis: where does the DMS-derived native
structure diverge from local thermodynamic folding potential, and
what biological class does each region belong to?

Per-position signal::

    P_model(i) = RNAplfold local pairing probability
    P_DMS(i)   = DMS dot-bracket paired fraction (smoothed)
    deviation(i) = P_model(i) − P_DMS(i)

Region calling on the smoothed deviation track yields intervals; each
region is classified into one of:

* ``model_high_dms_low``      thermodynamically foldable but DMS-open
                              (initiation, ribosome opening, RBP/helicase)
* ``model_low_dms_high``      DMS-paired/protected beyond local model
                              (RBP protection, nonlocal/RNP-stabilized pair)
* ``concordant_paired``       locally encoded paired structural element
* ``concordant_open``         intrinsically open / unstructured segment
* ``mixed_deviation``         signed deviation present, intermediate signal
* ``ambiguous``               weak signal — do not overinterpret

This module is the analysis layer; plotting lives in
``mtrnafeat.viz.structure_deviation`` and the CLI in
``mtrnafeat.commands.structure_deviation``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from mtrnafeat.analysis.local_probability import (
    _centered_rolling_mean,
    paired_binary_from_dotbracket,
    scan_one_gene,
)
from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core.shuffle import dinuc_shuffle
from mtrnafeat.core.stats import bh_fdr
from mtrnafeat.core.structure import extract_pairs
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import DbRecord, parse_db
from mtrnafeat.progress import progress, step
from mtrnafeat.rng import make_rng

# Pair-span cutoffs used to characterize the long-range fraction of
# DMS-derived pairs. These follow standard short/medium/long-range
# cutoffs in the RNA-structure literature and define the corresponding
# CSV column names (``Fraction_DMS_Pairs_Span_gt_50`` etc.); kept as
# module constants rather than config so the published table schema is
# stable across runs.
_LONG_PAIR_SPAN_CUTOFFS: tuple[int, ...] = (50, 100, 300)

# ──────────────────────── Result container ────────────────────────


@dataclass
class DeviationResult:
    species: str
    gene: str
    sequence: str
    dms_structure: str
    p_model_raw: np.ndarray
    p_model_smooth: np.ndarray
    p_dms_raw: np.ndarray  # 0/1 binary (paired in DMS dot-bracket)
    p_dms_smooth: np.ndarray
    deviation_raw: np.ndarray
    deviation_smooth: np.ndarray
    rolling_window: int
    rnaplfold_window: int
    rnaplfold_max_bp_span: int
    rnaplfold_cutoff: float


# ──────────────────────── Per-gene compute ────────────────────────


def compute_one_gene(species: str, gene: str, sequence: str,
                     dms_structure: str, *, cfg: Config) -> DeviationResult:
    """Run RNAplfold + DMS overlay and build the smoothed deviation track."""
    if len(sequence) != len(dms_structure):
        raise ValueError(
            f"{species} {gene}: sequence length {len(sequence)} "
            f"!= DMS structure length {len(dms_structure)}"
        )
    res = scan_one_gene(
        species, gene, sequence,
        window=cfg.rnaplfold_window,
        max_bp_span=cfg.rnaplfold_max_bp_span,
        cutoff=cfg.rnaplfold_cutoff,
        dms_structure=dms_structure,
    )
    p_model_raw = res.p_paired
    p_dms_raw = paired_binary_from_dotbracket(dms_structure).astype(float)
    smooth_w = int(cfg.rolling_window)
    p_model_smooth = _centered_rolling_mean(p_model_raw, smooth_w)
    p_dms_smooth = _centered_rolling_mean(p_dms_raw, smooth_w)
    deviation_raw = p_model_raw - p_dms_raw
    deviation_smooth = p_model_smooth - p_dms_smooth
    return DeviationResult(
        species=species,
        gene=canonical_gene(gene),
        sequence=sequence,
        dms_structure=dms_structure,
        p_model_raw=p_model_raw,
        p_model_smooth=p_model_smooth,
        p_dms_raw=p_dms_raw,
        p_dms_smooth=p_dms_smooth,
        deviation_raw=deviation_raw,
        deviation_smooth=deviation_smooth,
        rolling_window=smooth_w,
        rnaplfold_window=int(cfg.rnaplfold_window),
        rnaplfold_max_bp_span=int(cfg.rnaplfold_max_bp_span),
        rnaplfold_cutoff=float(cfg.rnaplfold_cutoff),
    )


# ──────────────────────── Per-position table ────────────────────────


def _local_composition(sequence: str, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Centered rolling GC and AU fraction at every position."""
    n = len(sequence)
    gc = np.fromiter(((1.0 if c in "GC" else 0.0) for c in sequence),
                     dtype=float, count=n)
    au = np.fromiter(((1.0 if c in "AU" else 0.0) for c in sequence),
                     dtype=float, count=n)
    return _centered_rolling_mean(gc, width), _centered_rolling_mean(au, width)


def per_position_table(result: DeviationResult,
                       annot: dict | None = None) -> pd.DataFrame:
    n = len(result.sequence)
    pos_1 = np.arange(1, n + 1)
    pos_0 = np.arange(0, n)

    if annot is not None:
        cds_start_1 = int(annot["l_utr5"]) + 1
        cds_end_1 = int(annot["l_utr5"]) + int(annot["l_cds"])
    else:
        cds_start_1 = None
        cds_end_1 = None

    transcript_region: list[str] = []
    distance_to_start: list[float] = []
    distance_to_stop: list[float] = []
    cds_phase: list[float] = []
    for p in pos_1:
        if cds_start_1 is None:
            transcript_region.append("unknown")
            distance_to_start.append(np.nan)
            distance_to_stop.append(np.nan)
            cds_phase.append(np.nan)
            continue
        if p < cds_start_1:
            transcript_region.append("5UTR")
        elif p <= cds_end_1:
            transcript_region.append("CDS")
        else:
            transcript_region.append("3UTR")
        distance_to_start.append(int(p) - cds_start_1)
        distance_to_stop.append(int(p) - cds_end_1)
        if cds_start_1 <= p <= cds_end_1:
            cds_phase.append((p - cds_start_1) % 3)
        else:
            cds_phase.append(np.nan)

    gc_smooth, au_smooth = _local_composition(result.sequence, result.rolling_window)

    return pd.DataFrame({
        "Species": result.species,
        "Gene": result.gene,
        "Position_1based": pos_1,
        "Position_0based": pos_0,
        "Nucleotide": list(result.sequence),
        "Transcript_Region": transcript_region,
        "Distance_To_Start": distance_to_start,
        "Distance_To_Stop": distance_to_stop,
        "CDS_Phase": cds_phase,
        "DMS_Paired_Raw": result.p_dms_raw.astype(np.int8),
        "DMS_Paired_Smooth": result.p_dms_smooth,
        "Model_Paired_Raw": result.p_model_raw,
        "Model_Paired_Smooth": result.p_model_smooth,
        "Deviation_Raw": result.deviation_raw,
        "Deviation_Smooth": result.deviation_smooth,
        "Abs_Deviation_Smooth": np.abs(result.deviation_smooth),
        "Local_GC_Fraction": gc_smooth,
        "Local_AU_Fraction": au_smooth,
    })


# ──────────────────────── Region calling ────────────────────────


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Half-open intervals [s, e) where ``mask`` is True."""
    if mask.size == 0:
        return []
    edges = np.diff(np.concatenate(([0], mask.astype(int), [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    return list(zip(starts.tolist(), ends.tolist()))


def _merge_runs(runs: list[tuple[int, int]], gap: int) -> list[tuple[int, int]]:
    if not runs:
        return runs
    merged = [runs[0]]
    for s, e in runs[1:]:
        ps, pe = merged[-1]
        if s - pe <= int(gap):
            merged[-1] = (ps, e)
        else:
            merged.append((s, e))
    return merged


def _filter_length(runs: list[tuple[int, int]], min_len: int) -> list[tuple[int, int]]:
    return [(s, e) for s, e in runs if (e - s) >= int(min_len)]


def call_regions(deviation_smooth: np.ndarray,
                 *, threshold: float, min_region_length: int,
                 merge_gap: int) -> list[tuple[int, int, str]]:
    """Return [(start_0, end_0, sign), ...] half-open intervals.

    ``sign`` is ``"pos"`` for ``deviation >= +threshold`` and ``"neg"``
    for ``deviation <= -threshold``. Same-sign intervals separated by
    ``<= merge_gap`` are merged before length filtering.
    """
    pos = _runs(deviation_smooth >= float(threshold))
    neg = _runs(deviation_smooth <= -float(threshold))
    pos = _filter_length(_merge_runs(pos, merge_gap), min_region_length)
    neg = _filter_length(_merge_runs(neg, merge_gap), min_region_length)
    out: list[tuple[int, int, str]] = []
    for s, e in pos:
        out.append((s, e, "pos"))
    for s, e in neg:
        out.append((s, e, "neg"))
    out.sort(key=lambda r: r[0])
    return out


# ──────────────────────── Classification ────────────────────────


def classify_region(mean_model: float, mean_dms: float, mean_deviation: float,
                    *, high: float, low: float,
                    dev_thresh: float) -> tuple[str, str]:
    """Return (region_class, interpretation_label).

    Called regions (passed the deviation threshold) end up as
    ``model_high_dms_low``, ``model_low_dms_high``, or ``mixed_deviation``.
    The ``concordant_*`` and ``ambiguous`` classes are reachable when this
    function is reused on bin-level / gene-level summaries that haven't
    been thresholded.
    """
    abs_dev = abs(mean_deviation)
    if abs_dev >= dev_thresh:
        if mean_deviation > 0 and mean_model >= high and mean_dms <= low:
            return ("model_high_dms_low",
                    "Thermodynamically foldable but experimentally open")
        if mean_deviation < 0 and mean_model <= low and mean_dms >= high:
            return ("model_low_dms_high",
                    "Experimentally paired/protected beyond local model")
        return ("mixed_deviation",
                "Intermediate model-DMS deviation; inspect manually")
    if mean_model >= high and mean_dms >= high:
        return ("concordant_paired",
                "Locally encoded paired structural element")
    if mean_model <= low and mean_dms <= low:
        return ("concordant_open",
                "Intrinsically open/unstructured segment")
    return ("ambiguous", "Weak / mixed signal — do not overinterpret")


# ──────────────────────── Null model ────────────────────────


def compute_null_max_stats(
    sequence: str,
    p_dms_smooth: np.ndarray,
    *,
    rolling_window: int,
    rnaplfold_window: int,
    rnaplfold_max_bp_span: int,
    rnaplfold_cutoff: float,
    n_null: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Per-gene null distribution of max-|deviation| under dinucleotide shuffle.

    For each null draw, dinucleotide-shuffle the input sequence (Altschul-
    Erikson, preserves dinucleotide composition), recompute RNAplfold
    P(paired), smooth, and report the maximum absolute deviation against
    the FIXED observed DMS-derived smoothed paired-fraction track.

    The DMS track is treated as the experimental observation and held
    constant across nulls; only the thermodynamic prior is permuted.
    Per-region p-values use the resulting per-gene max distribution as
    the reference (Westfall-Young max-statistic correction), which
    controls family-wise error across regions within a gene without
    needing a separate multiple-testing pass.
    """
    if n_null <= 0:
        return np.empty(0, dtype=float)
    out = np.empty(int(n_null), dtype=float)
    for k in range(int(n_null)):
        shuffled = dinuc_shuffle(sequence, rng)
        res_null = scan_one_gene(
            "_null", "_null", shuffled,
            window=int(rnaplfold_window),
            max_bp_span=int(rnaplfold_max_bp_span),
            cutoff=float(rnaplfold_cutoff),
            dms_structure=None,
        )
        p_model_null_smooth = _centered_rolling_mean(
            res_null.p_paired, int(rolling_window)
        )
        # Same minus operation as the primary deviation track.
        dev_null = p_model_null_smooth - p_dms_smooth
        finite = dev_null[np.isfinite(dev_null)]
        out[k] = float(np.max(np.abs(finite))) if finite.size else float("nan")
    return out


def _empirical_p_one_sided_max(observed: float,
                               null_stats: np.ndarray) -> float:
    """One-sided empirical p: P(null max >= observed) with +1 smoothing."""
    if null_stats.size == 0 or not np.isfinite(observed):
        return float("nan")
    finite = null_stats[np.isfinite(null_stats)]
    if finite.size == 0:
        return float("nan")
    return float((np.sum(finite >= observed) + 1) / (finite.size + 1))


# ──────────────────────── Region annotation ────────────────────────


def _region_majority(transcript_region: list[str], lo: int, hi: int) -> str:
    if hi <= lo or not transcript_region:
        return "unknown"
    span = transcript_region[lo:hi]
    counts: dict[str, int] = {}
    for r in span:
        counts[r] = counts.get(r, 0) + 1
    return max(counts.items(), key=lambda kv: kv[1])[0]


def _composition(sequence: str, lo: int, hi: int) -> dict[str, float]:
    seg = sequence[lo:hi]
    n = len(seg)
    if n == 0:
        return {
            "GC_Fraction": np.nan, "AU_Fraction": np.nan,
            "A_Fraction": np.nan, "U_Fraction": np.nan,
            "G_Fraction": np.nan, "C_Fraction": np.nan,
        }
    a = seg.count("A") / n
    u = seg.count("U") / n
    g = seg.count("G") / n
    c = seg.count("C") / n
    return {
        "GC_Fraction": g + c, "AU_Fraction": a + u,
        "A_Fraction": a, "U_Fraction": u, "G_Fraction": g, "C_Fraction": c,
    }


def _dms_pair_stats(dms_structure: str, lo: int, hi: int) -> dict[str, float]:
    """Mean / median / long-pair fraction for DMS pairs intersecting [lo, hi)."""
    pairs = [
        (i, j) for (i, j) in extract_pairs(dms_structure)
        if (lo <= i < hi) or (lo <= j < hi)
    ]
    if not pairs:
        out: dict[str, float] = dict(
            n_dms_pairs_in_region=0,
            mean_dms_pair_span=np.nan,
            median_dms_pair_span=np.nan,
        )
        for cutoff in _LONG_PAIR_SPAN_CUTOFFS:
            out[f"fraction_dms_pairs_span_gt_{cutoff}"] = np.nan
        return out
    spans = np.asarray([j - i for (i, j) in pairs], dtype=float)
    out = dict(
        n_dms_pairs_in_region=int(len(pairs)),
        mean_dms_pair_span=float(np.mean(spans)),
        median_dms_pair_span=float(np.median(spans)),
    )
    for cutoff in _LONG_PAIR_SPAN_CUTOFFS:
        out[f"fraction_dms_pairs_span_gt_{cutoff}"] = float(np.mean(spans > cutoff))
    return out


def _tis_window_0based(annot: dict | None, n: int,
                       up: int, down: int) -> tuple[int, int] | None:
    if annot is None:
        return None
    cds_start_0 = int(annot["l_utr5"])
    lo = max(0, cds_start_0 - int(up))
    hi = min(n, cds_start_0 + int(down))
    return (lo, hi) if hi > lo else None


def _stop_window_0based(annot: dict | None, n: int,
                        up: int, down: int) -> tuple[int, int] | None:
    if annot is None:
        return None
    cds_end_0 = int(annot["l_utr5"]) + int(annot["l_cds"])  # half-open end
    lo = max(0, cds_end_0 - int(up))
    hi = min(n, cds_end_0 + int(down))
    return (lo, hi) if hi > lo else None


def annotate_regions(result: DeviationResult,
                     regions: list[tuple[int, int, str]],
                     annot: dict | None,
                     transcript_region: list[str],
                     *, cfg: Config,
                     null_max_stats: np.ndarray | None = None) -> pd.DataFrame:
    """Build the per-region table from raw intervals.

    When ``null_max_stats`` is supplied (the per-gene dinucleotide-shuffle
    null-max distribution from :func:`compute_null_max_stats`), each
    region gets a one-sided empirical p-value and a max-statistic-based
    significance label. Without it, the p-value is NaN and the
    ``Statistical_Support_Label`` reads ``"effect_size_only"``.
    """
    n = len(result.sequence)
    if annot is not None:
        cds_start_1 = int(annot["l_utr5"]) + 1
        cds_end_1 = int(annot["l_utr5"]) + int(annot["l_cds"])
    else:
        cds_start_1 = None
        cds_end_1 = None
    tis_win = _tis_window_0based(
        annot, n,
        cfg.structure_deviation_tis_upstream,
        cfg.structure_deviation_tis_downstream,
    )
    stop_win = _stop_window_0based(
        annot, n,
        cfg.structure_deviation_stop_upstream,
        cfg.structure_deviation_stop_downstream,
    )

    rows: list[dict] = []
    for idx, (lo, hi, sign) in enumerate(regions):
        seg_model = result.p_model_smooth[lo:hi]
        seg_dms = result.p_dms_smooth[lo:hi]
        seg_dev = result.deviation_smooth[lo:hi]
        mean_model = float(np.mean(seg_model)) if seg_model.size else float("nan")
        mean_dms = float(np.mean(seg_dms)) if seg_dms.size else float("nan")
        mean_dev = float(np.mean(seg_dev)) if seg_dev.size else float("nan")
        max_abs_dev = float(np.max(np.abs(seg_dev))) if seg_dev.size else float("nan")
        integrated = float(np.sum(seg_dev))
        cls, interp = classify_region(
            mean_model, mean_dms, mean_dev,
            high=cfg.structure_deviation_high_threshold,
            low=cfg.structure_deviation_low_threshold,
            dev_thresh=cfg.structure_deviation_threshold,
        )
        comp = _composition(result.sequence, lo, hi)
        pair_stats = _dms_pair_stats(result.dms_structure, lo, hi)
        midpoint_0 = lo + (hi - lo - 1) / 2.0
        overlaps_tis = (tis_win is not None
                        and not (hi <= tis_win[0] or lo >= tis_win[1]))
        overlaps_stop = (stop_win is not None
                         and not (hi <= stop_win[0] or lo >= stop_win[1]))
        if cds_start_1 is not None:
            d_start_min = (lo + 1) - cds_start_1
            d_start_mid = (midpoint_0 + 1) - cds_start_1
            d_start_max = hi - cds_start_1
            d_stop_mid = (midpoint_0 + 1) - cds_end_1
        else:
            d_start_min = d_start_mid = d_start_max = d_stop_mid = np.nan
        long_pair_cols = {
            f"Fraction_DMS_Pairs_Span_gt_{cutoff}":
                pair_stats[f"fraction_dms_pairs_span_gt_{cutoff}"]
            for cutoff in _LONG_PAIR_SPAN_CUTOFFS
        }
        if null_max_stats is not None and null_max_stats.size > 0:
            empirical_p = _empirical_p_one_sided_max(max_abs_dev, null_max_stats)
            n_null = int(null_max_stats.size)
            null_label = str(getattr(cfg, "structure_deviation_null_model", "dinuc"))
            support_label = (
                "max_stat_significant" if (np.isfinite(empirical_p) and empirical_p < 0.05)
                else "max_stat_nonsignificant"
            )
        else:
            empirical_p = float("nan")
            n_null = 0
            null_label = "none"
            support_label = "effect_size_only"
        rows.append({
            "Species": result.species,
            "Gene": result.gene,
            "Region_ID": f"{result.species}_{result.gene}_{idx + 1:02d}",
            "Start_1based": lo + 1,
            "End_1based": hi,
            "Length_nt": hi - lo,
            "Midpoint_1based": midpoint_0 + 1,
            "Transcript_Region_Majority": _region_majority(transcript_region, lo, hi),
            "Overlaps_TIS_Window": bool(overlaps_tis),
            "Distance_To_Start_Min": d_start_min,
            "Distance_To_Start_Midpoint": d_start_mid,
            "Distance_To_Start_Max": d_start_max,
            "Overlaps_Stop_Window": bool(overlaps_stop),
            "Distance_To_Stop_Midpoint": d_stop_mid,
            "Mean_Model_Paired": mean_model,
            "Mean_DMS_Paired": mean_dms,
            "Mean_Deviation": mean_dev,
            "Abs_Mean_Deviation": abs(mean_dev),
            "Max_Abs_Deviation": max_abs_dev,
            "Integrated_Deviation": integrated,
            "Region_Class": cls,
            "Interpretation_Label": interp,
            "Sign": sign,
            **comp,
            "N_DMS_Pairs_In_Region": pair_stats["n_dms_pairs_in_region"],
            "Mean_DMS_Pair_Span": pair_stats["mean_dms_pair_span"],
            "Median_DMS_Pair_Span": pair_stats["median_dms_pair_span"],
            **long_pair_cols,
            "Empirical_P": empirical_p,
            "Q_Value": np.nan,  # filled in cross-gene at scan_all level
            "Null_Model": null_label,
            "N_Null": n_null,
            "Statistical_Support_Label": support_label,
        })
    return pd.DataFrame(rows)


# ──────────────────────── Gene summary + matrix ────────────────────────


def _bin_intervals(annot: dict, n: int, *, cfg: Config) -> dict[str, tuple[int, int]]:
    """Half-open 0-based intervals for the heatmap bins.

    Bins may overlap or be empty for short transcripts; downstream code
    skips empty bins and lets each bin be reported independently. The
    important property is that each bin's definition is gene-relative
    (not absolute) so multi-gene heatmaps stay comparable.
    """
    cds_start = int(annot["l_utr5"])
    cds_end = cds_start + int(annot["l_cds"])
    tis_up = int(cfg.structure_deviation_tis_upstream)
    tis_down = int(cfg.structure_deviation_tis_downstream)
    stop_up = int(cfg.structure_deviation_stop_upstream)
    stop_down = int(cfg.structure_deviation_stop_downstream)
    early_nt = int(cfg.structure_deviation_early_cds_nt)

    cds_len = max(0, cds_end - cds_start)
    mid_lo_frac = float(getattr(cfg, "structure_deviation_mid_cds_lo", 0.30))
    mid_hi_frac = float(getattr(cfg, "structure_deviation_mid_cds_hi", 0.70))
    late_cds_window = int(getattr(cfg, "structure_deviation_late_cds_window_nt", 300))
    mid_lo = cds_start + int(round(cds_len * mid_lo_frac))
    mid_hi = cds_start + int(round(cds_len * mid_hi_frac))

    bins: dict[str, tuple[int, int]] = {}
    bins["5_end"] = (0, max(0, cds_start - tis_up))
    bins["TIS"] = (max(0, cds_start - tis_up), min(n, cds_start + tis_down))
    bins["early_CDS"] = (
        min(n, cds_start + tis_down),
        min(n, cds_start + early_nt),
    )
    # mid_CDS: configured central fraction of the CDS (default middle 40%).
    bins["mid_CDS"] = (max(0, mid_lo), min(n, mid_hi))
    bins["late_CDS"] = (
        max(0, cds_end - late_cds_window - stop_up),
        max(0, cds_end - stop_up),
    )
    bins["stop_proximal"] = (max(0, cds_end - stop_up), min(n, cds_end + stop_down))
    bins["3_end"] = (min(n, cds_end + stop_down), n)
    return bins


def gene_region_matrix(result: DeviationResult,
                       annot: dict | None,
                       *, cfg: Config) -> pd.DataFrame:
    """One row per (gene, bin) in long format. Empty bins are skipped."""
    n = len(result.sequence)
    if annot is None:
        return pd.DataFrame()
    rows: list[dict] = []
    for name, (lo, hi) in _bin_intervals(annot, n, cfg=cfg).items():
        if hi <= lo:
            continue
        seg_model = result.p_model_smooth[lo:hi]
        seg_dms = result.p_dms_smooth[lo:hi]
        seg_dev = result.deviation_smooth[lo:hi]
        rows.append({
            "Species": result.species,
            "Gene": result.gene,
            "Region_Bin": name,
            "N_Positions": int(hi - lo),
            "Mean_Model_Paired": float(np.mean(seg_model)),
            "Mean_DMS_Paired": float(np.mean(seg_dms)),
            "Mean_Deviation": float(np.mean(seg_dev)),
            "Mean_Abs_Deviation": float(np.mean(np.abs(seg_dev))),
        })
    return pd.DataFrame(rows)


def gene_summary(result: DeviationResult,
                 regions_df: pd.DataFrame,
                 annot: dict | None,
                 *, cfg: Config) -> dict:
    n = len(result.sequence)
    dev = result.deviation_smooth
    abs_dev = np.abs(dev)

    def _class_count(cls: str) -> int:
        if regions_df.empty:
            return 0
        return int((regions_df["Region_Class"] == cls).sum())

    if not regions_df.empty and "Q_Value" in regions_df.columns:
        n_sig = int((regions_df["Q_Value"] < 0.05).fillna(False).sum())
    else:
        n_sig = 0

    summary: dict = {
        "Species": result.species,
        "Gene": result.gene,
        "Length_nt": int(n),
        "Mean_Abs_Deviation": float(np.mean(abs_dev)),
        "Median_Abs_Deviation": float(np.median(abs_dev)),
        "Max_Abs_Deviation": float(np.max(abs_dev)),
        "Deviation_Burden": float(np.mean(abs_dev)),
        "Positive_Deviation_Burden": float(np.mean(np.maximum(dev, 0.0))),
        "Negative_Deviation_Burden": float(np.mean(np.abs(np.minimum(dev, 0.0)))),
        "N_Regions_Total": int(len(regions_df)),
        "N_Significant_Regions_Q05": n_sig,
        "N_Model_High_DMS_Low": _class_count("model_high_dms_low"),
        "N_Model_Low_DMS_High": _class_count("model_low_dms_high"),
        "N_Concordant_Paired": _class_count("concordant_paired"),
        "N_Concordant_Open": _class_count("concordant_open"),
        "N_Mixed_Deviation": _class_count("mixed_deviation"),
        "N_Ambiguous": _class_count("ambiguous"),
    }

    # Region windows for TIS / early-CDS / stop-proximal
    if annot is not None:
        cds_start_0 = int(annot["l_utr5"])
        cds_end_0 = cds_start_0 + int(annot["l_cds"])
        tis_lo = max(0, cds_start_0 - int(cfg.structure_deviation_tis_upstream))
        tis_hi = min(n, cds_start_0 + int(cfg.structure_deviation_tis_downstream))
        ecds_lo = min(n, cds_start_0 + int(cfg.structure_deviation_tis_downstream))
        ecds_hi = min(n, cds_start_0 + int(cfg.structure_deviation_early_cds_nt))
        stop_lo = max(0, cds_end_0 - int(cfg.structure_deviation_stop_upstream))
        stop_hi = min(n, cds_end_0 + int(cfg.structure_deviation_stop_downstream))

        def _seg(arr: np.ndarray, lo: int, hi: int) -> float:
            return float(np.mean(arr[lo:hi])) if hi > lo else float("nan")

        summary.update({
            "TIS_Mean_Model_Paired": _seg(result.p_model_smooth, tis_lo, tis_hi),
            "TIS_Mean_DMS_Paired": _seg(result.p_dms_smooth, tis_lo, tis_hi),
            "TIS_Mean_Deviation": _seg(result.deviation_smooth, tis_lo, tis_hi),
            "TIS_DMS_Openness": (
                1.0 - _seg(result.p_dms_smooth, tis_lo, tis_hi)
                if tis_hi > tis_lo else float("nan")
            ),
            "TIS_Model_Discrepancy": (
                _seg(result.p_model_smooth, tis_lo, tis_hi)
                - _seg(result.p_dms_smooth, tis_lo, tis_hi)
                if tis_hi > tis_lo else float("nan")
            ),
            "Early_CDS_Mean_Model_Paired": _seg(result.p_model_smooth, ecds_lo, ecds_hi),
            "Early_CDS_Mean_DMS_Paired": _seg(result.p_dms_smooth, ecds_lo, ecds_hi),
            "Early_CDS_Mean_Deviation": _seg(result.deviation_smooth, ecds_lo, ecds_hi),
            "Stop_Mean_Model_Paired": _seg(result.p_model_smooth, stop_lo, stop_hi),
            "Stop_Mean_DMS_Paired": _seg(result.p_dms_smooth, stop_lo, stop_hi),
            "Stop_Mean_Deviation": _seg(result.deviation_smooth, stop_lo, stop_hi),
        })

        # Long-range pair-span metrics on the DMS structure
        dms_pairs = extract_pairs(result.dms_structure)
        if dms_pairs:
            spans = np.asarray([j - i for (i, j) in dms_pairs], dtype=float)
            for cutoff in _LONG_PAIR_SPAN_CUTOFFS:
                summary[f"DMS_Fraction_Pairs_gt_{cutoff}"] = float(
                    np.mean(spans > cutoff)
                )
            summary["DMS_Mean_Pair_Span"] = float(np.mean(spans))
        else:
            for cutoff in _LONG_PAIR_SPAN_CUTOFFS:
                summary[f"DMS_Fraction_Pairs_gt_{cutoff}"] = float("nan")
            summary["DMS_Mean_Pair_Span"] = float("nan")

    return summary


# ──────────────────────── Driver ────────────────────────


def scan_all(cfg: Config) -> dict:
    """Walk every (species, gene) target and emit all dataframes plus the
    raw per-gene results (for plotting). Returns a dict keyed by output
    name so the command layer just consumes it.

    When ``cfg.structure_deviation_null_model == "dinuc"`` and
    ``cfg.structure_deviation_n_null > 0``, a per-gene dinucleotide-
    shuffle null distribution of max-|deviation| is built (Westfall-Young
    max-statistic correction within gene), and a Benjamini-Hochberg pass
    converts the resulting per-region empirical p-values into ``Q_Value``
    across all genes pooled. Defaults are conservative
    (``null_model: "none"``, ``n_null: 0``); enable by setting both in
    your YAML config.
    """
    targets = {canonical_gene(g) for g in (cfg.target_genes or ())}
    null_model = str(getattr(cfg, "structure_deviation_null_model", "none"))
    n_null = int(getattr(cfg, "structure_deviation_n_null", 0))
    use_null = null_model == "dinuc" and n_null > 0
    null_rng = make_rng(int(getattr(cfg, "seed", 42)) + 17) if use_null else None
    step("structure-deviation: per-gene scan"
         + (f" (with dinuc null, n={n_null})" if use_null else ""))
    results: list[DeviationResult] = []
    per_pos_frames: list[pd.DataFrame] = []
    region_frames: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    matrix_frames: list[pd.DataFrame] = []
    for species, fname in cfg.db_files.items():
        recs = [r for r in parse_db(cfg.data_dir / fname)
                if (not targets or canonical_gene(r.gene) in targets)]
        for rec in progress(recs, desc=f"{species} deviation", unit="gene"):
            try:
                annot = annotation_for(species, rec.gene)
            except KeyError:
                annot = None
            result = compute_one_gene(
                species, rec.gene, rec.sequence, rec.structure, cfg=cfg,
            )
            results.append(result)
            pp = per_position_table(result, annot)
            per_pos_frames.append(pp)
            transcript_region = pp["Transcript_Region"].tolist()
            regions = call_regions(
                result.deviation_smooth,
                threshold=cfg.structure_deviation_threshold,
                min_region_length=cfg.structure_deviation_min_region_length,
                merge_gap=cfg.structure_deviation_merge_gap,
            )
            null_max_stats: np.ndarray | None = None
            if use_null and null_rng is not None and len(result.sequence) >= 3:
                null_max_stats = compute_null_max_stats(
                    result.sequence,
                    result.p_dms_smooth,
                    rolling_window=result.rolling_window,
                    rnaplfold_window=result.rnaplfold_window,
                    rnaplfold_max_bp_span=result.rnaplfold_max_bp_span,
                    rnaplfold_cutoff=result.rnaplfold_cutoff,
                    n_null=n_null,
                    rng=null_rng,
                )
            reg_df = annotate_regions(
                result, regions, annot, transcript_region, cfg=cfg,
                null_max_stats=null_max_stats,
            )
            region_frames.append(reg_df)
            summary_rows.append(gene_summary(result, reg_df, annot, cfg=cfg))
            matrix_frames.append(gene_region_matrix(result, annot, cfg=cfg))

    regions_all = (pd.concat(region_frames, ignore_index=True)
                   if region_frames else pd.DataFrame())
    # Cross-gene BH-FDR over the per-region empirical p-values. Each region
    # already carries a max-stat-corrected p-value within its gene, so this
    # second-stage BH controls the false-discovery rate across the global
    # set of called regions reported in one run.
    if not regions_all.empty and "Empirical_P" in regions_all.columns:
        regions_all["Q_Value"] = bh_fdr(
            regions_all["Empirical_P"].to_numpy(dtype=float)
        )

    return {
        "per_position": (pd.concat(per_pos_frames, ignore_index=True)
                         if per_pos_frames else pd.DataFrame()),
        "regions": regions_all,
        "gene_summary": pd.DataFrame(summary_rows),
        "gene_region_matrix": (pd.concat(matrix_frames, ignore_index=True)
                               if matrix_frames else pd.DataFrame()),
        "results": results,
    }
