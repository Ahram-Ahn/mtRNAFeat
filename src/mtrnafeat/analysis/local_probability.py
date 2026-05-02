"""Per-position local pair-probability scan via ViennaRNA's RNAplfold,
with optional overlay against the DMS-derived dot-bracket from the matching
.db record.

For each transcript:

* ``p_paired[i]`` — Boltzmann probability that position ``i`` is paired
  with *some* partner within ``max_bp_span`` nt, averaged across all
  overlapping windows of length ``window``. RNAplfold is sequence-only;
  DMS reactivity is **not** fed in (re-fitting DMS into RNAplfold would
  re-introduce the pseudo-energy step we're trying to avoid; that path is
  RNAstructure's ``Fold -dms``, already used by the ``window`` command).
* ``p_unpaired_u1[i] = 1 − p_paired[i]`` — probability that position ``i``
  is single-stranded.
* When the matching .db record is available (and length matches the
  RNAplfold input), three additional comparisons are produced:
  per-position ``DMS_Paired_Binary`` / smoothed track, sliding-window
  agreement (mean RNAplfold P(paired) vs DMS paired fraction), and a
  TIS-vs-CDS-background summary with a circular-shift empirical p-value
  (preserves autocorrelation without refolding).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import RNA  # type: ignore
except ImportError as e:
    RNA = None  # type: ignore
    _RNA_IMPORT_ERROR: ImportError | None = e
else:
    _RNA_IMPORT_ERROR = None

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core.structure import extract_pairs, pair_table
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import DbRecord, parse_db
from mtrnafeat.progress import progress, step


def _require_rna() -> None:
    if RNA is None:
        raise ImportError(
            "ViennaRNA Python bindings not installed; "
            "needed for local_probability. "
            f"(original: {_RNA_IMPORT_ERROR})"
        )


@dataclass
class LocalProbResult:
    species: str
    gene: str
    sequence: str
    p_paired: np.ndarray   # length N, in [0, 1]
    window: int
    max_bp_span: int
    cutoff: float
    # DMS-derived companion data (optional). Populated by ``scan_all`` when
    # the matching .db record's length agrees with the RNAplfold input;
    # left as ``None`` when no .db record was found or the lengths
    # disagree. ``dms_paired_binary`` is the materialized 0/1 vector used
    # downstream so consumers don't re-parse the dot-bracket.
    dms_structure: str | None = None
    dms_paired_binary: np.ndarray | None = None


def paired_binary_from_dotbracket(structure: str) -> np.ndarray:
    """0/1 per position: 1 if the position is paired in the dot-bracket."""
    pt = pair_table(structure)
    arr = np.zeros(len(structure), dtype=np.int8)
    for i in pt:
        arr[i] = 1
    return arr


def _centered_rolling_mean(values: np.ndarray, width: int) -> np.ndarray:
    """Length-preserving centered rolling mean with min_periods=1 semantics."""
    if width <= 1:
        return values.astype(float, copy=True)
    return (
        pd.Series(values, dtype=float)
        .rolling(int(width), center=True, min_periods=1)
        .mean()
        .to_numpy()
    )


def scan_one_gene(species: str, gene: str, sequence: str,
                   window: int, max_bp_span: int, cutoff: float,
                   dms_structure: str | None = None) -> LocalProbResult:
    """Run RNAplfold-style folding on one transcript and return p_paired track.

    If ``dms_structure`` is given and its length matches ``sequence``, the
    materialized DMS paired binary is attached to the result. Mismatched
    lengths drop the DMS overlay silently (the caller should already have
    logged the mismatch).
    """
    _require_rna()
    n = len(sequence)
    plist = RNA.pfl_fold(sequence, int(window), int(max_bp_span), float(cutoff))
    p_paired = np.zeros(n, dtype=float)
    # plist entries are 1-based (i, j, p); each pair touches both i and j.
    for entry in plist:
        i, j, p = int(entry.i), int(entry.j), float(entry.p)
        if 1 <= i <= n:
            p_paired[i - 1] += p
        if 1 <= j <= n and j != i:
            p_paired[j - 1] += p
    # Probabilities can drift fractionally above 1 in edge cases; clip.
    p_paired = np.clip(p_paired, 0.0, 1.0)
    dms_binary = None
    if dms_structure is not None and len(dms_structure) == n:
        dms_binary = paired_binary_from_dotbracket(dms_structure)
    return LocalProbResult(
        species=species,
        gene=canonical_gene(gene),
        sequence=sequence,
        p_paired=p_paired,
        window=int(window),
        max_bp_span=int(max_bp_span),
        cutoff=float(cutoff),
        dms_structure=dms_structure if dms_binary is not None else None,
        dms_paired_binary=dms_binary,
    )


def per_position_table(result: LocalProbResult,
                       annot: dict | None = None,
                       *, smooth: int = 25,
                       tis_upstream: int = 50,
                       tis_downstream: int = 50) -> pd.DataFrame:
    """Per-nt DataFrame with paired / unpaired probabilities, region tag,
    and (when DMS is available) the DMS-paired binary, smoothed tracks,
    Δ track, and an ``In_TIS_Window`` boolean."""
    n = len(result.sequence)
    if annot is not None:
        cds_start_1 = annot["l_utr5"] + 1
        cds_end_1 = annot["l_utr5"] + annot["l_cds"]
        region: list[str] = []
        for pos in range(1, n + 1):
            if pos < cds_start_1:
                region.append("5'UTR")
            elif pos <= cds_end_1:
                region.append("CDS")
            else:
                region.append("3'UTR")
    else:
        region = [""] * n

    p = result.p_paired
    p_smooth = _centered_rolling_mean(p, smooth)
    cols: dict[str, object] = {
        "Species": result.species,
        "Gene": result.gene,
        "Position_1based": np.arange(1, n + 1),
        "Nucleotide": list(result.sequence),
        "Region": region,
        "P_Paired": p,
        "P_Paired_Smoothed": p_smooth,
        "P_Unpaired_u1": 1.0 - p,
        "Window_nt": result.window,
        "Max_BP_Span_nt": result.max_bp_span,
        "Cutoff": result.cutoff,
    }

    if result.dms_paired_binary is not None:
        dms_bin = result.dms_paired_binary.astype(float)
        dms_smooth = _centered_rolling_mean(dms_bin, smooth)
        cols["DMS_Paired_Binary"] = dms_bin.astype(np.int8)
        cols["DMS_Paired_Smoothed"] = dms_smooth
        cols["Delta_Ppaired_minus_DMS_Smoothed"] = p_smooth - dms_smooth
    else:
        nan_arr = np.full(n, np.nan, dtype=float)
        cols["DMS_Paired_Binary"] = nan_arr
        cols["DMS_Paired_Smoothed"] = nan_arr
        cols["Delta_Ppaired_minus_DMS_Smoothed"] = nan_arr

    if annot is not None:
        cds_start_0 = int(annot["l_utr5"])
        tis_lo = max(0, cds_start_0 - int(tis_upstream))
        tis_hi = min(n, cds_start_0 + int(tis_downstream))
        in_tis = np.zeros(n, dtype=bool)
        if tis_hi > tis_lo:
            in_tis[tis_lo:tis_hi] = True
        cols["In_TIS_Window"] = in_tis
    else:
        cols["In_TIS_Window"] = np.zeros(n, dtype=bool)

    return pd.DataFrame(cols)


def _region_majority(region_strs: list[str], lo: int, hi: int) -> str:
    if not region_strs or hi <= lo:
        return ""
    # lo/hi are 0-based; region_strs is parallel to positions.
    span = region_strs[lo:hi]
    if not span:
        return ""
    counts: dict[str, int] = {}
    for r in span:
        counts[r] = counts.get(r, 0) + 1
    return max(counts.items(), key=lambda kv: kv[1])[0]


def per_window_agreement_table(result: LocalProbResult,
                               annot: dict | None,
                               *, win: int, step_nt: int) -> pd.DataFrame:
    """Sliding-window mean P(paired) and DMS paired fraction per gene.

    Windows are 0-based half-open ``[lo, hi)`` of length ``win`` at stride
    ``step_nt``; a final tail window is appended if the regular grid
    leaves the last position uncovered (mirrors ``analysis.cotrans``)."""
    n = len(result.sequence)
    if win <= 0 or step_nt <= 0:
        raise ValueError("win and step_nt must be positive")
    if win > n:
        intervals = [(0, n)] if n > 0 else []
    else:
        intervals = [(s, s + win) for s in range(0, n - win + 1, step_nt)]
        if intervals and intervals[-1][1] < n:
            intervals.append((n - win, n))

    # Region per position (1-based labels rendered into a parallel list)
    region_strs: list[str] = []
    if annot is not None:
        cds_start_1 = annot["l_utr5"] + 1
        cds_end_1 = annot["l_utr5"] + annot["l_cds"]
        for pos in range(1, n + 1):
            if pos < cds_start_1:
                region_strs.append("5'UTR")
            elif pos <= cds_end_1:
                region_strs.append("CDS")
            else:
                region_strs.append("3'UTR")

    p = result.p_paired
    dms_bin = result.dms_paired_binary
    pairs_dms = extract_pairs(result.dms_structure) if result.dms_structure else []

    rows: list[dict] = []
    for lo, hi in intervals:
        seg = p[lo:hi]
        rnaplfold_mean = float(np.mean(seg)) if seg.size else float("nan")
        rnaplfold_median = float(np.median(seg)) if seg.size else float("nan")
        rnaplfold_low_frac = float(np.mean(seg < 0.25)) if seg.size else float("nan")

        if dms_bin is not None:
            dms_seg = dms_bin[lo:hi]
            dms_pf = float(np.mean(dms_seg)) if dms_seg.size else float("nan")
            inwin_pairs = [(i, j) for (i, j) in pairs_dms
                           if lo <= i < hi and lo <= j < hi]
            if inwin_pairs:
                spans = [j - i for (i, j) in inwin_pairs]
                dms_mean_span = float(np.mean(spans))
                dms_long = float(
                    sum(1 for s in spans if s > 50) / len(spans)
                )
            else:
                dms_mean_span = float("nan")
                dms_long = float("nan")
            agreement_signed = rnaplfold_mean - dms_pf
            agreement_abs = abs(agreement_signed)
        else:
            dms_pf = float("nan")
            dms_mean_span = float("nan")
            dms_long = float("nan")
            agreement_signed = float("nan")
            agreement_abs = float("nan")

        rows.append({
            "Species": result.species,
            "Gene": result.gene,
            "Window_Start_1based": lo + 1,
            "Window_End_1based": hi,
            "Window_Center_1based": lo + 1 + (hi - lo - 1) / 2.0,
            "Window_Size": hi - lo,
            "Region_Majority": _region_majority(region_strs, lo, hi),
            "RNAplfold_Mean_Ppaired": rnaplfold_mean,
            "RNAplfold_Median_Ppaired": rnaplfold_median,
            "RNAplfold_Frac_Low_Ppaired_0p25": rnaplfold_low_frac,
            "DMS_Paired_Fraction": dms_pf,
            "DMS_Mean_Pair_Span": dms_mean_span,
            "DMS_Long_Pair_Fraction_gt_50": dms_long,
            "Agreement_Signed_Delta": agreement_signed,
            "Agreement_Abs_Delta": agreement_abs,
        })

    return pd.DataFrame(rows)


def _classify_concordance(rnaplfold: float, dms: float,
                          *, low: float = 0.25, high: float = 0.5) -> str:
    if not (np.isfinite(rnaplfold) and np.isfinite(dms)):
        return "unknown"
    rl = rnaplfold < low
    dl = dms < low
    rh = rnaplfold > high
    dh = dms > high
    if rl and dl:
        return "concordant_open"
    if rh and dh:
        return "concordant_paired"
    if dl and rh:
        return "DMS_open_only"
    if rl and dh:
        return "RNAplfold_open_only"
    return "discordant"


def _circular_shift_p_low(track: np.ndarray, lo: int, hi: int,
                          n_shifts: int, rng: np.random.Generator) -> float:
    """Empirical p-value: P(mean(shifted track[lo:hi]) ≤ observed mean).

    A circular shift preserves the autocorrelation structure of the track,
    which is what we want — RNAplfold P(paired) and the DMS paired binary
    are heavily autocorrelated, so an i.i.d. permutation null would
    overstate significance. Returns NaN when the input track is empty,
    constant, or the window is degenerate.
    """
    n = track.size
    if n == 0 or hi <= lo or n_shifts <= 0:
        return float("nan")
    finite = track[np.isfinite(track)]
    if finite.size == 0:
        return float("nan")
    observed = float(np.nanmean(track[lo:hi]))
    if not np.isfinite(observed):
        return float("nan")
    shifts = rng.integers(low=1, high=n, size=int(n_shifts))
    le = 0
    for k in shifts:
        rolled = np.roll(track, int(k))
        m = float(np.nanmean(rolled[lo:hi]))
        if np.isfinite(m) and m <= observed:
            le += 1
    return (le + 1) / (int(n_shifts) + 1)


def tis_summary_row(result: LocalProbResult,
                    annot: dict,
                    *, upstream: int, downstream: int,
                    n_circ_shifts: int,
                    rng: np.random.Generator) -> dict:
    """One-row TIS summary: TIS vs CDS-background effect size + circular-
    shift empirical p-values for "is the start region unusually open?".

    Background: same-size windows whose center lies in the CDS but at
    least 100 nt away from the start codon (matches the report's spec).
    """
    n = len(result.sequence)
    cds_start_0 = int(annot["l_utr5"])
    l_cds = int(annot["l_cds"])
    win_size = int(upstream) + int(downstream)
    upstream_avail = min(int(upstream), cds_start_0)
    has_full_upstream = upstream_avail == int(upstream)
    tis_lo = max(0, cds_start_0 - int(upstream))
    tis_hi = min(n, cds_start_0 + int(downstream))

    p = result.p_paired
    dms_bin = (result.dms_paired_binary.astype(float)
               if result.dms_paired_binary is not None else None)

    tis_rnaplfold = (float(np.nanmean(p[tis_lo:tis_hi]))
                     if tis_hi > tis_lo else float("nan"))
    tis_dms = (float(np.nanmean(dms_bin[tis_lo:tis_hi]))
               if (dms_bin is not None and tis_hi > tis_lo) else float("nan"))

    # CDS background: windows of size win_size whose CENTER lies in the
    # CDS interior, excluding ±100 nt around the start codon.
    cds_lo = cds_start_0
    cds_hi = cds_start_0 + l_cds
    bg_means_rnaplfold: list[float] = []
    bg_means_dms: list[float] = []
    if win_size > 0 and cds_hi - cds_lo >= win_size:
        for s in range(cds_lo, cds_hi - win_size + 1):
            center = s + win_size / 2.0
            if abs(center - cds_start_0) < 100:
                continue
            seg_p = p[s:s + win_size]
            bg_means_rnaplfold.append(float(np.nanmean(seg_p)))
            if dms_bin is not None:
                bg_means_dms.append(float(np.nanmean(dms_bin[s:s + win_size])))
    cds_bg_rnaplfold = (float(np.nanmean(bg_means_rnaplfold))
                        if bg_means_rnaplfold else float("nan"))
    cds_bg_dms = (float(np.nanmean(bg_means_dms))
                  if bg_means_dms else float("nan"))

    p_emp_rnaplfold = _circular_shift_p_low(
        p, tis_lo, tis_hi, int(n_circ_shifts), rng)
    p_emp_dms = (
        _circular_shift_p_low(dms_bin, tis_lo, tis_hi,
                              int(n_circ_shifts), rng)
        if dms_bin is not None else float("nan")
    )

    return {
        "Species": result.species,
        "Gene": result.gene,
        "TIS_Upstream_Nt_Requested": int(upstream),
        "TIS_Downstream_Nt_Requested": int(downstream),
        "TIS_Upstream_Nt_Available": int(upstream_avail),
        "Has_Full_Upstream_Context": bool(has_full_upstream),
        "TIS_Window_Start_1based": int(tis_lo + 1),
        "TIS_Window_End_1based": int(tis_hi),
        "TIS_RNAplfold_Mean_Ppaired": tis_rnaplfold,
        "TIS_DMS_Paired_Fraction": tis_dms,
        "CDS_Background_RNAplfold_Mean_Ppaired": cds_bg_rnaplfold,
        "CDS_Background_DMS_Paired_Fraction": cds_bg_dms,
        "TIS_vs_CDS_Delta_RNAplfold": tis_rnaplfold - cds_bg_rnaplfold
            if np.isfinite(tis_rnaplfold) and np.isfinite(cds_bg_rnaplfold)
            else float("nan"),
        "TIS_vs_CDS_Delta_DMS": tis_dms - cds_bg_dms
            if np.isfinite(tis_dms) and np.isfinite(cds_bg_dms)
            else float("nan"),
        "Empirical_P_TIS_Low_Ppaired_vs_CDS_Windows": p_emp_rnaplfold,
        "Empirical_P_TIS_Low_DMS_Pairfrac_vs_CDS_Windows": p_emp_dms,
        "TIS_Concordance_Class": _classify_concordance(tis_rnaplfold, tis_dms),
        "N_Circular_Shifts": int(n_circ_shifts),
    }


def scan_all(cfg: Config, window: int, max_bp_span: int,
             cutoff: float,
             *, smooth: int | None = None,
             tis_upstream: int | None = None,
             tis_downstream: int | None = None
             ) -> tuple[pd.DataFrame, list[LocalProbResult]]:
    """Walk every (species, gene) in cfg.target_genes and emit a long-form
    per-position DataFrame plus the raw per-gene results (for plotting).

    The matching .db record's dot-bracket is attached to each result when
    its length agrees with the RNAplfold input, which unlocks the DMS
    overlay columns in ``per_position_table``."""
    targets = {canonical_gene(g) for g in (cfg.target_genes or ())}
    smooth_w = int(smooth if smooth is not None else cfg.rolling_window)
    tis_up = int(tis_upstream if tis_upstream is not None else cfg.tis_upstream_nt)
    tis_dn = int(tis_downstream if tis_downstream is not None else cfg.tis_downstream_nt)
    step(f"local-probability scan (W={window}, L={max_bp_span}, cutoff={cutoff})")
    results: list[LocalProbResult] = []
    for species, fname in cfg.db_files.items():
        recs = [r for r in parse_db(cfg.data_dir / fname)
                if (not targets or canonical_gene(r.gene) in targets)]
        for rec in progress(recs, desc=f"{species} plfold", unit="gene"):
            res = scan_one_gene(
                species, rec.gene, rec.sequence,
                window=window, max_bp_span=max_bp_span, cutoff=cutoff,
                dms_structure=rec.structure,
            )
            results.append(res)
    if not results:
        return pd.DataFrame(), []
    frames = []
    for res in results:
        try:
            annot = annotation_for(res.species, res.gene)
        except KeyError:
            annot = None
        frames.append(per_position_table(
            res, annot,
            smooth=smooth_w,
            tis_upstream=tis_up,
            tis_downstream=tis_dn,
        ))
    return pd.concat(frames, ignore_index=True), results
