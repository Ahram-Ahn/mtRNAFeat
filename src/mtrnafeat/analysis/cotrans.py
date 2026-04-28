"""Per-gene cotranscriptional / sliding-window signal scan.

For each gene, fold a series of subsequences and extract Vienna's structural
signals so we can spot positions where the structure rearranges sharply —
the kind of "folding event" the legacy ``signal_analysis.py`` script chased.

Two scan modes:

* ``"prefix"`` — the cotranscriptional one. Fold ``seq[0:L]`` for
  ``L = step, 2·step, …, len(seq)``. Each prefix represents the nascent
  RNA at a given transcription length. Useful for asking "when does the
  emerging structure commit?"
* ``"sliding"`` — fixed-size windows of length ``window`` at step intervals.
  Faster (every fold is small); equivalent to ScanFold's per-window MFE but
  enriched with diversity + paired-fraction tracks.

For each sample we record:

* ``MFE`` and ``MFE_per_nt``
* ``Diversity`` — ViennaRNA's mean ensemble base-pair distance, an
  uncertainty proxy. Lower → the ensemble is committed to one structure.
* ``Paired_Fraction`` — fraction of nt paired in the MFE structure.
* ``BP_Dist_To_Prev`` — base-pair distance to the previous sample's
  structure (only valid when both samples have the same length, i.e. in
  ``sliding`` mode; ``NaN`` for ``prefix``).

Then ``Delta_*`` columns hold first differences. ``z_score`` and
``find_signal_peaks`` work on those.

This module deliberately doesn't fold or compare against the .db DMS
structure — that's the job of the ``window`` command. Here we're purely
asking "what does the engine think this sequence wants to do, and where
does its mind change?"
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.structure import paired_fraction
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step

VALID_MODES = ("prefix", "sliding")


@dataclass
class ScanParams:
    mode: str = "sliding"
    window: int = 120          # only used in sliding mode
    step: int = 30
    min_len: int = 30          # skip prefixes / windows shorter than this
    z_threshold: float = 2.0
    smooth_window: int = 5     # rolling mean width for delta-signal smoothing

    def validate(self) -> None:
        if self.mode not in VALID_MODES:
            raise ValueError(f"mode must be one of {VALID_MODES}; got {self.mode!r}")
        if self.step <= 0 or self.window <= 0 or self.min_len <= 0:
            raise ValueError("step, window, min_len must be positive")


def _scan_intervals(n: int, params: ScanParams) -> list[tuple[int, int]]:
    """Return half-open (start, end) windows in 0-based coords."""
    if params.mode == "prefix":
        ends = list(range(params.step, n + 1, params.step))
        if ends and ends[-1] < n:
            ends.append(n)
        return [(0, e) for e in ends if e >= params.min_len]
    # sliding
    if params.window > n:
        return [(0, n)] if n >= params.min_len else []
    starts = list(range(0, n - params.window + 1, params.step))
    out = [(s, s + params.window) for s in starts]
    if out and out[-1][1] < n:
        out.append((n - params.window, n))
    return out


def _structures_comparable(a: str, b: str) -> bool:
    return len(a) == len(b)


def gene_signals(species: str, gene: str, sequence: str,
                  params: ScanParams,
                  rec_label: str | None = None) -> pd.DataFrame:
    """Slide / grow over one gene and emit the signal DataFrame."""
    params.validate()
    intervals = _scan_intervals(len(sequence), params)
    rows: list[dict] = []
    prev_struct: str | None = None
    label = rec_label or f"{species} {gene}"
    bar = progress(intervals, desc=label, unit="win", leave=False)
    for s, e in bar:
        sub = sequence[s:e]
        if len(sub) < params.min_len:
            continue
        struct, mfe = thermo.fold_mfe(sub)
        try:
            div = thermo.ensemble_diversity(sub)
        except Exception:
            div = float("nan")
        pf = paired_fraction(struct)
        bp_dist = float("nan")
        if prev_struct is not None and _structures_comparable(prev_struct, struct):
            try:
                bp_dist = float(thermo.bp_distance(prev_struct, struct))
            except Exception:
                bp_dist = float("nan")
        center = s + 1 + (len(sub) - 1) / 2.0
        rows.append({
            "Species": species,
            "Gene": canonical_gene(gene),
            "Mode": params.mode,
            "Window_Start_1based": s + 1,
            "Window_End_1based": e,
            "Window_Length_nt": len(sub),
            "Window_Center_1based": center,
            "MFE": float(mfe),
            "MFE_per_nt": float(mfe) / len(sub),
            "Diversity": float(div),
            "Paired_Fraction": float(pf),
            "BP_Dist_To_Prev": bp_dist,
        })
        prev_struct = struct

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # First differences for the "event detection" tracks. For sliding mode
    # the windows are equal-length so deltas are clean; for prefix mode the
    # diff is between prefixes so MFE diff has a length-bias — we use MFE/nt
    # which mostly cancels that out.
    df["Delta_MFE_per_nt"] = df["MFE_per_nt"].diff().fillna(0.0)
    df["Delta_Diversity"] = df["Diversity"].diff().fillna(0.0)
    df["Delta_Paired_Fraction"] = df["Paired_Fraction"].diff().fillna(0.0)

    # Smooth the delta signals before z-scoring so single-window noise doesn't
    # generate spurious peaks.
    w = max(1, int(params.smooth_window))
    for col in ("Delta_MFE_per_nt", "Delta_Diversity", "Delta_Paired_Fraction"):
        df[f"{col}_Smooth"] = df[col].rolling(w, center=True, min_periods=1).mean()

    return df


def z_score(values: np.ndarray | pd.Series) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    mu = float(np.nanmean(arr))
    sd = float(np.nanstd(arr, ddof=1)) if arr.size > 1 else 0.0
    if sd == 0.0 or not np.isfinite(sd):
        return np.zeros_like(arr)
    return (arr - mu) / sd


def add_z_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Append Z_<col> columns for every smoothed delta signal."""
    if df.empty:
        return df
    df = df.copy()
    for col in ("Delta_MFE_per_nt_Smooth", "Delta_Diversity_Smooth",
                 "Delta_Paired_Fraction_Smooth"):
        if col in df.columns:
            df[f"Z_{col}"] = z_score(df[col].values)
    return df


def find_signal_peaks(df: pd.DataFrame, threshold: float = 2.0,
                      distance: int = 5) -> dict[str, list[int]]:
    """Locate indices where each smoothed Z-signal exceeds ±threshold.

    Returns a dict: signal name → list of row indices that are peaks.
    Both directions are tracked: very negative MFE deltas (sudden stabilization),
    very positive diversity deltas (sudden uncertainty), etc.
    """
    out: dict[str, list[int]] = {}
    if df.empty:
        return out
    for col, sign in (
        ("Z_Delta_MFE_per_nt_Smooth", -1),       # very negative = big stabilization
        ("Z_Delta_Diversity_Smooth", -1),         # very negative = sudden commitment
        ("Z_Delta_Paired_Fraction_Smooth", 1),    # very positive = pairing surge
    ):
        if col not in df.columns:
            continue
        z = df[col].to_numpy()
        if sign < 0:
            peaks, _ = find_peaks(-z, height=threshold, distance=distance)
        else:
            peaks, _ = find_peaks(z, height=threshold, distance=distance)
        out[col] = peaks.tolist()
    return out


def per_gene_cotrans_scan(cfg: Config, params: ScanParams) -> pd.DataFrame:
    """Walk every (species, gene) in cfg.target_genes and emit a long-form
    DataFrame with all per-window signals + z-scores.
    """
    params.validate()
    targets = {canonical_gene(g) for g in (cfg.target_genes or ())}
    step(f"cotrans scan (mode={params.mode}, step={params.step}"
         + (f", window={params.window}" if params.mode == "sliding" else "")
         + ")")
    frames: list[pd.DataFrame] = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        recs = [r for r in parse_db(path)
                if (not targets or canonical_gene(r.gene) in targets)]
        for rec in progress(recs, desc=f"{species} cotrans", unit="gene"):
            sig = gene_signals(species, rec.gene, rec.sequence, params)
            if sig.empty:
                continue
            sig = add_z_columns(sig)
            frames.append(sig)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)
