"""Yeast↔human COX1 codon-aligned comparative analysis.

Reads the PAL2NAL alignment, classifies every column as identical /
synonymous / non-synonymous, locates substitutions by codon position
(1st / 2nd / wobble), and projects local Vienna ΔG difference onto the
codon frame.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.core import thermo
from mtrnafeat.io.alignment import CodonAlignment, codon_position_changes, parse_pal2nal
from mtrnafeat.io.codons import HUMAN_MT_CODON_TABLE, YEAST_MT_CODON_TABLE
from mtrnafeat.io.db_parser import get_record


def classify_column(codon_y: str, codon_h: str) -> dict:
    aa_y = YEAST_MT_CODON_TABLE.get(codon_y, "?")
    aa_h = HUMAN_MT_CODON_TABLE.get(codon_h, "?")
    same_aa = aa_y == aa_h
    diffs = codon_position_changes(codon_y, codon_h)
    return {
        "yeast_codon": codon_y,
        "human_codon": codon_h,
        "yeast_aa": aa_y,
        "human_aa": aa_h,
        "same_aa": bool(same_aa),
        "n_diffs": len(diffs),
        "diff_positions": ",".join(str(p) for p in diffs),
        "is_synonymous": bool(same_aa and len(diffs) > 0),
        "is_nonsynonymous": bool((not same_aa) and len(diffs) > 0 and "-" not in (codon_y + codon_h)),
        "is_gap": bool("-" in (codon_y + codon_h)),
    }


def alignment_table(cfg: Config) -> pd.DataFrame:
    alignment = parse_pal2nal(cfg.data_dir / cfg.alignment_file)
    rows = []
    for col, (cy, ch, ay, ah) in enumerate(zip(
        alignment.yeast_codons, alignment.human_codons,
        alignment.yeast_aa, alignment.human_aa,
    ), start=1):
        rows.append({"col": col, **classify_column(cy, ch)})
    return pd.DataFrame(rows)


def substitution_summary(table: pd.DataFrame) -> pd.DataFrame:
    """Per-position substitution direction tally (e.g. A->C at position 1)."""
    bases = list("ACGTU")
    rows = []
    for _, r in table.iterrows():
        if r["is_gap"] or r["n_diffs"] == 0:
            continue
        cy = r["yeast_codon"]; ch = r["human_codon"]
        for pos in [int(p) for p in r["diff_positions"].split(",") if p]:
            rows.append({
                "Position": pos,
                "Yeast_Base": cy[pos - 1],
                "Human_Base": ch[pos - 1],
                "Same_AA": r["same_aa"],
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.groupby(["Position", "Yeast_Base", "Human_Base", "Same_AA"]).size().reset_index(name="Count")


def cox1_local_dG_track(cfg: Config, window: int = 60) -> pd.DataFrame:
    """For each codon column, compute Vienna local ΔG over a +/- window/2-nt
    span around the column's nucleotide position in each species. Also
    reports the difference (yeast − human)."""
    table = alignment_table(cfg)
    rec_y = get_record(cfg.data_dir / cfg.db_files["Yeast"], "COX1")
    rec_h = get_record(cfg.data_dir / cfg.db_files["Human"], "COX1")

    annot_y = _cox1_offsets(cfg, "Yeast")
    annot_h = _cox1_offsets(cfg, "Human")

    seq_y = rec_y.sequence
    seq_h = rec_h.sequence

    pos_y = annot_y  # 0-based offset into yeast COX1 sequence at start of CDS
    pos_h = annot_h
    half = window // 2
    rows = []
    for _, r in table.iterrows():
        col = int(r["col"])
        if r["is_gap"]:
            rows.append({"col": col, "yeast_dG": np.nan, "human_dG": np.nan, "delta_dG_yeast_minus_human": np.nan})
            continue
        ny = pos_y + (col - 1) * 3
        nh = pos_h + (col - 1) * 3
        yseg = _safe_slice(seq_y, ny - half, ny + half)
        hseg = _safe_slice(seq_h, nh - half, nh + half)
        try:
            _, dg_y = thermo.fold_mfe(yseg) if yseg else (None, np.nan)
        except Exception:
            dg_y = np.nan
        try:
            _, dg_h = thermo.fold_mfe(hseg) if hseg else (None, np.nan)
        except Exception:
            dg_h = np.nan
        rows.append({"col": col, "yeast_dG": dg_y, "human_dG": dg_h,
                       "delta_dG_yeast_minus_human": dg_y - dg_h if (yseg and hseg) else np.nan})
    return pd.DataFrame(rows)


def _safe_slice(seq: str, start: int, end: int) -> str:
    if start < 0 or end > len(seq) or end <= start:
        return ""
    return seq[start:end]


def _cox1_offsets(cfg: Config, species: str) -> int:
    from mtrnafeat.io.annotations import annotation_for
    annot = annotation_for(species, "COX1")
    return int(annot["l_utr5"])


# --- Directional substitution flux (was legacy/base_substitution/02.*) ---

_BASES = ("A", "C", "G", "T")


def _normalize(b: str) -> str:
    return "T" if b == "U" else b.upper()


def _benjamini_hochberg(pvals: list[float]) -> list[float]:
    n = len(pvals)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvals), key=lambda kv: kv[1])
    bh = [0.0] * n
    cur_min = 1.0
    for rank, (orig_idx, p) in enumerate(reversed(indexed)):
        k = n - rank
        adj = min(1.0, p * n / k)
        cur_min = min(cur_min, adj)
        bh[orig_idx] = cur_min
    return bh


def directional_flux_table(cfg: Config) -> pd.DataFrame:
    """For each codon position (1, 2, 3) and each ordered (from, to) pair of
    nucleotides (A↔T excluded as redundant), count yeast→human substitutions.

    Returns one row per (Position, From, To) with Count, expected (under
    matched-base background), binomial p-value, and Benjamini-Hochberg q.
    """
    try:
        from scipy.stats import binomtest as _binom

        def _binom_pval(k: int, n: int, p: float) -> float:
            return float(_binom(k, n=n, p=p, alternative="greater").pvalue)
    except ImportError:
        from scipy.stats import binom_test as _legacy  # type: ignore[attr-defined]

        def _binom_pval(k: int, n: int, p: float) -> float:
            return float(_legacy(k, n=n, p=p, alternative="greater"))

    table = alignment_table(cfg)
    rows: list[dict] = []
    for pos in (1, 2, 3):
        # Background: per-base frequency at this codon position in yeast.
        base_counts: dict[str, int] = {b: 0 for b in _BASES}
        for _, r in table.iterrows():
            if r["is_gap"]:
                continue
            base_counts[_normalize(r["yeast_codon"][pos - 1])] = (
                base_counts.get(_normalize(r["yeast_codon"][pos - 1]), 0) + 1
            )
        total = sum(base_counts.values())
        if total == 0:
            continue
        for from_b in _BASES:
            for to_b in _BASES:
                if from_b == to_b:
                    continue
                count = 0
                from_total = 0
                for _, r in table.iterrows():
                    if r["is_gap"] or r["n_diffs"] == 0:
                        continue
                    cy = r["yeast_codon"]; ch = r["human_codon"]
                    if str(pos) not in r["diff_positions"].split(","):
                        continue
                    yb = _normalize(cy[pos - 1])
                    hb = _normalize(ch[pos - 1])
                    if yb == from_b:
                        from_total += 1
                    if yb == from_b and hb == to_b:
                        count += 1
                expected_p = base_counts.get(to_b, 0) / total
                rows.append({
                    "Position": pos,
                    "From": from_b,
                    "To": to_b,
                    "Count": count,
                    "From_Total": from_total,
                    "Expected_p": expected_p,
                })
    flux = pd.DataFrame(rows)
    if flux.empty:
        return flux
    pvals = []
    for _, r in flux.iterrows():
        if r["From_Total"] > 0:
            try:
                p = _binom_pval(int(r["Count"]), n=int(r["From_Total"]),
                                 p=float(r["Expected_p"]))
            except Exception:
                p = 1.0
        else:
            p = 1.0
        pvals.append(float(p))
    flux["Binomial_p_one_sided"] = pvals
    flux["BH_q"] = _benjamini_hochberg(pvals)
    flux["Bias_Score"] = flux.apply(
        lambda r: (r["Count"] / r["From_Total"]) - r["Expected_p"]
                   if r["From_Total"] > 0 else 0.0,
        axis=1,
    )
    return flux


def transition_transversion_summary(table: pd.DataFrame) -> pd.DataFrame:
    """Ti/Tv ratio per codon position. Transitions: A↔G or C↔T."""
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    rows = []
    for pos in (1, 2, 3):
        ti = 0; tv = 0
        for _, r in table.iterrows():
            if r["is_gap"] or r["n_diffs"] == 0:
                continue
            if str(pos) not in r["diff_positions"].split(","):
                continue
            cy = _normalize(r["yeast_codon"][pos - 1])
            ch = _normalize(r["human_codon"][pos - 1])
            if (cy, ch) in transitions:
                ti += 1
            else:
                tv += 1
        rows.append({"Position": pos, "Transitions": ti, "Transversions": tv,
                     "Ti_over_Tv": ti / tv if tv > 0 else float("nan")})
    return pd.DataFrame(rows)
