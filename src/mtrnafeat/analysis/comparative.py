"""Yeast↔human COX1 codon-aligned comparative analysis.

Reads the PAL2NAL alignment, classifies every column as identical /
synonymous / non-synonymous, and locates substitutions by codon position
(1st / 2nd / wobble). The downstream tally feeds the substitution
heatmap and directional-flux outputs of the `compare` stage.
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.io.alignment import codon_position_changes, parse_pal2nal
from mtrnafeat.io.codons import HUMAN_MT_CODON_TABLE, YEAST_MT_CODON_TABLE


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
