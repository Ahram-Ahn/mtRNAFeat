"""Per-transcript statistics: length, MFE, foldedness, sequence GC,
pair-typed composition (G-C / A-U / G-U), centroid + MEA bp-distance.

Replaces legacy/03.1.generate_statistics.py and the broken pair counting in
legacy/03_integrate_sim_dms.py (bug #7).
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.core.structure import paired_fraction
from mtrnafeat.io.db_parser import DbRecord, parse_db


def sequence_gc_pct(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    return 100.0 * (s.count("G") + s.count("C")) / len(s)


def paired_composition(seq: str, structure: str) -> tuple[float, float, float]:
    """Return (paired_GC%, paired_AU%, paired_GU%) — pair-typed.

    Bug fix: this is the pair-typed counterpart of legacy 03.1's
    extract_pair_composition. The other legacy module
    (03_integrate_sim_dms.calc_paired_composition) counted bases instead of
    pairs and put G-U wobble's G into "GC" and U into "AT", which is wrong.
    """
    stack: list[int] = []
    gc = au = gu = 0
    total = 0
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                continue
            j = stack.pop()
            pair = {seq[i].upper(), seq[j].upper()}
            if pair == {"G", "C"}:
                gc += 1
            elif pair == {"A", "U"}:
                au += 1
            elif pair == {"G", "U"}:
                gu += 1
            total += 1
    if total == 0:
        return 0.0, 0.0, 0.0
    return 100.0 * gc / total, 100.0 * au / total, 100.0 * gu / total


def transcript_stats(rec: DbRecord, condition: str) -> dict:
    length = len(rec.sequence)
    fold_pct = 100.0 * paired_fraction(rec.structure)
    gc_pair, au_pair, gu_pair = paired_composition(rec.sequence, rec.structure)
    return {
        "Condition": condition,
        "Gene": rec.gene,
        "Length": length,
        "MFE": rec.mfe,
        "Normalized_MFE_per_nt": rec.mfe / length if length else 0.0,
        "Foldedness_Pct": fold_pct,
        "Sequence_GC_Pct": sequence_gc_pct(rec.sequence),
        "Paired_GC_Pct": gc_pair,
        "Paired_AU_Pct": au_pair,
        "Paired_GU_Pct": gu_pair,
    }


def per_transcript_stats(db_path, condition: str) -> pd.DataFrame:
    """Run transcript_stats over every record in a .db file."""
    rows = [transcript_stats(rec, condition) for rec in parse_db(db_path)]
    return pd.DataFrame(rows)


def add_centroid_distances(df: pd.DataFrame, db_path) -> pd.DataFrame:
    """Augment a stats DataFrame with bp_distance(DMS, MFE) and
    bp_distance(DMS, centroid) and bp_distance(DMS, MEA)."""
    from mtrnafeat.core import thermo
    rec_by_gene = {r.gene: r for r in parse_db(db_path)}
    out = df.copy()
    bp_mfe, bp_centroid, bp_mea = [], [], []
    for gene in out["Gene"]:
        rec = rec_by_gene[gene]
        try:
            mfe_struct, _ = thermo.fold_mfe(rec.sequence)
            cen_struct, _ = thermo.centroid_structure(rec.sequence)
            mea_struct, _ = thermo.mea_structure(rec.sequence)
            bp_mfe.append(thermo.bp_distance(rec.structure, mfe_struct))
            bp_centroid.append(thermo.bp_distance(rec.structure, cen_struct))
            bp_mea.append(thermo.bp_distance(rec.structure, mea_struct))
        except Exception:
            bp_mfe.append(float("nan"))
            bp_centroid.append(float("nan"))
            bp_mea.append(float("nan"))
    out["BPdist_DMS_vs_MFE"] = bp_mfe
    out["BPdist_DMS_vs_Centroid"] = bp_centroid
    out["BPdist_DMS_vs_MEA"] = bp_mea
    return out
