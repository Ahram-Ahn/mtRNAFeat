"""Whole-transcript fold-and-compare for the `mtrnafeat window` command.

For each transcript we compute three position-resolved paired-fraction tracks:

* **DMS** — the dot-bracket from the .db file (DMS-MaPseq + RNAstructure
  upstream pipeline). Energy is *recalculated* under ViennaRNA so the
  legend ΔG is on the same model as the engine fold below.
* **Vienna full** — ``RNA.fold`` over the entire transcript with no span
  cap. The pure-thermodynamic baseline, engine-agnostic.
* **Engine span** — whole-transcript fold with ``max_bp_span = cfg.max_bp_span``
  (default 300) under the engine selected by ``cfg.fold_engine``:
    - ``"vienna"``: ``RNA.fold_compound`` with ``md.max_bp_span``.
    - ``"rnastructure"``: ``Fold -md <max_bp_span>``.

The per-position 0/1 paired vector is smoothed by a centered rolling mean
of width ``cfg.rolling_window`` (default 25 nt) before plotting; the raw
0/1 vectors are also written out so downstream code can re-smooth.

``sliding_intervals`` is kept here for backwards compatibility — it's
imported by ``analysis.significance`` and ``analysis.cofold_sweep``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.structure import filter_max_bp_span
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import get_record
from mtrnafeat.progress import step

_VALID_ENGINES = ("rnastructure", "vienna")


def sliding_intervals(n: int, window: int, step: int, include_terminal: bool = True) -> list[tuple[int, int]]:
    """Yield ``(start, end)`` half-open intervals for a sliding-window scan.

    Kept here for ``analysis.significance`` and ``analysis.cofold_sweep`` —
    the new ``window`` command does whole-transcript folds and no longer
    calls this directly.
    """
    if window <= 0 or step <= 0:
        raise ValueError("window and step must be positive")
    if window > n:
        return [(0, n)]
    starts = list(range(0, n - window + 1, step))
    intervals = [(s, s + window) for s in starts]
    if include_terminal and intervals[-1][1] < n:
        intervals.append((n - window, n))
    seen = set()
    out: list[tuple[int, int]] = []
    for s, e in intervals:
        if (s, e) not in seen:
            out.append((s, e))
            seen.add((s, e))
    return out


@dataclass
class TranscriptFoldResult:
    species: str
    gene: str
    sequence: str
    annot: dict
    max_bp_span: int
    engine: str
    # Structures (dot-brackets, all length == len(sequence))
    dms_structure: str
    dms_span_structure: str
    dms_pairs_removed: int
    vienna_full_structure: str
    engine_span_structure: str
    # Energies (kcal/mol, all evaluated under Vienna for apples-to-apples)
    dms_recalc_mfe: float
    dms_span_mfe: float
    vienna_full_mfe: float
    engine_span_mfe_native: float    # the engine's own MFE (Mathews '04 for RNAstructure)
    engine_span_mfe_vienna_eval: float  # same structure, scored by Vienna for direct comparison


def _engine_fold(seq: str, max_bp_span: int, engine: str) -> tuple[str, float]:
    """Whole-transcript fold under the requested engine. Returns (structure, native_MFE)."""
    if engine == "vienna":
        return thermo.fold_mfe(seq, max_bp_span=max_bp_span)
    if engine == "rnastructure":
        from mtrnafeat.engines import rnastructure
        res = rnastructure.fold(seq, max_bp_span=max_bp_span)
        return res.dot_bracket, float(res.mfe)
    raise ValueError(
        f"Unknown fold_engine: {engine!r} (expected one of {_VALID_ENGINES})"
    )


def fold_transcript(species: str, db_path: str | Path, gene: str, cfg: Config,
                    max_bp_span: int | None = None,
                    engine: str | None = None) -> TranscriptFoldResult:
    """Whole-transcript four-way fold compare for one gene.

    Reads sequence + DMS dot-bracket from the .db file, then folds with the
    requested engine at full-transcript scale. Returns dot-brackets, ΔGs,
    and the count of long-range DMS pairs removed by the span filter.
    """
    if max_bp_span is None:
        max_bp_span = cfg.max_bp_span
    if engine is None:
        engine = cfg.fold_engine
    if engine not in _VALID_ENGINES:
        raise ValueError(
            f"Unknown fold_engine: {engine!r} (expected one of {_VALID_ENGINES})"
        )
    rec = get_record(db_path, gene)
    annot = annotation_for(species, gene)

    # DMS structure from the .db, with span-sanitized variant.
    dms_struct = rec.structure
    dms_span_struct, removed = filter_max_bp_span(dms_struct, max_bp_span)
    dms_recalc_mfe = thermo.eval_structure(rec.sequence, dms_struct)
    dms_span_mfe = thermo.eval_structure(rec.sequence, dms_span_struct, max_bp_span=max_bp_span)

    # Vienna-full baseline (no span limit).
    vfull_struct, vfull_mfe = thermo.fold_mfe(rec.sequence)

    # Engine-span fold (the headline track).
    engine_span_struct, engine_span_mfe_native = _engine_fold(rec.sequence, max_bp_span, engine)
    # Score the same structure under Vienna for an apples-to-apples ΔG vs the DMS recalc.
    engine_span_mfe_vienna = thermo.eval_structure(
        rec.sequence, engine_span_struct, max_bp_span=max_bp_span
    )

    return TranscriptFoldResult(
        species=species,
        gene=canonical_gene(gene),
        sequence=rec.sequence,
        annot=annot,
        max_bp_span=int(max_bp_span),
        engine=engine,
        dms_structure=dms_struct,
        dms_span_structure=dms_span_struct,
        dms_pairs_removed=removed,
        vienna_full_structure=vfull_struct,
        engine_span_structure=engine_span_struct,
        dms_recalc_mfe=dms_recalc_mfe,
        dms_span_mfe=dms_span_mfe,
        vienna_full_mfe=vfull_mfe,
        engine_span_mfe_native=engine_span_mfe_native,
        engine_span_mfe_vienna_eval=engine_span_mfe_vienna,
    )


def _paired_vector(structure: str) -> list[int]:
    return [1 if ch in "()" else 0 for ch in structure]


def _rolling_mean(values: list[int], window: int) -> list[float]:
    """Centered rolling mean with shrinking window at the edges."""
    if window <= 1 or not values:
        return [float(v) for v in values]
    s = pd.Series(values, dtype=float)
    return s.rolling(window=int(window), center=True, min_periods=1).mean().tolist()


def per_position_table(res: TranscriptFoldResult, rolling_window: int) -> pd.DataFrame:
    """Per-nucleotide table with raw 0/1 paired vectors and their rolling means."""
    n = len(res.sequence)
    cds_start_1 = res.annot["l_utr5"] + 1
    cds_end_1 = res.annot["l_utr5"] + res.annot["l_cds"]
    region: list[str] = []
    cds_rel: list[int | None] = []
    for pos in range(1, n + 1):
        if pos < cds_start_1:
            region.append("5'UTR")
            cds_rel.append(None)
        elif pos <= cds_end_1:
            region.append("CDS")
            cds_rel.append(pos - cds_start_1 + 1)
        else:
            region.append("3'UTR")
            cds_rel.append(None)

    dms = _paired_vector(res.dms_structure)
    dms_span = _paired_vector(res.dms_span_structure)
    vfull = _paired_vector(res.vienna_full_structure)
    espan = _paired_vector(res.engine_span_structure)
    span = res.max_bp_span
    eng = res.engine.capitalize() if res.engine == "vienna" else "RNAstructure"
    return pd.DataFrame({
        "Species": res.species,
        "Gene": res.gene,
        "Position_1based": list(range(1, n + 1)),
        "Nucleotide": list(res.sequence),
        "Region": region,
        "CDS_Position_1based": cds_rel,
        "DMS_Paired": dms,
        f"DMS_Span{span}_Paired": dms_span,
        "ViennaFull_Paired": vfull,
        f"{eng}Span{span}_Paired": espan,
        "DMS_RollingPairedFrac": _rolling_mean(dms, rolling_window),
        f"DMS_Span{span}_RollingPairedFrac": _rolling_mean(dms_span, rolling_window),
        "ViennaFull_RollingPairedFrac": _rolling_mean(vfull, rolling_window),
        f"{eng}Span{span}_RollingPairedFrac": _rolling_mean(espan, rolling_window),
        "Rolling_Window_nt": int(rolling_window),
        "Fold_Engine": res.engine,
    })


def summarize_transcript(res: TranscriptFoldResult) -> pd.DataFrame:
    """One-row summary per (species, gene) — lengths and ΔGs for every track."""
    n = len(res.sequence)
    cds_end = min(n, res.annot["l_utr5"] + res.annot["l_cds"])
    span = res.max_bp_span
    eng = res.engine.capitalize() if res.engine == "vienna" else "RNAstructure"
    return pd.DataFrame([{
        "Species": res.species,
        "Gene": res.gene,
        "Fold_Engine": res.engine,
        "Max_BP_Span": span,
        "Transcript_Length_Observed": n,
        "Transcript_Length_Annotated": int(res.annot["l_tr"]),
        "UTR5_Length": int(res.annot["l_utr5"]),
        "CDS_Length": int(res.annot["l_cds"]),
        "CDS_Start_1based": res.annot["l_utr5"] + 1,
        "CDS_End_1based": cds_end,
        "DMS_Recalc_MFE": res.dms_recalc_mfe,
        f"DMS_Span{span}_MFE": res.dms_span_mfe,
        f"DMS_Pairs_Removed_For_Span{span}": res.dms_pairs_removed,
        "Vienna_Full_MFE": res.vienna_full_mfe,
        f"{eng}_Span{span}_MFE": res.engine_span_mfe_native,
        f"{eng}_Span{span}_MFE_ViennaEval": res.engine_span_mfe_vienna_eval,
    }])


def step_log(species: str, gene: str) -> None:
    """Lightweight progress hook so commands/window.py can log via the package style."""
    step(f"transcript fold: {species} {canonical_gene(gene)}")
