"""Sliding-window local-folding scan with a single max_bp_span constraint.

Replaces legacy 11. Each window is folded by plain Vienna MFE (with
max_bp_span). The DMS-projected pair-table is energy-evaluated under the
same model; pairs whose distance exceeds max_bp_span are stripped from
the projection so DMS and Vienna are compared apples-to-apples.

CoFold is no longer scored here; the dedicated `cofold` stage explores
its parameter space against DMS in a separate analysis.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.projection import project_structure_to_window
from mtrnafeat.core.structure import filter_max_bp_span, paired_fraction
from mtrnafeat.io.annotations import annotation_for, classify_region
from mtrnafeat.io.db_parser import get_record
from mtrnafeat.progress import progress, step


def sliding_intervals(n: int, window: int, step: int, include_terminal: bool = True) -> list[tuple[int, int]]:
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
class WindowScanResult:
    species: str
    gene: str
    sequence: str
    dms_structure: str
    annot: dict
    window_nt: int
    step_nt: int
    max_bp_span: int | None
    windows: pd.DataFrame


def scan_one_gene(species: str, db_path: str | Path, gene: str, cfg: Config,
                  max_bp_span: int | None = None) -> WindowScanResult:
    """Slide windows over the gene; per-window: DMS-eval'd energy + Vienna MFE.

    Both use max_bp_span (defaulting to cfg.max_bp_span) as the hard cutoff —
    this gives apples-to-apples energies because the same long-range pairs
    are excluded from the DMS projection and from the Vienna optimum.
    """
    if max_bp_span is None:
        max_bp_span = cfg.max_bp_span
    rec = get_record(db_path, gene)
    annot = annotation_for(species, gene)
    intervals = sliding_intervals(len(rec.sequence), cfg.window_nt, cfg.step_nt)
    rows = []
    bar = progress(list(enumerate(intervals, start=1)),
                    desc=f"{species} {canonical_gene(gene)} window",
                    unit="win", leave=False)
    for idx, (s, e) in bar:
        seq_w = rec.sequence[s:e]
        dms_proj = project_structure_to_window(rec.structure, s, e)
        dms_proj_span, removed = filter_max_bp_span(dms_proj, max_bp_span)
        e_dms = thermo.eval_structure(seq_w, dms_proj_span, max_bp_span=max_bp_span)
        v_struct, v_mfe = thermo.fold_mfe(seq_w, max_bp_span=max_bp_span)
        center = s + 1 + (len(seq_w) - 1) / 2.0
        rows.append({
            "Species": species,
            "Gene": canonical_gene(gene),
            "Window_Index": idx,
            "Window_Start_1based": s + 1,
            "Window_End_1based": e,
            "Window_Length_nt": len(seq_w),
            "Window_Center_1based": center,
            "Region_At_Center": classify_region(annot["l_utr5"], annot["l_cds"], int(center)),
            "DMS_Window_Energy": e_dms,
            "DMS_Window_PairedFraction": paired_fraction(dms_proj_span),
            "DMS_Pairs_Removed_For_Span": removed,
            "Vienna_Window_MFE": v_mfe,
            "Vienna_Window_PairedFraction": paired_fraction(v_struct),
            "Max_BP_Span": int(max_bp_span),
        })
    df = pd.DataFrame(rows)
    return WindowScanResult(
        species=species, gene=canonical_gene(gene),
        sequence=rec.sequence, dms_structure=rec.structure, annot=annot,
        window_nt=cfg.window_nt, step_nt=cfg.step_nt, max_bp_span=max_bp_span,
        windows=df,
    )


def scan_for_gene(species: str, db_path, gene: str, cfg: Config) -> pd.DataFrame:
    """Single-span window scan for one gene."""
    step(f"window scan: {species} {canonical_gene(gene)}")
    res = scan_one_gene(species, db_path, gene, cfg, cfg.max_bp_span)
    return res.windows


def summarize_windows(df: pd.DataFrame) -> pd.DataFrame:
    """Per (species, gene) summary: mean/min/max of energy and paired fraction."""
    grp = df.groupby(["Species", "Gene"])
    return grp.agg(
        N_Windows=("Window_Index", "count"),
        DMS_Energy_Mean=("DMS_Window_Energy", "mean"),
        DMS_Energy_Min=("DMS_Window_Energy", "min"),
        DMS_Energy_Max=("DMS_Window_Energy", "max"),
        Vienna_MFE_Mean=("Vienna_Window_MFE", "mean"),
        Vienna_MFE_Min=("Vienna_Window_MFE", "min"),
        Vienna_MFE_Max=("Vienna_Window_MFE", "max"),
        DMS_PairedFrac_Mean=("DMS_Window_PairedFraction", "mean"),
        Vienna_PairedFrac_Mean=("Vienna_Window_PairedFraction", "mean"),
    ).reset_index()
