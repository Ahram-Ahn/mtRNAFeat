"""Per-position local pair-probability scan via ViennaRNA's RNAplfold.

For each transcript:

* ``p_paired[i]`` — Boltzmann probability that position ``i`` is paired
  with *some* partner within ``max_bp_span`` nt, averaged across all
  overlapping windows of length ``window``.
* ``p_unpaired_u1[i] = 1 − p_paired[i]`` — probability that position ``i``
  is single-stranded.

DMS reactivity is **not** fed in — this is the purely-thermodynamic local
prior that complements the .db DMS-guided structure. Re-fitting DMS into
RNAplfold would re-introduce the very pseudo-energy fitting we're trying
to avoid; that path is RNAstructure's ``Fold -dms`` (already used by the
``window`` command).

Multi-length accessibility (probability a length-u window is fully unpaired
for u ∈ {5, 10, 20}) is intentionally deferred. ``RNA.pfl_fold`` gives us
the pair-prob matrix in a single fast call; the multi-u variant requires
the callback-based ``probs_window`` API and is worth a separate slice.
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
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import parse_db
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


def scan_one_gene(species: str, gene: str, sequence: str,
                   window: int, max_bp_span: int, cutoff: float) -> LocalProbResult:
    """Run RNAplfold-style folding on one transcript and return p_paired track."""
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
    return LocalProbResult(
        species=species,
        gene=canonical_gene(gene),
        sequence=sequence,
        p_paired=p_paired,
        window=int(window),
        max_bp_span=int(max_bp_span),
        cutoff=float(cutoff),
    )


def per_position_table(result: LocalProbResult, annot: dict | None = None) -> pd.DataFrame:
    """Per-nt DataFrame with paired / unpaired probabilities and region tag."""
    n = len(result.sequence)
    if annot is not None:
        cds_start = annot["l_utr5"] + 1
        cds_end = annot["l_utr5"] + annot["l_cds"]
        region: list[str] = []
        for pos in range(1, n + 1):
            if pos < cds_start:
                region.append("5'UTR")
            elif pos <= cds_end:
                region.append("CDS")
            else:
                region.append("3'UTR")
    else:
        region = [""] * n
    return pd.DataFrame({
        "Species": result.species,
        "Gene": result.gene,
        "Position_1based": np.arange(1, n + 1),
        "Nucleotide": list(result.sequence),
        "Region": region,
        "P_Paired": result.p_paired,
        "P_Unpaired_u1": 1.0 - result.p_paired,
        "Window_nt": result.window,
        "Max_BP_Span_nt": result.max_bp_span,
        "Cutoff": result.cutoff,
    })


def scan_all(cfg: Config, window: int, max_bp_span: int,
             cutoff: float) -> tuple[pd.DataFrame, list[LocalProbResult]]:
    """Walk every (species, gene) in cfg.target_genes and emit a long-form
    per-position DataFrame plus the raw per-gene results (for plotting)."""
    targets = {canonical_gene(g) for g in (cfg.target_genes or ())}
    step(f"local-probability scan (W={window}, L={max_bp_span}, cutoff={cutoff})")
    results: list[LocalProbResult] = []
    for species, fname in cfg.db_files.items():
        recs = [r for r in parse_db(cfg.data_dir / fname)
                if (not targets or canonical_gene(r.gene) in targets)]
        for rec in progress(recs, desc=f"{species} plfold", unit="gene"):
            res = scan_one_gene(species, rec.gene, rec.sequence,
                                window=window, max_bp_span=max_bp_span,
                                cutoff=cutoff)
            results.append(res)
    if not results:
        return pd.DataFrame(), []
    frames = []
    for res in results:
        try:
            annot = annotation_for(res.species, res.gene)
        except KeyError:
            annot = None
        frames.append(per_position_table(res, annot))
    return pd.concat(frames, ignore_index=True), results
