"""Structural-element decomposition + region stratification.

Replaces legacy 10. Adds 5'UTR / CDS / 3'UTR stratification (the legacy
script ignored region annotations).
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.shuffle import random_gc_sequence
from mtrnafeat.core.structure import extract_pairs, parse_element_sizes
from mtrnafeat.io.annotations import annotation_for, classify_region
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step
from mtrnafeat.rng import make_rng


def extract_motifs_from_record(species: str, source: str, gene: str, structure: str,
                                max_artifact: int) -> list[dict]:
    elements = parse_element_sizes(structure, max_artifact)
    rows = []
    for motif, sizes in (
        ("Macro_Helix", elements.macro_helix),
        ("Hairpin", elements.hairpin),
        ("Bulge", elements.bulge),
        ("Internal_Loop", elements.internal_loop),
        ("Multi_Branch", elements.multi_branch),
    ):
        for size in sizes:
            if size > 0:
                rows.append({"Species": species, "Type": source, "Gene": gene, "Motif": motif, "Size": int(size)})
    return rows


def extract_spans_from_record(species: str, source: str, gene: str, structure: str) -> list[dict]:
    return [
        {"Species": species, "Type": source, "Gene": gene, "Span": int(j - i)}
        for i, j in extract_pairs(structure)
    ]


def features_dms(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    motifs, spans = [], []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        for rec in parse_db(path):
            motifs.extend(extract_motifs_from_record(
                species, "DMS", rec.gene, rec.structure, cfg.max_loop_artifact_size
            ))
            spans.extend(extract_spans_from_record(species, "DMS", rec.gene, rec.structure))
    return pd.DataFrame(motifs), pd.DataFrame(spans)


def features_simulated(cfg: Config, conditions: list[tuple[str, float]] | None = None,
                        n_per_condition: int = 500, length: int = 300) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Simulate per-species GC-matched controls. Species column is the real species
    name (Human/Yeast); the Type column ("Sim") is the disambiguator from DMS.
    """
    if conditions is None:
        conditions = [("Human", 0.46), ("Yeast", 0.20)]
    rng = make_rng(cfg.seed + 2)
    motifs, spans = [], []
    step(f"simulating features for {len(conditions)} conditions × {n_per_condition} sequences")
    for species, gc in conditions:
        for k in progress(range(n_per_condition), desc=f"sim {species}", unit="seq"):
            seq = random_gc_sequence(length, gc, rng)
            struct, _ = thermo.fold_mfe(seq)
            sim_id = f"Sim_{species}_{k}"
            motifs.extend(extract_motifs_from_record(species, "Sim", sim_id, struct, cfg.max_loop_artifact_size))
            spans.extend(extract_spans_from_record(species, "Sim", sim_id, struct))
    return pd.DataFrame(motifs), pd.DataFrame(spans)


def region_stratified_elements(cfg: Config) -> pd.DataFrame:
    """For each DMS structure, classify each base-pair span by the region of
    its left endpoint and report its (motif, size). Uses io.annotations."""
    rows = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        for rec in parse_db(path):
            try:
                annot = annotation_for(species, rec.gene)
            except KeyError:
                continue
            elements = parse_element_sizes(rec.structure, cfg.max_loop_artifact_size)
            pairs = extract_pairs(rec.structure)
            # Map every pair to a region by the left index (1-based).
            for i, j in pairs:
                region = classify_region(annot["l_utr5"], annot["l_cds"], i + 1)
                rows.append({
                    "Species": species,
                    "Gene": canonical_gene(rec.gene),
                    "Region": region,
                    "Pair_Span": int(j - i),
                    "Pair_Left_1based": int(i + 1),
                    "Pair_Right_1based": int(j + 1),
                })
            for motif, sizes in (
                ("Macro_Helix", elements.macro_helix),
                ("Hairpin", elements.hairpin),
                ("Bulge", elements.bulge),
                ("Internal_Loop", elements.internal_loop),
                ("Multi_Branch", elements.multi_branch),
            ):
                for size in sizes:
                    rows.append({
                        "Species": species,
                        "Gene": canonical_gene(rec.gene),
                        "Motif": motif,
                        "Size": int(size),
                    })
    return pd.DataFrame(rows)
