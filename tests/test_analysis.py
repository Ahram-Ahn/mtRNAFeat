"""Tests that exercise the analysis layer without ViennaRNA."""
from __future__ import annotations

import pytest

from mtrnafeat.analysis import features, landscape, tis
from mtrnafeat.analysis.statistics import paired_composition, sequence_gc_pct, transcript_stats
from mtrnafeat.config import Config
from mtrnafeat.io.db_parser import DbRecord, parse_db


def test_sequence_gc_pct():
    assert sequence_gc_pct("AAGCAA") == pytest.approx(2 / 6 * 100)
    assert sequence_gc_pct("") == 0.0


def test_paired_composition_pair_typed():
    seq2 = "GCAUUUUGC"
    struct2 = "(((...)))"
    gc, au, gu = paired_composition(seq2, struct2)
    # Pairs from outside-in: (G, C) -> GC, (C, G) -> GC, (A, U) -> AU
    assert gc == pytest.approx(200.0 / 3)  # 2/3
    assert au == pytest.approx(100.0 / 3)
    assert gu == 0.0


def test_transcript_stats_columns(mini_human_db):
    rec = parse_db(mini_human_db)[0]
    row = transcript_stats(rec, condition="X")
    expected = {"Condition", "Gene", "Length", "MFE", "Normalized_MFE_per_nt",
                "Foldedness_Pct", "Sequence_GC_Pct", "Paired_GC_Pct",
                "Paired_AU_Pct", "Paired_GU_Pct"}
    assert expected.issubset(row.keys())
    assert row["Length"] == 538
    assert row["MFE"] == -26.3


def test_features_simulated_uses_configured_species_nulls(monkeypatch):
    seen_sequences = []

    def fake_fold(seq):
        seen_sequences.append(seq)
        return "(" * (len(seq) // 2) + "." * (len(seq) % 2) + ")" * (len(seq) // 2), -1.0

    monkeypatch.setattr(features.thermo, "fold_mfe", fake_fold)
    cfg = Config(
        db_files={"Test": "unused.db"},
        sim_freqs_per_species={"Test": {"A": 1.0, "U": 0.0, "G": 0.0, "C": 0.0}},
        sim_num_sequences=3,
        sim_seq_length=7,
    )

    motifs, spans = features.features_simulated(cfg)

    assert len(seen_sequences) == 3
    assert set(seen_sequences) == {"AAAAAAA"}
    assert set(motifs["Species"]) == {"Test"}
    assert set(spans["Species"]) == {"Test"}


def test_landscape_experimental_overlay_recalculates_db_energy(monkeypatch, tmp_path):
    db = tmp_path / "sample.db"
    db.write_text(">G1: -99.0 kcal/mol\nGCAU\n(())\n")
    monkeypatch.setattr(landscape.thermo, "eval_structure", lambda seq, struct: -4.5)
    cfg = Config(data_dir=tmp_path, db_files={"Sample": "sample.db"})

    df = landscape.experimental_overlay(cfg)

    row = df.iloc[0]
    assert row["MFE"] == -4.5
    assert row["Header_MFE"] == -99.0
    assert "eval_structure" in row["MFE_Source"]


def test_region_overlay_uses_sample_annotation_mapping(monkeypatch, tmp_path):
    db = tmp_path / "ko.db"
    db.write_text(">COX2: -1.0 kcal/mol\n" + "A" * 30 + "\n" + "." * 30 + "\n")
    monkeypatch.setattr(landscape.thermo, "eval_structure", lambda seq, struct: -float(len(seq)))
    monkeypatch.setattr(landscape, "annotation_for_sample", lambda sample, gene, mapping: {
        "l_tr": 30,
        "l_utr5": 5,
        "l_cds": 20,
        "l_utr3": 5,
    })
    monkeypatch.setattr(landscape, "infer_annotation_species", lambda sample, mapping: mapping.get(sample))
    cfg = Config(
        data_dir=tmp_path,
        db_files={"KO1": "ko.db"},
        sample_annotation_species={"KO1": "Yeast"},
    )

    df = landscape.experimental_overlay_regions(cfg)

    assert set(df["Region"]) == {"5'UTR", "CDS", "3'UTR/tail"}
    assert set(df["Annotation_Species"]) == {"Yeast"}
    assert df.set_index("Region").loc["CDS", "MFE"] == -20.0


def test_tis_downstream_window_stays_inside_cds(monkeypatch, tmp_path):
    rec = DbRecord(
        gene="COX1",
        raw_gene="COX1",
        mfe=-1.0,
        sequence="A" * 40,
        structure="." * 40,
    )

    monkeypatch.setattr(tis, "parse_db", lambda _path: [rec])
    monkeypatch.setattr(tis, "annotation_for_sample", lambda _species, _gene, _mapping: {
        "l_tr": 40,
        "l_utr5": 18,
        "l_cds": 5,
        "l_utr3": 17,
    })
    monkeypatch.setattr(tis.thermo, "eval_structure", lambda seq, struct: -float(len(seq)))
    monkeypatch.setattr(tis.thermo, "fold_mfe", lambda seq: ("." * len(seq), -float(len(seq))))
    cfg = Config(
        data_dir=tmp_path,
        db_files={"Test": "unused.db"},
        target_genes=("COX1",),
    )

    df = tis.tis_table(cfg, upstream_nt=18, downstream_nt=50)

    row = df.iloc[0]
    assert row["Window_Start_1based"] == 1
    assert row["Window_End_1based"] == 23
    assert row["Window_Length"] == 23
    assert row["L_5UTR_in_window"] == 18
    assert row["L_CDS_in_window"] == 5
