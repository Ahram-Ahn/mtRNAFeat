"""Tests for io/ layer: parser round-trip, gene normalization, alignment."""
from __future__ import annotations

import pytest

from mtrnafeat.constants import canonical_gene
from mtrnafeat.io.alignment import codon_position_changes, parse_pal2nal
from mtrnafeat.io.annotations import annotation_for, annotation_df
from mtrnafeat.io.db_parser import get_record, list_genes, parse_db


def test_parse_db_records_match_lengths(mini_human_db, mini_yeast_db):
    for path in (mini_human_db, mini_yeast_db):
        records = parse_db(path)
        assert len(records) == 1
        rec = records[0]
        assert len(rec.sequence) == len(rec.structure)
        assert "T" not in rec.sequence  # T -> U normalized
        assert set(rec.structure) <= set(".()")


def test_db_parser_balanced_brackets(mini_human_db):
    rec = parse_db(mini_human_db)[0]
    depth = 0
    for ch in rec.structure:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            assert depth >= 0
    assert depth == 0


def test_canonical_gene_aliases():
    assert canonical_gene("ATP86") == "ATP8/ATP6"
    assert canonical_gene("ND4L4") == "ND4L/ND4"
    assert canonical_gene("COX1") == "COX1"


def test_get_record_by_alias(mini_human_db):
    rec = get_record(mini_human_db, "ND6")
    assert rec.gene == "ND6"


def test_list_genes_canonical(mini_yeast_db):
    genes = list_genes(mini_yeast_db)
    assert "ATP9" in genes


def test_annotation_for_handles_alias():
    annot = annotation_for("Human", "ATP86")
    assert annot["l_tr"] == 843
    assert annot["l_utr5"] == 1
    annot2 = annotation_for("Human", "ATP8/ATP6")
    assert annot == annot2


def test_annotation_df_round_trip():
    h = annotation_df("Human")
    y = annotation_df("Yeast")
    assert "transcript" in h.columns
    assert (h["l_tr"] == h["l_utr5"] + h["l_cds"] + h["l_utr3"]).all()
    assert (y["l_tr"] == y["l_utr5"] + y["l_cds"] + y["l_utr3"]).all()


def test_parse_pal2nal_columns_aligned(tmp_path):
    pal = tmp_path / "tiny.txt"
    pal.write_text(
        "PAL2NAL output\n\n"
        "                M   F   A\n"
        "S.cerevisiae    --- ATG GTA\n"
        "                M   L   A\n"
        "H.sapiens       ATG TTC GCC\n"
    )
    aln = parse_pal2nal(pal)
    assert aln.n_columns == 3
    assert aln.yeast_codons == ["---", "ATG", "GTA"]
    assert aln.human_codons == ["ATG", "TTC", "GCC"]
    assert aln.yeast_aa == ["-", "M", "V"]
    assert aln.human_aa == ["M", "F", "A"]


def test_codon_position_changes():
    assert codon_position_changes("ATG", "ACG") == [2]
    assert codon_position_changes("---", "ATG") == []
    assert codon_position_changes("AAA", "GGG") == [1, 2, 3]
