"""TIS zoom: −50 / +50 nt around the start codon (clamped at the 5' end).

For each (species, gene):
  - The downstream half is always +50 nt of CDS starting from the
    annotated TIS (or sequence end, whichever comes first).
  - The upstream half is min(L_5UTR, upstream_nt) nt of 5'UTR — i.e. as
    much real 5'UTR as exists, up to `upstream_nt` (default 50). Genes
    with no 5'UTR contribute 0 nt upstream and the window collapses to
    50 nt of CDS.

So genes with a full 5'UTR get the textbook −50/+50 = 100 nt window;
genes with a partial 5'UTR get (L_5UTR + 50) nt; genes with no 5'UTR get
50 nt. The CSV records `L_5UTR_in_window` and `L_CDS_in_window` so the
plot can stratify or annotate.
"""
from __future__ import annotations

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.projection import project_structure_to_window
from mtrnafeat.core.structure import paired_fraction
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import parse_db


def tis_table(cfg: Config, upstream_nt: int = 50, downstream_nt: int = 50) -> pd.DataFrame:
    rows = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        rec_by_gene = {r.gene: r for r in parse_db(path)}
        for gene in cfg.target_genes:
            gene_canon = canonical_gene(gene)
            if gene_canon not in rec_by_gene:
                continue
            rec = rec_by_gene[gene_canon]
            try:
                annot = annotation_for(species, gene)
            except KeyError:
                continue
            l_utr5 = int(annot.get("l_utr5", 0))
            tis_pos = l_utr5  # 0-based start of CDS
            up = min(l_utr5, upstream_nt)
            down = min(len(rec.sequence) - tis_pos, downstream_nt)
            start = tis_pos - up
            end = tis_pos + down
            if end - start < 20:
                continue
            seq_w = rec.sequence[start:end]
            dms_w = project_structure_to_window(rec.structure, start, end)
            try:
                dms_e = thermo.eval_structure(seq_w, dms_w)
            except Exception:
                dms_e = float("nan")
            try:
                vienna_struct, vienna_e = thermo.fold_mfe(seq_w)
            except Exception:
                vienna_struct, vienna_e = "", float("nan")
            rows.append({
                "Species": species,
                "Gene": gene_canon,
                "TIS_Pos_1based": tis_pos + 1,
                "Window_Start_1based": start + 1,
                "Window_End_1based": end,
                "Window_Length": end - start,
                "L_5UTR_total": l_utr5,
                "L_5UTR_in_window": up,
                "L_CDS_in_window": down,
                "Has_5UTR": bool(l_utr5 > 0),
                "Has_Full_5UTR_Context": bool(l_utr5 >= upstream_nt),
                "DMS_TIS_Energy": dms_e,
                "DMS_TIS_PairedFrac": paired_fraction(dms_w),
                "Vienna_TIS_MFE": vienna_e,
                "Vienna_TIS_PairedFrac": paired_fraction(vienna_struct),
                "Delta_E_DMS_minus_MFE": dms_e - vienna_e if (pd.notna(dms_e) and pd.notna(vienna_e)) else float("nan"),
            })
    return pd.DataFrame(rows)
