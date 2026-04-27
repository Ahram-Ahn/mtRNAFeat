"""5'-end zoom: first `window_len` nt of every transcript, DMS vs Vienna MFE.

The previous −50/+50 framing degraded to "50 nt of CDS" for the 8/11 human
genes whose annotated 5'UTR is 0 nt. This version takes the unconditional
first `window_len` nt of the transcript (default 50). Genes with a real
5'UTR get 5'UTR + (window_len − 5'UTR) of CDS; genes with no 5'UTR get
pure CDS at the start codon. The summary CSV preserves the biological
context via `Has_5UTR` and `L_5UTR` columns so the user can stratify the
plot afterwards.
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


def tis_table(cfg: Config, window_len: int = 50) -> pd.DataFrame:
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
            start = 0
            end = min(len(rec.sequence), window_len)
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
            l_utr5 = int(annot.get("l_utr5", 0))
            rows.append({
                "Species": species,
                "Gene": gene_canon,
                "Window_Start_1based": start + 1,
                "Window_End_1based": end,
                "Window_Length": end - start,
                "Has_5UTR": bool(l_utr5 > 0),
                "L_5UTR": l_utr5,
                "DMS_5end_Energy": dms_e,
                "DMS_5end_PairedFrac": paired_fraction(dms_w),
                "Vienna_5end_MFE": vienna_e,
                "Vienna_5end_PairedFrac": paired_fraction(vienna_struct),
                "Delta_E_DMS_minus_MFE": dms_e - vienna_e if (pd.notna(dms_e) and pd.notna(vienna_e)) else float("nan"),
            })
    return pd.DataFrame(rows)
