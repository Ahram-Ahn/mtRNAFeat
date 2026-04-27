"""`mtrnafeat gene-panel` — composition + paired-pair + local-foldedness per gene.

Generalizes the legacy ND6-only panel to every (species, gene) the .db
files contain. Outputs one figure per gene under `gene_panels/`.
"""
from __future__ import annotations

from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene, file_safe_gene
from mtrnafeat.io.annotations import annotation_for
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.viz import gene_panel
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "gene_panels"
    out.mkdir(parents=True, exist_ok=True)
    target_set = {canonical_gene(g) for g in cfg.target_genes}
    plotted = 0
    for species, fname in cfg.db_files.items():
        for rec in parse_db(cfg.data_dir / fname):
            target = canonical_gene(rec.gene)
            if target_set and target not in target_set:
                continue
            try:
                annot = annotation_for(species, rec.gene)
            except KeyError:
                annot = None
            gene_panel.plot_gene(
                rec,
                species=species,
                annot=annot,
                out_path=plot_path(out, f"panel_{species}_{file_safe_gene(rec.gene)}", cfg.plot_format),
                dpi=cfg.dpi,
            )
            plotted += 1
    print(f"[mtrnafeat] gene-panel: wrote {plotted} panel(s) to {out}")
    return 0
