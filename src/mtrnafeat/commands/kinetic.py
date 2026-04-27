"""`mtrnafeat kinetic` — DrTransformer kinetic-folding for selected genes."""
from __future__ import annotations

from mtrnafeat.analysis import kinetic
from mtrnafeat.config import Config
from mtrnafeat.constants import file_safe_gene
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import kinetic_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "kinetic"
    out.mkdir(parents=True, exist_ok=True)
    if not kinetic.has_drtransformer():
        print("[mtrnafeat] DrTransformer not on PATH. Install with `pip install drtransformer`.")
        return 2
    genes = ["COX1", "ND6", "ATP9"]
    if args:
        for i, a in enumerate(args):
            if a == "--genes" and i + 1 < len(args):
                genes = args[i + 1].split(",")
                break
    summary, traj = kinetic.run_kinetic_for_genes(cfg, genes)
    canonical_csv(summary, out / "kinetic_summary.csv")
    canonical_csv(traj, out / "kinetic_trajectory.csv")
    for (gene, species), _ in traj.groupby(["Gene", "Species"]):
        kinetic_plot.plot_kinetic_trajectory(traj, gene, species,
                                              plot_path(out, f"kinetic_{species}_{file_safe_gene(gene)}", cfg.plot_format), dpi=cfg.dpi)
    return 0
