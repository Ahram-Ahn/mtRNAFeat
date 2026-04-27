"""`mtrnafeat compare` — yeast↔human COX1 codon-aligned comparative analysis.

Skips gracefully if the prerequisites are missing (alignment file or COX1
record in either species), so the pipeline keeps moving on test fixtures
that don't include COX1.
"""
from __future__ import annotations

from mtrnafeat.analysis import comparative
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv, tables_csv
from mtrnafeat.viz import comparative_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "compare"
    out.mkdir(parents=True, exist_ok=True)

    alignment_path = cfg.data_dir / cfg.alignment_file
    if not alignment_path.exists():
        print(f"[mtrnafeat] compare: skipping — alignment file not found at {alignment_path}")
        return 0
    try:
        from mtrnafeat.io.db_parser import get_record
        get_record(cfg.data_dir / cfg.db_files.get("Yeast", ""), "COX1")
        get_record(cfg.data_dir / cfg.db_files.get("Human", ""), "COX1")
    except (KeyError, FileNotFoundError) as e:
        print(f"[mtrnafeat] compare: skipping — COX1 prerequisite missing ({e})")
        return 0

    table = comparative.alignment_table(cfg)
    summary = comparative.substitution_summary(table)
    track = comparative.cox1_local_dG_track(cfg)
    flux = comparative.directional_flux_table(cfg)
    titv = comparative.transition_transversion_summary(table)

    canonical_csv(table, out / "cox1_alignment_table.csv")
    canonical_csv(summary, out / "cox1_substitution_summary.csv")
    canonical_csv(track, out / "cox1_local_dG_track.csv")
    tables_csv(summary, cfg.outdir, "cox1_substitution_summary")
    tables_csv(track, cfg.outdir, "cox1_local_dG_track")
    if not flux.empty:
        canonical_csv(flux, out / "cox1_directional_flux.csv")
        tables_csv(flux, cfg.outdir, "cox1_directional_flux")
        comparative_plot.plot_directional_flux(
            flux, plot_path(out, "cox1_directional_flux_heatmap", cfg.plot_format), dpi=cfg.dpi
        )
    if not titv.empty:
        canonical_csv(titv, out / "cox1_transition_transversion.csv")
        tables_csv(titv, cfg.outdir, "cox1_transition_transversion")

    comparative_plot.plot_dG_track_species(
        track, "Yeast", plot_path(out, "cox1_dG_track_yeast", cfg.plot_format), dpi=cfg.dpi,
    )
    comparative_plot.plot_dG_track_species(
        track, "Human", plot_path(out, "cox1_dG_track_human", cfg.plot_format), dpi=cfg.dpi,
    )
    comparative_plot.plot_substitution_summary(
        summary, plot_path(out, "cox1_substitution_heatmap", cfg.plot_format), dpi=cfg.dpi,
    )
    return 0
