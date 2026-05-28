"""`mtrnafeat landscape` — simulated GC-gradient + experimental overlay + plots."""
from __future__ import annotations

from mtrnafeat.analysis import landscape
from mtrnafeat.config import Config
from mtrnafeat.constants import file_safe_sample
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import landscape_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "landscape"
    out.mkdir(parents=True, exist_ok=True)

    sim_df = landscape.simulate_specific_conditions(cfg)
    grad_df = landscape.simulate_gradient(cfg)
    biased_df = landscape.simulate_biased_gradient(cfg)
    biased_heavy_df = landscape.simulate_biased_gradient(cfg, heavy_strand=True)
    exp_df = landscape.experimental_overlay(cfg)
    exp_region_df = landscape.experimental_overlay_regions(cfg)

    canonical_csv(sim_df, out / "specific_conditions.csv")
    canonical_csv(grad_df, out / "gc_gradient.csv")
    canonical_csv(biased_df, out / "gc_gradient_biased.csv")
    canonical_csv(biased_heavy_df, out / "gc_gradient_biased_heavy.csv")
    canonical_csv(exp_df, out / "experimental_overlay.csv")
    if not exp_region_df.empty:
        canonical_csv(exp_region_df, out / "experimental_overlay_regions.csv")

    landscape_plot.landscape_overlay(sim_df, exp_df, plot_path(out, "landscape_overlay", cfg.plot_format), dpi=cfg.dpi)
    landscape_plot.gradient_curves(grad_df, plot_path(out, "gc_gradient_curves", cfg.plot_format), dpi=cfg.dpi)
    landscape_plot.pairing_bias(grad_df, exp_df, plot_path(out, "pairing_bias_GC", cfg.plot_format),
                                  y_col="Paired_GC_Pct", ylabel="Paired G-C (%)", include_yx_line=True, dpi=cfg.dpi)
    landscape_plot.pairing_bias(grad_df, exp_df, plot_path(out, "pairing_bias_AU", cfg.plot_format),
                                  y_col="Paired_AU_Pct", ylabel="Paired A-U (%)", include_yx_line=False, dpi=cfg.dpi)
    landscape_plot.pairing_bias(grad_df, exp_df, plot_path(out, "pairing_bias_GU", cfg.plot_format),
                                  y_col="Paired_GU_Pct", ylabel="Paired G-U wobble (%)", include_yx_line=False, dpi=cfg.dpi)

    # Per-species H-strand-biased pairing curves: baseline preserves empirical
    # C/(G+C) and A/(A+U) ratios from H-strand transcripts only (drops human
    # ND6, which is L-strand-encoded and reverses the skew). Symmetric baseline
    # is overlaid as a dashed reference so the bias-induced shift is visible.
    heavy_freqs = landscape.species_freqs_for_pipeline(cfg, heavy_strand=True)
    for species in cfg.db_files:
        sp_lc = file_safe_sample(species).lower()
        landscape_plot.pairing_bias_species(
            biased_heavy_df, exp_df,
            plot_path(out, f"pairing_bias_GC_{sp_lc}", cfg.plot_format),
            species=species, y_col="Paired_GC_Pct", ylabel="Paired G-C (%)",
            include_yx_line=True, symmetric_gradient_df=grad_df,
            species_freqs=heavy_freqs.get(species), dpi=cfg.dpi,
        )
        landscape_plot.pairing_bias_species(
            biased_heavy_df, exp_df,
            plot_path(out, f"pairing_bias_AU_{sp_lc}", cfg.plot_format),
            species=species, y_col="Paired_AU_Pct", ylabel="Paired A-U (%)",
            include_yx_line=False, symmetric_gradient_df=grad_df,
            species_freqs=heavy_freqs.get(species), dpi=cfg.dpi,
        )
        # Sequence-fraction view: shows the foldedness drop driven by the
        # composition ceiling 2·min(G,C) and 2·min(A,U), which the pair-type
        # ratio metric above hides.
        landscape_plot.paired_nt_fractions(
            grad_df, biased_heavy_df, exp_df,
            plot_path(out, f"paired_nt_fractions_{sp_lc}", cfg.plot_format),
            species=species, species_freqs=heavy_freqs.get(species), dpi=cfg.dpi,
        )

    # Per-species separate overlay files (species-filtered simulation cloud
    # plus GC reference contours overlaid).
    for species in cfg.db_files:
        landscape_plot.landscape_overlay_one(
            sim_df, exp_df,
            plot_path(out, f"landscape_overlay_{file_safe_sample(species).lower()}", cfg.plot_format),
            species=species, dpi=cfg.dpi,
        )
        if not exp_region_df.empty and (
            not exp_region_df[(exp_region_df["Species"] == species)
                              & (exp_region_df["Annotation_Species"] == "Yeast")].empty
        ):
            landscape_plot.landscape_overlay_regions(
                sim_df, exp_region_df,
                plot_path(out, f"landscape_overlay_{file_safe_sample(species).lower()}_regions", cfg.plot_format),
                species=species, dpi=cfg.dpi,
            )

    # Per-nucleotide composition bias plot (replaces the previous violin).
    for species in cfg.db_files:
        landscape_plot.per_base_composition_bias(
            biased_df, exp_df,
            plot_path(out, f"nucleotide_bias_{file_safe_sample(species).lower()}", cfg.plot_format),
            species=species, dpi=cfg.dpi,
        )
    return 0
