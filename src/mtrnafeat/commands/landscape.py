"""`mtrnafeat landscape` — simulated GC-gradient + experimental overlay + plots."""
from __future__ import annotations

from mtrnafeat.analysis import landscape
from mtrnafeat.config import Config
from mtrnafeat.io.writers import canonical_csv
from mtrnafeat.viz import landscape_plot
from mtrnafeat.viz.style import plot_path


def run(cfg: Config, args: list[str] | None = None) -> int:
    out = cfg.outdir / "landscape"
    out.mkdir(parents=True, exist_ok=True)

    sim_df = landscape.simulate_specific_conditions(cfg)
    grad_df = landscape.simulate_gradient(cfg)
    biased_df = landscape.simulate_biased_gradient(cfg)
    exp_df = landscape.experimental_overlay(cfg)

    canonical_csv(sim_df, out / "specific_conditions.csv")
    canonical_csv(grad_df, out / "gc_gradient.csv")
    canonical_csv(biased_df, out / "gc_gradient_biased.csv")
    canonical_csv(exp_df, out / "experimental_overlay.csv")

    landscape_plot.landscape_overlay(sim_df, exp_df, plot_path(out, "landscape_overlay", cfg.plot_format), dpi=cfg.dpi)
    landscape_plot.gradient_curves(grad_df, plot_path(out, "gc_gradient_curves", cfg.plot_format), dpi=cfg.dpi)
    landscape_plot.pairing_bias(grad_df, exp_df, plot_path(out, "pairing_bias_GC", cfg.plot_format),
                                  y_col="Paired_GC_Pct", ylabel="Paired G-C (%)", include_yx_line=True, dpi=cfg.dpi)
    landscape_plot.pairing_bias(grad_df, exp_df, plot_path(out, "pairing_bias_AU", cfg.plot_format),
                                  y_col="Paired_AU_Pct", ylabel="Paired A-U (%)", include_yx_line=False, dpi=cfg.dpi)
    landscape_plot.pairing_bias(grad_df, exp_df, plot_path(out, "pairing_bias_GU", cfg.plot_format),
                                  y_col="Paired_GU_Pct", ylabel="Paired G-U wobble (%)", include_yx_line=False, dpi=cfg.dpi)

    # Per-species separate overlay files (species-filtered simulation cloud
    # plus GC reference contours overlaid).
    for species in ["Human", "Yeast"]:
        landscape_plot.landscape_overlay_one(
            sim_df, exp_df,
            plot_path(out, f"landscape_overlay_{species.lower()}", cfg.plot_format),
            species=species, dpi=cfg.dpi,
        )

    # Per-nucleotide composition bias plot (replaces the previous violin).
    for species in ["Human", "Yeast"]:
        landscape_plot.per_base_composition_bias(
            biased_df, exp_df,
            plot_path(out, f"nucleotide_bias_{species.lower()}", cfg.plot_format),
            species=species, dpi=cfg.dpi,
        )
    return 0
