# Figure outputs

Reference for every figure that `mtrnafeat` writes. Each entry lists the
output filename, what is plotted, how to read each axis, and which CSV
holds the underlying numeric data (so figures can be re-rendered or
re-styled without rerunning the analysis).

All plots default to `.svg` (editable in Inkscape / Illustrator). Switch
to `png` or `pdf` by setting `plot_format` in the YAML config.

---

## stats

### `stats/stats_summary.svg`
Per-transcript boxplot grid summarizing length, MFE (kcal/mol),
foldedness (paired-fraction in the DMS-derived structure), GC fraction,
AU fraction, and median paired-pair span.

- **x-axis (per panel):** species (Human, Yeast).
- **y-axis (per panel):** the metric named in the panel title.
- **Each point:** one transcript.
- **CSV:** [`stats/per_transcript_statistics.csv`](../runs/) — one row
  per (Species, Gene).

---

## landscape

### `landscape/landscape_overlay.svg`
Per-sample simulation overlay against experimental DMS data.

- **x-axis:** normalized de novo DMS / simulated MFE
  (kcal/mol/nt).
- **y-axis:** structured percentage (% nucleotides paired).
- **Filled clouds:** simulated null structures for each sample.
- **Markers:** DMS-derived transcripts; ΔG is recomputed by ViennaRNA
  `eval_structure`, not taken from the `.db` header.
- **CSV:** `landscape/specific_conditions.csv`,
  `landscape/experimental_overlay.csv`.

### `landscape/landscape_overlay_{sample}.svg`
Single-sample version of the overlay above, with sample-specific null
cloud plus configured GC-reference clouds.

- **CSV:** `landscape/specific_conditions.csv`,
  `landscape/experimental_overlay.csv`.

### `landscape/landscape_overlay_{yeast_sample}_regions.svg`
Region-level overlay for yeast-annotated samples.

- **x-axis:** normalized de novo DMS ΔG for the projected region
  structure.
- **y-axis:** structured percentage in that region.
- **Contours:** UTR and CDS region-specific composition nulls when
  `landscape/region_mode_nulls.csv` is available.
- **Markers:** separate 5'UTR, CDS, and 3'UTR/tail points per gene.
- **CSV:** `landscape/experimental_overlay_regions.csv`,
  `landscape/region_mode_nulls.csv`.

### `landscape/yeast_region_folding_modes_{yeast_sample}.svg`
Focused yeast architecture figure.

- **Top row:** region-specific thermodynamic landscape
  (normalized DMS ΔG vs structured percentage).
- **Bottom row:** nucleotide-composition contour
  (region GC% vs C/(G+C)).
- **Contours:** simulated region nulls preserving each region's empirical
  A/U/G/C composition; sampled lengths are capped at `sim_seq_length`.
- **Dots:** per-gene DMS 5'UTR, CDS, and 3'UTR/tail regions.
- **CSV:** `landscape/region_mode_nulls.csv`,
  `landscape/experimental_overlay_regions.csv`.

### `landscape/pairing_bias_{GC,AU,GU}.svg`
Per-base-pair-type composition relative to background.

- **x-axis:** simulated GC fraction.
- **y-axis:** fraction of paired positions whose partner is the named
  base-pair class.
- **CSV:** `landscape/gc_gradient.csv`, `landscape/experimental_overlay.csv`.

### `landscape/paired_nt_fractions_{sample}.svg`
Three-panel view of absolute paired-nucleotide budget.

- **Panel A:** foldedness (% nt paired).
- **Panel B:** G-C paired nucleotides as % of sequence plus theoretical
  composition ceiling.
- **Panel C:** A-U paired nucleotides as % of sequence plus ceiling.
- **Legend:** outside the plot boundary.
- **CSV:** `landscape/gc_gradient.csv`,
  `landscape/gc_gradient_biased_heavy.csv`,
  `landscape/experimental_overlay.csv`.

---

## features

### `features/heatmap_size_ratios.svg`
Element-size composition heatmap.

- **x-axis:** Simulated and DMS motif columns for each structural element.
- **y-axis:** element size, capped by `max_heatmap_size`.
- **Color:** percent enrichment within each sample/type/motif group.
- **CSV:** `features/raw_motifs.csv`.

### `features/phase_space_contour.svg`
2-D contour of helix-size vs loop-size feature space.

- **x-axis:** average macro-helix size per transcript.
- **y-axis:** average loop size per transcript
  (hairpin + bulge + internal loop).
- **Contours/clouds:** simulated sample-specific nulls.
- **Markers:** DMS-derived transcripts labelled by gene.
- **CSV:** derived from `features/raw_motifs.csv`.

### `features/span_boxplot.svg`
Base-pairing distance ECDF.

- **x-axis:** base-pairing distance `|i-j|` in nt, log scale.
- **y-axis:** cumulative fraction of pairs.
- **Lines:** Sim vs DMS per sample; median guide lines are shown.
- **CSV:** `features/raw_spans.csv`.

---

## window

### `window/window_{species}_{gene}.svg`
Whole-transcript fold-and-compare trace per gene.

- **x-axis:** transcript position (1-based).
- **Top panel:** smoothed paired fraction for DMS-derived dot-bracket
  vs the configured span-limited folding engine.
- **Bottom strip:** transcript architecture (5′UTR / CDS / 3′UTR) when
  annotations are available.
- **CSV:** `window/window_per_position.csv` (long-format per gene).

---

## local-probability

### `local_probability/local_probability_{species}_{gene}.svg`
Four-panel per-gene plot when DMS overlay is available; falls back to a
two-panel layout without DMS.

- **Panel 1 (RNAplfold P(paired)):** raw + smoothed local pair-probability
  track from RNAplfold. y in [0, 1].
- **Panel 2 (DMS paired fraction):** DMS-derived binary paired indicator
  smoothed with the same window. y in [0, 1].
- **Panel 3 (signed Δ):** per-window mean Δ = mean P(paired) − mean DMS
  paired fraction at the configured `local_probability_scan_window_nt`
  scale. Positive = thermodynamically more pairable than experimentally
  observed; negative = experimentally more paired than locally
  thermodynamic.
- **Panel 4 (architecture):** 5′UTR / CDS / 3′UTR strip with the TIS
  shaded; vertical lines at start and stop codons. The TIS shading
  spans all four panels.
- **CSVs:** `local_probability_per_position.csv` (panels 1–3 raw),
  `local_probability_per_window.csv` (panel 3 windowed).

### `local_probability/local_probability_TIS_summary.csv` (no figure)
TIS vs CDS-background effect size table with circular-shift empirical
p-values and Benjamini-Hochberg q-values across genes. See README §1.1.

### `local_probability/local_probability_TIS_sensitivity.csv` (no figure)
Same metrics at multiple TIS context widths (default ±30, ±50, ±100,
±200, ±500). BH-FDR is applied within each `TIS_Window_Tag` (each
context width is its own hypothesis family). Use for robustness checks
when the standard ±50/±50 window may hide signal in transcripts with
unusually long 5′UTRs (e.g. yeast COX1).

---

## structure-deviation

### `structure_deviation/structure_deviation_{species}_{gene}.svg`
Four-panel per-gene plot showing where the global Vienna MFE model
disagrees with the DMS-derived structure.

- **Panel 1 (global-MFE paired fraction):** smoothed paired-binary
  track from the whole-transcript Vienna MFE structure.
- **Panel 2 (DMS paired fraction):** smoothed paired-binary track from
  the DMS dot-bracket.
- **Panel 3 (Δ track):** signed deviation `P_model − P_DMS` at the
  same per-window scale used by the plot/table aggregation. Horizontal band at
  ±`structure_deviation_threshold` (default ±0.25).
- **Panel 4 (called regions + architecture):** per-region rectangles colored by
  `Region_Class` (`model_high_dms_low` = magenta open-but-foldable;
  `model_low_dms_high` = teal protected-beyond-model;
  `concordant_paired` / `concordant_open` = neutral; `mixed_deviation`,
  `ambiguous` = gray), overlaid on the transcript architecture strip.
- **CSVs:** `structure_deviation_per_position.csv` (panels 1–3),
  `structure_deviation_regions.csv` (panel 4).

### `structure_deviation/structure_deviation_lollipop_{species}.svg`
Per-species summary of every called region.

- **x-axis:** transcript position (1-based) of the region midpoint.
- **y-axis:** signed mean deviation. Stem heights are `Mean_Deviation`;
  positive stems are `model_high_dms_low`, negative stems are
  `model_low_dms_high`.
- **Color:** `Region_Class`.
- **Faceted by gene:** one row per gene.
- **CSV:** `structure_deviation_regions.csv`.

### `structure_deviation/structure_deviation_heatmap.svg`
Cross-gene architectural-bin heatmap.

- **x-axis:** architectural bin (`5_end`, `TIS`, `early_CDS`,
  `mid_CDS`, `late_CDS`, `stop_proximal`, `3_end`).
- **y-axis:** species × gene.
- **Color:** mean signed deviation per bin (red = positive, blue =
  negative). Empty bins (e.g. transcripts without UTRs) are masked.
- **CSV:** `structure_deviation_gene_region_matrix.csv`.

### Statistical interpretation of the regions table
When `structure_deviation_null_model: dinuc` and
`structure_deviation_n_null > 0` are set in the YAML config (defaults
are conservative `none` / `0` to keep `run-all` fast), each region in
`structure_deviation_regions.csv` carries:

- `Empirical_P` — one-sided permutation p-value against the per-gene
  null distribution of max-|deviation| under Altschul-Erikson
  dinucleotide shuffle (Westfall-Young max-statistic correction within
  gene; family-wise error controlled at the per-gene level).
- `Q_Value` — Benjamini-Hochberg FDR q-value across all called regions
  pooled across genes.
- `Statistical_Support_Label` — `max_stat_significant` (p < 0.05) or
  `max_stat_nonsignificant`. Falls back to `effect_size_only` when no
  null was run.

---

## tis

### `tis/tis_zoom_grid.svg`
Per-gene TIS energy comparison.

- **x-axis:** gene, with `*` when the requested upstream context is
  truncated by a short 5'UTR.
- **y-axis:** ΔG (kcal/mol).
- **Bars:** DMS-derived projected TIS window vs Vienna MFE on the same
  sequence window.
- **CSV:** `tis/tis_dms_vs_mfe.csv`.

---

## substitution

### `substitution/substitution_kde_panels_{species}.svg`
Per-species KDE panels of ΔG distributions for each null pool.

- **x-axis:** ΔG (kcal/mol).
- **y-axis:** kernel density.
- **Lines:** three null pools (flat-ACGU, positional-ACGU, synonymous).
- **Vertical line:** wild-type observed ΔG.
- **Faceted by gene:** one panel per gene.
- **CSV:** `substitution/substitution_thermo_distribution.csv`.

### `substitution/substitution_z_heatmap_{species}.svg`
Z-score heatmap of (gene × null pool).

- **Cell value:** standardized z-score of the wild-type ΔG against the
  null pool (more negative = wild-type more stable than expected).
- **CSV:** `tables/substitution_thermo_summary.csv`.

### `substitution/substitution_effect_shift_{species}.svg`
Effect-size companion to the density plots.

- **x-axis:** ΔΔG per nt versus the null-pool mean.
- **y-axis:** gene.
- **Panels:** flat-ACGU, positional-ACGU, and synonymous nulls.
- **Markers:** WT Vienna MFE and DMS-structure ΔG; connecting segments
  show how far the realized DMS structure sits from the WT MFE reference.
- **Shading:** negative ΔΔG, where the reference is more stable than the
  null mean.
- **CSV:** `tables/substitution_thermo_summary.csv`.

---

## cofold

### `cofold/cofold_gap_closure.svg`
Gap-closure summary for the CoFold (α, τ) sweep.

- **x-axis:** α penalty strength.
- **y-axis:** fraction of the plain-Vienna-to-DMS ΔG gap closed
  (0 = plain Vienna, 1 = matches DMS-evaluated ΔG).
- **Lines:** τ values; shading = ±1 SD across genes.
- **CSV:** `cofold/cofold_grid.csv`; per-gene best in
  `cofold/cofold_best_per_gene.csv`.

### `cofold/cofold_parameter_landscape.svg`
Species-level α/τ heatmap of the CoFold sweep.

- **x-axis:** α (penalty strength, kcal/mol).
- **y-axis:** τ (decay length, nt).
- **Color:** mean |CoFold − DMS| ΔG gap (kcal/mol), aggregated over genes.
- **CSV:** `cofold/cofold_grid.csv`.

### `cofold/cofold_parameter_landscape_{species}.svg`
Per-gene α/τ heatmaps of the CoFold sweep.

- **x-axis:** τ (decay constant, nt).
- **y-axis:** α (penalty strength, kcal/mol).
- **Color:** |ΔG_cofold − ΔG_DMS| (kcal/mol).
- **CSV:** `cofold/cofold_grid.csv`.

`cofold/cofold_per_window_corr.csv` is a CSV-only diagnostic containing
per-window Pearson correlations; no per-window cofold figure is currently
emitted.

---

## compare

### `compare/cox1_substitution_heatmap.svg`
Yeast↔human COX1 codon-aligned substitution heatmap.

- **x-axis:** alignment column (codon-level).
- **y-axis:** substitution type (synonymous, non-synonymous, gap).
- **Cell value:** count.
- **CSV:** `compare/cox1_alignment_table.csv`.

### `compare/cox1_directional_flux_heatmap.svg`
Directional substitution flux Yeast→Human and Human→Yeast.

- **Cell value:** signed flux per codon position class.
- **CSV:** `compare/cox1_directional_flux.csv`.

---

## gene-panel

### `gene_panels/panel_{species}_{gene}.svg`
Per-gene composite panel.

- **Panel A:** base composition counts.
- **Panel B:** paired-pair composition of the DMS-derived structure.
- **Panel C:** rolling local paired fraction along the transcript, with
  an architecture strip when annotations are available.
- **CSV:** built directly from the input `.db` records.

---

## kinetic (opt-in)

### `kinetic/kinetic_{species}_{gene}.svg`
DrTransformer co-transcriptional folding trajectory plot.

- **x-axis:** transcript length (nt) during transcription.
- **y-axis:** macrostate occupancy (0–1).
- **Lines:** per-macrostate occupancy over time.
- **CSV:** `kinetic/kinetic_summary.csv`.

This stage requires a working DrTransformer runtime and is never auto-run
by `run-all`.
