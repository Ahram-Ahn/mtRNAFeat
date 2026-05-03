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
GC-gradient simulation overlay against experimental DMS data.

- **x-axis:** simulated GC fraction (or per-species empirical GC).
- **y-axis:** mean MFE per nucleotide (kcal/mol/nt).
- **Lines:** simulated null pools at each GC step (median + IQR band).
- **Markers:** observed Human / Yeast transcripts.
- **CSV:** `landscape/gc_gradient.csv`.

### `landscape/pairing_bias_{GC,AU,GU}.svg`
Per-base-pair-type composition relative to background.

- **x-axis:** simulated GC fraction.
- **y-axis:** fraction of paired positions whose partner is the named
  base-pair class.
- **CSV:** `landscape/pairing_bias.csv`.

---

## features

### `features/heatmap_size_ratios.svg`
Element-size composition heatmap.

- **x-axis:** structural element type (helix length, hairpin loop,
  internal loop, multi-loop, bulge, free 5′/3′).
- **y-axis:** transcript region (5′UTR / CDS / 3′UTR) × species.
- **Color:** size ratio against the genome-wide pool (log2).
- **CSV:** `features/raw_motifs.csv` (raw counts), aggregated by region
  in the plotting code.

### `features/phase_space_contour.svg`
2-D contour of MFE vs paired-fraction.

- **x-axis:** paired fraction (0–1).
- **y-axis:** MFE per nucleotide (kcal/mol/nt).
- **Contours:** kernel density of all transcripts; overlay of region
  centroids (5′UTR / CDS / 3′UTR).
- **CSV:** derived from `features/raw_motifs.csv`.

### `features/span_boxplot.svg`
Pair-span distribution by region.

- **x-axis:** transcript region.
- **y-axis:** base-pair span in nucleotides (log scale).
- **CSV:** `features/raw_motifs.csv`.

---

## window

### `window/window_{species}_{gene}.svg`
Whole-transcript fold-and-compare trace per gene.

- **x-axis:** transcript position (1-based).
- **y-axis (top panel):** smoothed paired fraction; three lines for
  DMS-derived, Vienna full-fold, and configured-span engine fold.
- **y-axis (bottom panel):** signed Δ(paired fraction) between the
  engine-span fold and the DMS reference, with a band at ±0.25.
- **Architecture strip (under both panels):** 5′UTR (blue), CDS (gray),
  3′UTR (orange) with the start codon marked.
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
Four-panel per-gene plot showing where the local thermodynamic prior
disagrees with the DMS-derived structure.

- **Panel 1 (RNAplfold P(paired) + DMS paired fraction):** smoothed
  tracks for both, allowing visual concordance check.
- **Panel 2 (Δ track):** signed deviation `P_model − P_DMS` at the
  same per-window scale used by `local-probability` (so the two figures
  tell visually identical stories). Horizontal band at
  ±`structure_deviation_threshold` (default ±0.25).
- **Panel 3 (called regions):** per-region rectangles colored by
  `Region_Class` (`model_high_dms_low` = magenta open-but-foldable;
  `model_low_dms_high` = teal protected-beyond-model;
  `concordant_paired` / `concordant_open` = neutral; `mixed_deviation`,
  `ambiguous` = gray). The top-N regions by |Δ| are labeled with
  `Region_ID`.
- **Panel 4 (architecture + TIS):** as in `local-probability`.
- **CSVs:** `structure_deviation_per_position.csv` (panels 1–2),
  `structure_deviation_regions.csv` (panel 3).

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
Per-gene TIS −50/+50 nt zoom grid.

- **x-axis:** position relative to start codon (−50 to +50 nt).
- **y-axis:** smoothed paired fraction; one line for DMS-derived, one
  for engine MFE, both at the configured `rolling_window`.
- **Vertical line at 0:** start codon.
- **Faceted by (species, gene):** one panel per gene.
- **CSV:** `tis/tis_dms_vs_mfe.csv` (one row per (species, gene,
  position)).

---

## substitution

### `substitution/substitution_kde_{species}.svg`
Per-species KDE panels of ΔG distributions for each null pool.

- **x-axis:** ΔG (kcal/mol).
- **y-axis:** kernel density.
- **Lines:** five null pools (flat-GC, flat-ACGU, positional-GC,
  positional-ACGU, synonymous).
- **Vertical line:** wild-type observed ΔG.
- **Faceted by gene:** one panel per gene.
- **CSV:** `substitution/substitution_thermo_distribution.csv`.

### `substitution/substitution_z_heatmap.svg`
Z-score heatmap of (gene × null pool).

- **Cell value:** standardized z-score of the wild-type ΔG against the
  null pool (more negative = wild-type more stable than expected).
- **CSV:** `tables/substitution_thermo_summary.csv`.

---

## cofold

### `cofold/cofold_gap_strip_{species}.svg`
|Gap| strip plot of CoFold (α, τ) sweep.

- **x-axis:** τ (decay constant, nt).
- **y-axis:** |ΔG_cofold − ΔG_DMS| (kcal/mol).
- **Color:** α (penalty strength, kcal/mol).
- **Faceted by gene.**
- **CSV:** `cofold/cofold_grid.csv`; per-gene best in
  `cofold/cofold_best_per_gene.csv`.

### `cofold/cofold_per_window_corr_{species}.svg`
Per-window CoFold-vs-DMS correlation curves.

- **x-axis:** window center (nt).
- **y-axis:** Spearman correlation between CoFold-evaluated paired
  fraction and DMS paired fraction over the window.
- **CSV:** `cofold/cofold_per_window_corr.csv`.

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

- **Top panel:** GC and AU composition tracks.
- **Middle panel:** smoothed DMS-paired fraction.
- **Bottom panel:** local foldedness (paired-pair density) with an
  architecture strip (5′UTR / CDS / 3′UTR).
- **CSV:** built directly from the input `.db` records and the
  per-position outputs of `stats` and `local-probability`.

---

## kinetic (opt-in)

### `kinetic/kinetic_{species}_{gene}.svg`
DrTransformer co-transcriptional folding trajectory plot.

- **x-axis:** transcript length (nt) during transcription.
- **y-axis:** macrostate occupancy (0–1).
- **Lines:** per-macrostate occupancy over time.
- **CSV:** `kinetic/kinetic_summary.csv`.

This stage requires DrTransformer on `PATH` and is never auto-run by
`run-all`.
