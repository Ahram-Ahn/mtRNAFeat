# Changelog

All notable changes to this project will be documented in this file. Format
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioned
per [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **`structure-deviation` dinucleotide-shuffle null model** (opt-in via
  `structure_deviation_null_model: dinuc` + `structure_deviation_n_null > 0`,
  also exposed as `--null dinuc --n-null N` on the CLI). Per-gene null
  distribution of max-|deviation| under Altschul-Erikson dinucleotide
  shuffle (Westfall-Young max-statistic correction within gene), then
  Benjamini-Hochberg q-values across all called regions pooled across
  genes. Each row in `structure_deviation_regions.csv` now carries
  `Empirical_P`, `Q_Value`, `Null_Model`, `N_Null`, and a
  `Statistical_Support_Label` of `max_stat_significant`,
  `max_stat_nonsignificant`, or `effect_size_only` (when null is off).
  Per-gene summary gains `N_Significant_Regions_Q05`. Defaults stay
  conservative (`null_model: "none"`, `n_null: 0`) so `run-all`
  wall-time is unchanged unless the user opts in.
- **Per-run reproducibility manifest** (`<outdir>/run_manifest.json`):
  written by every `mtrnafeat` invocation. Records the package version,
  resolved ViennaRNA and RNAstructure versions/paths, current git
  commit (with `(dirty)` marker when local edits are present), Python
  version, platform, the command name and its full argv, the seed,
  the config path, and an ISO-8601 UTC timestamp. Lets reviewers and
  downstream consumers reproduce a run with the same tool stack and
  audit drift between runs (e.g. silent ViennaRNA upgrades).
- **CI workflow** (`.github/workflows/test.yml`): two-job matrix that runs
  unit tests + ruff on Python 3.11/3.12 without ViennaRNA (skipping
  `needs_rna` tests), and a separate Conda-based job that installs
  ViennaRNA + RNAstructure via `environment.yml` and runs the full suite
  plus the `examples/01_smoke_mini.sh` smoke test.
- **`CITATION.cff`** at the repo root for machine-readable citation
  metadata (GitHub auto-renders this on the repo page).
- **`pyproject.toml` metadata**: `[project.urls]` (homepage, repository,
  issues, changelog), trove classifiers, keywords, and author email.
- **`Config.__post_init__` validation**: bounds and enum checks on every
  numeric / categorical field. Out-of-range values (negative thresholds,
  zero windows, `plot_format: jpg`, `fold_engine: mfold`, etc.) now
  raise `ValueError` at config load instead of silently producing
  NaN tracks downstream. Cross-field invariant checks include
  `rnaplfold_max_bp_span <= rnaplfold_window` and
  `structure_deviation_low_threshold <= structure_deviation_high_threshold`.
- **BH-FDR multiple-testing correction** on TIS empirical p-values
  (`mtrnafeat.core.stats.bh_fdr`). Each `Empirical_P_*` column in
  `local_probability_TIS_summary.csv` now ships with a parallel
  `Q_Value_*` column. The sensitivity-sweep table corrects within each
  TIS context width (each `(upstream, downstream)` pair is its own
  hypothesis family); the primary summary corrects across genes.
- **New config fields** for `structure-deviation` bin definitions:
  `structure_deviation_mid_cds_lo` (default 0.30),
  `structure_deviation_mid_cds_hi` (default 0.70), and
  `structure_deviation_late_cds_window_nt` (default 300). Lifted from
  hard-coded values inside `_bin_intervals` so the cross-gene heatmap
  bin definitions are now visible and tunable from YAML.

### Changed

- **`README.md`**: replaced the `<repo-url>` placeholders in the install
  block and the bibtex citation with the actual repository URL. Added a
  link to the new `docs/FIGURES.md`.
- **CLI `_parse()` helpers reject unknown flags.** `local-probability`,
  `structure-deviation`, `cofold`, `significance`, `substitution`,
  `window`, and `run-all` now raise `SystemExit` with a clear message on
  any unrecognized flag instead of silently ignoring it (a typo like
  `--threshhold 0.3` no longer slips through with the default value).

### Documentation

- **`docs/FIGURES.md`**: new reference describing every figure output
  (what is plotted, how to read each axis, source CSV) — required for
  publication supplementary methods and for users re-rendering plots
  from cached CSVs.

## [0.3.0] — 2026-05-02

### Added

- **`structure-deviation` stage**: region-discovery analysis on top of
  the RNAplfold-vs-DMS comparison. Calls intervals where the smoothed
  signed deviation `P_model − P_DMS` exceeds a threshold, classifies
  each region into one of `model_high_dms_low`, `model_low_dms_high`,
  `concordant_paired`, `concordant_open`, `mixed_deviation`, or
  `ambiguous`, and emits four publication-oriented CSVs (per-position,
  regions, gene-summary, gene-region-matrix) plus per-gene four-panel
  plots, per-species lollipop summaries, and a cross-gene
  architectural-bin heatmap. Replaces the biological-interpretation
  role of the legacy `significance --scan`. New CLI:
  `mtrnafeat structure-deviation --config configs/all.yaml --outdir runs/dev`.
  Wired into `run-all` alongside the existing stages.
- **`local-probability` plot polish**: TIS shading now runs vertically
  through all four panels; the per-window Δ panel replaces the
  per-position smoothed Δ (cleaner on long transcripts); a context
  subtitle (`5'UTR=… · CDS=… · TIS=…`) sits below every title; raw
  P(paired) fades on transcripts > 1000 nt.
- **`local-probability` TIS sensitivity sweep**: new
  `local_probability_TIS_sensitivity.csv` reporting the same TIS
  metrics at a configurable sweep of `(upstream, downstream)` widths
  (defaults `(30,30), (50,50), (100,100), (200,200), (500,500)`). A
  single fixed window hides signal in long-5'UTR yeast genes
  (Yeast COX1 is unusually open at `−30/+30`, P ≈ 0.04, but trivial
  at `−100/+100`); reviewers can read robustness from this table
  directly.
- **`local-probability` DMS overlay**: per-gene plots now render four
  panels (RNAplfold P(paired), DMS paired fraction, signed Δ,
  architecture+TIS) and emit `local_probability_per_window.csv`
  (windowed agreement) and `local_probability_TIS_summary.csv`
  (TIS-vs-CDS-background effect size with circular-shift empirical
  p-values). The per-position CSV now carries `DMS_Paired_*` columns
  whenever the matching `.db` record's length agrees with the
  RNAplfold input.

### Changed

- **`run-all` no longer invokes `significance`** by default.
  `structure-deviation` is the publication-facing biological-interpretation
  layer; `significance` is preserved as a standalone command for users
  who want the thermodynamic-null QC. Removing it from the default
  pipeline cuts wall-clock time substantially (the dinuc-shuffle null
  pool was the slowest stage). Run `mtrnafeat significance --config
  configs/all.yaml --outdir runs/sig -- --scan` directly when needed.
- **`structure-deviation` Δ panel harmonized with `local-probability`**.
  Both stages now display the deviation track at the same
  per-`local_probability_scan_window_nt`-window scale (default 120 nt,
  step 10) so the two figures tell visually identical stories. Region
  calling itself still runs on the per-position 25-nt rolling track in
  the analysis layer; the rectangles in the Δ panel carry that
  information.
- **`significance` plot/CSV relabeled** so within-gene Z-scores are
  no longer presented as p-values. The cotrans plot title now reads
  *"local structural-change scan (mode=…, n=N windows, candidate
  peaks |within-gene Z| ≥ T)"*; bottom-panel y-label is "Within-gene
  Z (smoothed Δ)"; peak-marker labels carry "candidate peak". Every
  row of `cotrans_per_window.csv` is tagged with
  `Z_score_type="within_gene_window_standardized_delta"` and
  `Is_statistical_pvalue=False`. The Workman-Krogh `z_per_gene.csv`
  is the only null-model-backed output of this command; layer-2
  scan is exploratory.
- **README + STAGES.md** rewritten to describe the new
  `structure-deviation` stage as the publication-facing
  interpretation layer; the legacy `significance` is documented as
  the optional thermodynamic-null QC.
- **`window` plot now shows two traces, not three.** The
  unconstrained ``Vienna full`` trace is computed and stored in the
  per-position CSV but no longer plotted; for long mRNAs the
  configured ``max_bp_span`` (default 300 nt) already captures every
  realistic contact. The legend title became "DMS vs <engine>" so the
  comparison reads directly.
- **`tis` plot refactored.** X-tick labels are now a single line
  (just the gene name with a ``*`` suffix when the 5'UTR was
  truncated); the per-gene 5'UTR length moved to a small italic
  annotation under the bar group instead of being baked into the
  tick label. Legend wording shortened ("DMS-derived",
  "Vienna prediction", "* 5'UTR truncated", "ΔG = 0").
- **`cofold` plot replaces the per-gene heatmap grid with a single
  per-species strip plot** of |CoFold − DMS| across the (α, τ) sweep.
  Genes are sorted by best-fit gap; each gene's dots are colored by
  τ; a red ring marks the best-fit dot per gene; a black square marks
  the published-default α=0.5 / τ=640. The full sweep is still in
  ``cofold_grid.csv`` for users who want the (α, τ) surface. Filename
  changed: ``cofold_gap_heatmap_{species}.{ext}`` →
  ``cofold_gap_strip_{species}.{ext}``. The
  ``gap_heatmap_panels`` Python entry-point is kept as a backward
  alias for ``gap_strip_panels``.

### Fixed

- **Substitution KDE panels: WT vertical lines could fall off the
  auto-scaled x-axis** when the wild-type ΔG was far from the null
  pool's mean (e.g. yeast ATP9: WT_MFE = −68.9 kcal/mol vs pool means
  near −38). The KDE x-limits are now expanded explicitly to include
  both WT lines so the WT-vs-pool comparison is visible for every
  gene.

## [0.2.0] — 2026-05-01

Phase 1 of the realignment toward a careful downstream analysis layer
(see `mtrnafeat_ai_refinement_instructions.txt`). User-facing changes
focus on environment / input pre-flight checks and documentation
terminology. No analysis-stage outputs change in this release.

### Added

- **Two new substitution null pools, `flat_acgu` and `positional_acgu`,**
  that preserve the full A/C/G/T frequency vector (and per-codon-position
  variant of it) instead of pooling G with C and A with T. Critical for
  transcripts with marked nucleotide-pair asymmetry — e.g. human ND6
  (L-strand) has G≈191 vs C≈37, which the overall-GC null hides
  entirely. Existing `flat_gc`, `positional_gc`, and `synonymous` pools
  are unchanged; the new pools draw from the same RNG after the
  existing pools so historical bit-stable outputs are preserved.
- **`viz.style.legend_outside()`** — shared helper that places legends
  outside the data axes (right / bottom / top) so they never occlude
  plotted data. Replaces ad-hoc `loc="best"` and `loc="upper right"`
  calls in `substitution`, `features`, `local_probability`, `cofold`,
  and `cotrans` plots.
- **`mtrnafeat doctor`** — environment diagnostics command. Reports
  Python version, ViennaRNA importability, RNAplfold and RNAstructure
  availability (the latter conditionally `ERROR` when
  `fold_engine: rnastructure`), `DATAPATH`, and DrTransformer
  (`OPTIONAL`). Supports `--json PATH` for machine-readable output.
- **`mtrnafeat validate-inputs`** — pre-flight check of config, `.db`
  files (length parity, bracket balance, ACGU alphabet), bundled
  annotations (UTR/CDS coordinates within transcript length, CDS
  divisibility by 3), and optional inputs (alignment, modifications
  table). Supports `--json PATH`.
- **`src/mtrnafeat/validation.py`** — reusable `ValidationIssue`
  dataclass plus pure-function checkers used by both new commands.
- README sections **"Recommended workflow"**, **"What mtrnafeat is
  not"**, and **"Interpretation cautions"**, framing the package's
  scope and limitations explicitly.

### Changed

- **Substitution `Pool` enum gained `flat_acgu` and `positional_acgu`
  values.** The CSV schema is unchanged (still long-format with `Pool`
  column); per-gene Z-heatmaps and KDE panels now show 5 pools instead
  of 3.
- **Plot legends are now placed outside the data axes** for every
  figure that would otherwise occlude plotted content
  (`substitution_kde_panels_*`, `phase_space_contour`, `span_boxplot`,
  `local_probability_*`, `cofold_per_window_corr_*`, `cotrans_*`).
  Plots that already used `bbox_to_anchor` (window, tis, landscape,
  kinetic, gene-panel) are unchanged. The two-species heatmap title
  in `heatmap_size_ratios` no longer overlaps thanks to a `\n` line
  break and explicit `subplots_adjust`.
- Documentation now refers to the `.db` dot-brackets as
  **DMS-derived** or **DMS-constrained model** rather than
  **"ground truth"**, reflecting that these structures are model
  outputs from DMS-MaPseq data, not direct experimental observations.
  Affected files: `README.md`, `docs/STAGES.md`, `docs/CONFIG.md`,
  `configs/all.yaml`, `configs/template.yaml`, `src/mtrnafeat/config.py`,
  `CHANGELOG.md`.
- README opening rewritten around the six analysis questions the
  package is meant to answer, with explicit non-goals up front.

## [0.1.0-pre-rename] — pre-0.2.0 development

Items previously listed under `[Unreleased]`:

### Added

- **`local-probability` stage** — new subcommand that runs ViennaRNA's
  RNAplfold per-position pair-probability scan over each transcript, with
  a smoothed plot overlay. Wired into `run-all` so it ships in the full
  pipeline. New config fields `rnaplfold_window` (default 80),
  `rnaplfold_max_bp_span` (default 50), `rnaplfold_cutoff` (default 0.001),
  and `rolling_window` (default 25; reused by the `window` plot smoother
  too).
- **Cotranscriptional / sliding scan in `significance --scan`** —
  per-gene scan that walks the transcript in `prefix` (nascent-style) or
  `sliding` mode and z-scores Vienna MFE / ensemble diversity / paired
  fraction. Emits `cotrans_per_window.csv` and one PNG per (species,
  gene). `--per-window` is preserved as a deprecated alias so `run-all`
  keeps working.
- **`WT_DMS_Eval_FullLength` reference column in substitution output** —
  full-length Vienna `eval_structure` of the `.db` dot-bracket reported
  alongside the chunk-truncated `WT_DMS_Eval`, plus a `DMS_dG_Source`
  audit column making provenance explicit.
- **Architecture strip in `gene-panel`** — third panel now uses a
  nested grid with a dedicated 5'UTR / CDS / 3'UTR strip below the
  paired-fraction trace; consistent label treatment (inside wide bands,
  leader-lined when narrow).
- **ΔG = 0 vs missing-data disambiguation in `tis`** — open-diamond
  marker plus literal `0` annotation when the projected DMS structure
  has no surviving pairs in the TIS window; explicit `n/a` annotation
  for true NaN.

### Changed

- **`features/span_boxplot`** redesigned as per-species ECDF on log-x
  axis with median dropdown lines. Heavy-tailed pairing-distance
  distributions are now legible.
- **`window` plot legend** moved outside the data axis (figure-coord
  bbox), so the three pairing-fraction traces are never occluded on
  long or peak-rich transcripts.
- **`compare/` no longer emits the per-codon local-ΔG track** (the
  per-species `cox1_dG_track_*.svg` and `cox1_local_dG_track.csv`).
  The `window` stage covers local-ΔG along the transcript at higher
  resolution and with both folding engines, so the duplicate panel
  was visually noisy without adding signal.
- **`fold_engine` default is `rnastructure`** (was `vienna` in the
  initial release). Matches how the upstream `.db` DMS-derived
  dot-brackets were produced. Override per run with `-- --engine vienna`
  if RNAstructure isn't on PATH.

### Fixed

- **Substitution stage: ΔG provenance was implicit and easy to misread.**
  Module docstring, plot title/legend, and the new `DMS_dG_Source`
  column now make it explicit that the wild-type DMS ΔG is recomputed
  by ViennaRNA's `eval_structure` on the `.db` dot-bracket — the `.db`
  header MFE value (e.g. `>COX1: -150.4 kcal/mol`) is intentionally
  bypassed.
- **Documentation drift across `docs/STAGES.md`, `docs/CONFIG.md`, and
  `README.md`** — wrong filenames for `window` (`window_per_position.csv`
  / `window_summary.csv`, not `window_scan_metrics.csv` / `_summary.csv`),
  `substitution` (`substitution_kde_panels_{species}.{ext}` / `_z_heatmap_
  {species}.{ext}`), `cofold` (`cofold_gap_heatmap_{species}.{ext}`),
  missing `local-probability` documentation, wrong claim that
  substitution pools fold under CoFold (they fold under plain Vienna),
  wrong claim that only ViennaRNA was wired through (RNAstructure has
  been the default for `window` for a while), missing `rolling_window`
  and `rnaplfold_*` config fields.

## [0.1.0-docs] — first published documentation set

### Added
- `configs/template.yaml` — annotated all-fields configuration template
  covering every field of the `Config` dataclass.
- `docs/CONFIG.md` — field-by-field configuration reference.
- `docs/STAGES.md` — per-subcommand reference (purpose, I/O, flags, notes).
- `CONTRIBUTING.md` — dev setup, tests, lint, conventions for adding
  config fields and stages.
- `CHANGELOG.md` — this file.

### Changed
- `README.md` rewritten with quickstart, stage overview table, and links
  into the new docs.

## [0.1.0] — initial release

### Pipeline stages

- `stats` — per-transcript statistics (length, MFE, foldedness, GC%,
  paired-pair composition) + boxplot.
- `landscape` — GC-gradient simulation against experimental DMS overlay,
  using empirical (A,U,G,C) frequencies by default and symmetric-GC nulls
  on configurable conditions.
- `features` — element decomposition (motifs, spans), region-stratified
  (5'UTR / CDS / 3'UTR) heatmaps, phase-space contour, span boxplots.
- `window` — sliding-window scan (default 120 nt / 10 nt step / 300 nt
  max bp span), DMS-vs-Vienna trace per gene.
- `significance` — dinucleotide-shuffle z-scores per gene, optionally
  per-window.
- `tis` — TIS zoom (−50/+50 nt around start codon). Window honors actual
  5'UTR length when present, clamps at the 5' end when absent. Emits
  `Has_5UTR` / `Has_Full_5UTR_Context` columns; truncated upstream is
  hatched in plots.
- `substitution` — synonymous-recoding thermodynamic permutation test
  (flat-GC, positional-GC, synonymous null pools) under plain Vienna MFE.
  *(The 0.1.0 release notes originally said "under CoFold" — that was
  inaccurate even at the time. Both the wild-type fold and the pool
  folds use `thermo.fold_mfe`. The CoFold soft-penalty exploration is
  the dedicated `cofold` sweep.)*
- `cofold` — CoFold (α, τ) parameter sweep against DMS ΔG, with optional
  per-window correlation pass.
- `compare` — yeast↔human COX1 codon-aligned comparative (substitution
  table, ΔG track, transition/transversion summary, directional flux).
- `gene-panel` — per-gene composition + paired-pair + local-foldedness panel.
- `kinetic` — DrTransformer kinetic folding (opt-in only; never auto-runs
  in `run-all`).
- `plot` — re-render plots from cached CSVs.
- `run-all` — orchestrates every independent stage; supports `--parallel`
  and `--skip`.

### Infrastructure

- Single `Config` dataclass + YAML loader; unknown keys raise.
- Deterministic CSV output (atomic write, fixed column order, fixed float
  format); SHA-pinned in `tests/test_determinism.py`.
- Per-run output layout: per-stage subdirectories plus a centralized
  `tables/` directory for cross-stage summary CSVs.
- Plot format selectable via YAML (`svg` default, `png`, `pdf`).
- Smoke test fixture under `test_data/` runs the full pipeline (minus
  `kinetic`) in under 2 minutes.

[Unreleased]: https://github.com/your-user/mtrnafeat/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/your-user/mtrnafeat/releases/tag/v0.1.0
