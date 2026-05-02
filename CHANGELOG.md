# Changelog

All notable changes to this project will be documented in this file. Format
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioned
per [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] — 2026-05-01

Phase 1 of the realignment toward a careful downstream analysis layer
(see `mtrnafeat_ai_refinement_instructions.txt`). User-facing changes
focus on environment / input pre-flight checks and documentation
terminology. No analysis-stage outputs change in this release.

### Added

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
