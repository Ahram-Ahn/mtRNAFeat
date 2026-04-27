# Changelog

All notable changes to this project will be documented in this file. Format
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioned
per [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
  (flat-GC, positional-GC, synonymous null pools) under CoFold.
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
