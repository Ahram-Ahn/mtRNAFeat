# mtrnafeat

`mtrnafeat` is a downstream analysis toolkit for mitochondrial mRNA
structural models. It consumes DMS-MaPseq–derived dot-bracket structures
and pairs them with thermodynamic predictions (ViennaRNA, RNAstructure,
RNAplfold) to help answer:

1. Which mt-mRNA regions are locally paired or accessible?
2. Do DMS-derived structures agree with thermodynamic predictions?
3. Are conclusions robust to folding-engine choice?
4. Are known RNA modifications likely to confound DMS interpretation?
5. Are coding sequences unusually structured relative to codon-aware
   synonymous nulls?
6. Do structural features differ by species, strand, gene, or region?

It is **not** a raw FASTQ-to-reactivity pipeline, a replacement for
ShapeMapper2 / RNA Framework / SEISMIC-RNA, a new thermodynamic folding
algorithm, or a validated co-transcriptional folding simulator for long
mRNAs. See [What mtrnafeat is not](#what-mtrnafeat-is-not) for the
boundary.

## Features

The pipeline ships **11 analysis stages** plus orchestration (`run-all`),
a re-render utility (`plot`), and two pre-flight helpers (`doctor`,
`validate-inputs`). See [docs/STAGES.md](docs/STAGES.md) for the per-stage
deep dive.

- **doctor** — environment diagnostics (Python, ViennaRNA, RNAplfold, RNAstructure, DrTransformer).
- **validate-inputs** — pre-flight check of config + `.db` files + annotations + optional alignment / modifications table.
- **stats** — per-transcript length, MFE, foldedness, GC%, paired-pair composition.
- **landscape** — simulated GC-gradient (empirical + symmetric-GC nulls) against the experimental DMS overlay.
- **features** — element decomposition, region-stratified heatmaps, phase-space contour, base-pairing-distance ECDFs.
- **window** — whole-transcript fold-and-compare trace per gene (DMS-derived vs configured-span engine fold; `--engine vienna|rnastructure`).
- **local-probability** — ViennaRNA RNAplfold per-position pair-probability track per gene, overlaid against the DMS-derived dot-bracket (smoothed paired-fraction track + signed Δ); emits per-window agreement and a TIS-vs-CDS-background summary with circular-shift empirical p-values.
- **structure-deviation** — region-discovery pass on top of the RNAplfold-vs-DMS comparison; calls intervals where `P_model − P_DMS` exceeds a threshold and classifies each region as `model_high_dms_low` (DMS-open / thermodynamically foldable), `model_low_dms_high` (DMS-protected / nonlocally paired), `concordant_paired`, `concordant_open`, `mixed_deviation`, or `ambiguous`. Emits per-gene four-panel plots, a per-species lollipop summary, and a cross-gene architectural-bin heatmap.
- **significance** — *(legacy; not run by default)* per-gene dinucleotide-shuffle z-score with empirical p-value; optional within-gene structural-change scan flagging candidate transition peaks (within-gene Z, not a statistical test). Biological interpretation has moved to `structure-deviation`; invoke `significance` directly only when you want the thermodynamic-null QC.
- **tis** — TIS −50/+50 nt zoom; honors 5'UTR when present, clamps at the 5'-end when not.
- **substitution** — synonymous-recoding ΔG permutation test (flat-GC / positional-GC / synonymous null pools, all folded under plain Vienna MFE so the wild-type vs pool comparison is apples-to-apples). The DMS reference ΔG is recomputed by Vienna `eval_structure` on the `.db` dot-bracket — the `.db` header MFE is intentionally bypassed.
- **cofold** — CoFold (α, τ) parameter sweep against the DMS ΔG. CoFold is the co-transcriptional folding model of [Proctor & Meyer (2013, *NAR*)](https://academic.oup.com/nar/article/41/19/9090/2411166), which adds a soft penalty `f(d) = α · (1 − exp(−d/τ))` to every candidate base pair of sequence distance `d`. **α** (alpha, kcal/mol) is the asymptotic penalty strength — larger values discourage long-range pairs more strongly. **τ** (tau, nt) is the decay constant — the distance at which the penalty reaches `α·(1 − 1/e) ≈ 0.63·α`; small `τ` makes the penalty bite at short range, large `τ` lets short loops form freely and only penalizes truly long-range contacts. CoFold's published defaults (α = 0.5, τ = 640 nt) correspond to a transcription speed of ~50 nt/s with a ~12.8 s pairing window. The sweep grids both axes and picks the `(α, τ)` that best matches the experimental DMS ΔG per gene.
- **compare** — yeast↔human COX1 codon-aligned comparative (substitution heatmap, transition/transversion summary, directional substitution flux).
- **gene-panel** — per-gene composition + paired-pair + local-foldedness panel with a 5'UTR / CDS / 3'UTR architecture strip.
- **kinetic** — DrTransformer kinetic folding (opt-in; never auto-runs).
- **plot** — re-render `landscape` or `features` plots from cached CSVs.
- **run-all** — orchestrate every independent stage; supports `--parallel` and `--skip stage1,stage2`.

## Requirements

- Python ≥ 3.11
- ViennaRNA Python bindings (`viennarna ≥ 2.5`) — required for every
  stage that folds RNA (i.e. all of them)
- (Optional) RNAstructure with `DATAPATH` set — needed if you keep
  `fold_engine: rnastructure` (the default) for the `window` stage. See
  [Folding engines](#folding-engines) below.
- (Optional) DrTransformer on `PATH` for the `kinetic` stage (opt-in;
  never auto-runs)

## Installation

```bash
git clone https://github.com/Ahram-Ahn/mtRNAFeat.git
cd mtRNAFeat
python -m venv .venv && source .venv/bin/activate
pip install -e .                         # core
```

Optional extras:

```bash
pip install -e '.[kinetic]'              # + DrTransformer for `mtrnafeat kinetic`
pip install -e '.[labels]'               # + adjustText for cleaner figure labels
pip install -e '.[dev]'                  # + pytest, ruff
```

If `pip install viennarna` fails on your platform, install via bioconda:

```bash
conda install -c bioconda viennarna
```

### Folding engines

Two engines are wired through. The `window` command can use either; every
other stage (significance, cofold, kinetic, …) folds via ViennaRNA.

| Engine | Selector (`fold_engine` YAML / `--engine` CLI) | Requires | Notes |
|--------|------------------------------------------------|----------|-------|
| RNAstructure | `rnastructure` (**default**) | the `RNAstructure` binary on PATH and the `DATAPATH` env var pointing at its `data_tables/` directory | matches how the upstream `.db` DMS-derived dot-brackets were produced (`Fold -md 350`) |
| ViennaRNA | `vienna` | the `viennarna` Python package (already required) | `RNA.fold_compound` with `md.max_bp_span` |

If you don't have RNAstructure installed, override per run with
`mtrnafeat window … -- --engine vienna`, or set `fold_engine: vienna` in
the YAML. The default `configs/all.yaml` uses RNAstructure; the smoke
fixture `test_data/mini.config.yaml` uses Vienna so the test suite
doesn't need RNAstructure.

## Quickstart

### 0. Sanity-check the environment and inputs

```bash
mtrnafeat doctor                                          # Python, ViennaRNA, optional binaries
mtrnafeat validate-inputs --config configs/all.yaml       # config + .db + annotations
```

`doctor` exits non-zero only on hard failures (missing ViennaRNA, missing
RNAstructure when `fold_engine: rnastructure`). Optional tools like
DrTransformer are reported but never block. `validate-inputs` is similar:
missing optional files (alignment, modifications table) are warnings.

### 1. Smoke test (under 2 minutes, no real data)

Runs every analysis stage (minus `kinetic`) on the bundled `test_data/`
fixtures (Human ND6 + Yeast ATP9, ~1.5 kb total):

```bash
bash examples/01_smoke_mini.sh
```

Or the parallel variant:

```bash
bash examples/02_smoke_parallel.sh
```

### 2. One stage on real data

```bash
mtrnafeat tis --config configs/all.yaml --outdir runs/tis
```

Outputs land in `runs/tis/tis/` (per-stage CSV + plot) and
`runs/tis/tables/tis_dms_vs_mfe.csv` (centralized copy).

Other common single-stage invocations:

```bash
# Whole-transcript fold-and-compare with the Vienna engine override:
mtrnafeat window --config configs/all.yaml --outdir runs/win -- --engine vienna

# RNAplfold per-position pair probability with custom window/cutoff:
mtrnafeat local-probability --config configs/all.yaml --outdir runs/lp -- --window 80 --cutoff 0.001

# Substitution permutation with 1000 nulls per pool, 600-nt chunk:
mtrnafeat substitution --config configs/all.yaml --outdir runs/sub -- --n 1000 --max-nt 600

# Significance + the optional cotranscriptional sliding-window scan:
mtrnafeat significance --config configs/all.yaml --outdir runs/sig -- --scan --mode sliding
```

### 3. Full pipeline (all stages, both species)

```bash
mtrnafeat run-all --parallel --config configs/all.yaml --outdir runs/all
```

Wall time: 30–60 min in parallel mode on a recent laptop.

### 4. Scoping to one gene

Override `target_genes` in the YAML to scope a quick exploratory run:

```yaml
# my-cox1-only.yaml
data_dir: data
target_genes: [COX1]
fold_engine: vienna
```

```bash
mtrnafeat run-all --config my-cox1-only.yaml --outdir runs/cox1-only
```

## Configuration

Every analysis parameter lives in one YAML file. The repo ships three:

| File | Purpose |
|------|---------|
| [configs/default.yaml](configs/default.yaml) | Minimal silent defaults — what `run-all` uses with no `--config`. |
| [configs/all.yaml](configs/all.yaml) | Canonical real run; targets every mt-mRNA in both species. |
| [configs/template.yaml](configs/template.yaml) | **Annotated, all-fields template — copy this and edit.** |

To customize, copy `configs/template.yaml`, edit, and pass via `--config`:

```bash
cp configs/template.yaml configs/my-run.yaml
# ... edit my-run.yaml ...
mtrnafeat run-all --config configs/my-run.yaml --outdir runs/my-run
```

Per-field documentation is in [docs/CONFIG.md](docs/CONFIG.md). Two
parameters can also be overridden on the CLI without editing YAML:

```bash
mtrnafeat <subcommand> --config configs/all.yaml --outdir runs/x --seed 7
```

Subcommand-specific flags go after a literal `--`:

```bash
mtrnafeat substitution --config configs/all.yaml --outdir runs/x -- --n 1000 --max-nt 600
```

## Recommended workflow

1. Generate DMS reactivity / DMS-constrained structures upstream
   (e.g. ShapeMapper2 → RNAstructure `Fold -d`). `mtrnafeat` consumes
   the resulting dot-brackets; it does not generate them.
2. `mtrnafeat doctor` to check the environment.
3. `mtrnafeat validate-inputs --config <your.yaml>` to catch input
   mistakes before expensive analyses.
4. `mtrnafeat run-all --config <your.yaml>` for the standard pipeline.
   Inspect `local_probability/`, `window/`, and `gene_panels/` first —
   these are the most directly interpretable outputs for long mRNAs.
5. Run extended modules (`substitution`, `cofold`) if their assumptions
   suit your question; treat their outputs as exploratory until
   validated.
6. Treat `kinetic` (DrTransformer) and `cofold` parameter sweeps as
   exploratory: they are parameter-sensitive and not validated for long
   mature mt-mRNAs.

A note on terminology: `.db` dot-brackets are **DMS-derived /
DMS-constrained models**, not direct experimental observations. Use the
language "DMS-derived dot-bracket" or "DMS-constrained model" rather
than "ground truth" — secondary-structure prediction under DMS
constraints is still a model.

## Stages at a glance

| Stage | What it does | Key outputs |
|-------|--------------|-------------|
| [doctor](docs/STAGES.md#doctor) | environment diagnostics | (stdout summary; optional `--json`) |
| [validate-inputs](docs/STAGES.md#validate-inputs) | pre-flight config + data validation | (stdout summary; optional `--json`) |
| [stats](docs/STAGES.md#stats) | per-transcript summary | `stats/per_transcript_statistics.csv`, `stats_summary.svg` |
| [landscape](docs/STAGES.md#landscape) | GC-gradient sim + experimental overlay | `landscape/landscape_overlay.svg`, `gc_gradient.csv`, `pairing_bias_{GC,AU,GU}.svg` |
| [features](docs/STAGES.md#features) | element decomposition + heatmaps + ECDFs | `features/raw_motifs.csv`, `heatmap_size_ratios.svg`, `phase_space_contour.svg`, `span_boxplot.svg` |
| [window](docs/STAGES.md#window) | whole-transcript fold-and-compare trace | `window/window_{species}_{gene}.svg`, `window_per_position.csv`, `window_summary.csv` |
| [local-probability](docs/STAGES.md#local-probability) | RNAplfold per-position pair probabilities + DMS overlay | `local_probability/local_probability_per_position.csv`, `local_probability_per_window.csv`, `local_probability_TIS_summary.csv`, `local_probability_TIS_sensitivity.csv`, per-gene `.svg` |
| [structure-deviation](docs/STAGES.md#structure-deviation) | RNAplfold-vs-DMS deviation regions (called, classified, summarized) | `structure_deviation/structure_deviation_per_position.csv`, `structure_deviation_regions.csv`, `structure_deviation_gene_summary.csv`, `structure_deviation_gene_region_matrix.csv`, per-gene SVG, per-species lollipop SVG, cross-gene heatmap SVG |
| [significance](docs/STAGES.md#significance) | dinuc-shuffle z-scores (+ optional cotrans scan) | `significance/z_per_gene.csv`, `cotrans_per_window.csv` (with `--scan`) |
| [tis](docs/STAGES.md#tis) | −50/+50 nt TIS zoom | `tis/tis_dms_vs_mfe.csv`, `tis_zoom_grid.svg` |
| [substitution](docs/STAGES.md#substitution) | synonymous-recoding ΔG perm test (Vienna MFE) | `substitution/substitution_thermo_distribution.csv`, `tables/substitution_thermo_summary.csv` |
| [cofold](docs/STAGES.md#cofold) | CoFold (α, τ) sweep | `cofold/cofold_best_per_gene.csv`, `cofold_grid.csv`, per-species |gap| strip plot |
| [compare](docs/STAGES.md#compare) | yeast↔human COX1 comparative | `compare/cox1_alignment_table.csv`, `cox1_substitution_heatmap.svg`, `cox1_directional_flux_heatmap.svg` |
| [gene-panel](docs/STAGES.md#gene-panel) | per-gene composition + architecture | `gene_panels/panel_{species}_{gene}.svg` |
| [kinetic](docs/STAGES.md#kinetic) | DrTransformer kinetic folding (opt-in) | `kinetic/kinetic_summary.csv`, per-gene trajectory plots |
| [plot](docs/STAGES.md#plot) | re-render `landscape` / `features` plots | (stage-specific, no recomputation) |
| [run-all](docs/STAGES.md#run-all) | orchestrate every analysis stage | every output above except `kinetic` |

## Inputs

The pipeline reads three kinds of files from `cfg.data_dir` (default `data/`):

**`.db` files** (per species). Three-line records:

```
>COX1: -150.4 kcal/mol
GUAGCUAUCAGCAUC...
((((....))))....
```

Header line, RNA sequence, dot-bracket structure. One record per gene.
Default filenames: `human_mt-mRNA_all.db`, `yeast_mt-mRNA_all.db` —
override via `db_files` in the config.

**Codon-aligned alignment** (PAL2NL format) for the `compare` stage.
Default: `PAL2NL_aa-dna_alignment_yeast_human.txt`. Optional — `compare`
skips silently if missing.

**The YAML config itself** (`--config`).

## Outputs

Each run writes under `cfg.outdir` (default `runs/default`):

```
runs/<your-run>/
├── stats/               per-transcript statistics CSV + boxplot
├── landscape/           gradient curves, overlay, pairing-bias plots
├── features/            motifs, spans, region table, heatmaps, phase-space contour, span ECDFs
├── window/              per-gene whole-transcript trace + per-position + summary CSVs
├── local_probability/   per-gene RNAplfold pair-probability tracks + DMS overlay (per-position, per-window, TIS summary CSVs)
├── structure_deviation/ per-position + region + gene-summary + gene-region-matrix CSVs; per-gene 4-panel SVG; per-species lollipop SVG; cross-gene heatmap SVG
├── significance/        z_per_gene.csv (+ cotrans_per_window.csv and per-gene plots when --scan)
├── tis/                 TIS −50/+50 zoom CSV + plot
├── substitution/        permutation distribution CSV + per-species KDE panels + per-species z heatmap
├── cofold/              (α, τ) sweep grid + per-species |gap| strip plot + per-window correlation curves
├── compare/             COX1 alignment table, substitution heatmap, directional-flux heatmap
├── gene_panels/         one panel figure per (species, gene)
├── kinetic/             DrTransformer summary + per-gene trajectory plots (only if you ran `kinetic`)
└── tables/              centralized copies of per-stage summary CSVs
                         (per_transcript_statistics, tis_dms_vs_mfe,
                          substitution_thermo_summary, cofold_best_per_gene,
                          cox1_substitution_summary, cox1_directional_flux,
                          cox1_transition_transversion)
```

Plots default to `.svg` (editable in Inkscape/Illustrator). Switch to
`png` or `pdf` by setting `plot_format` in the YAML.

## What mtrnafeat is not

- Not a raw FASTQ-to-reactivity pipeline. Use ShapeMapper2
  ([Busan & Weeks, 2018](https://doi.org/10.1261/rna.061945.117)) or
  RNA Framework ([Incarnato et al., 2018](https://doi.org/10.1093/nar/gky486))
  upstream.
- Not a replacement for ShapeMapper2, RNA Framework, or SEISMIC-RNA /
  DREEM-style ensemble deconvolution. `mtrnafeat` consumes their outputs.
- Not a new thermodynamic folding algorithm. ViennaRNA
  ([Lorenz et al., 2011](https://doi.org/10.1186/1748-7188-6-26))
  and RNAstructure do the folding; `mtrnafeat` summarizes and compares.
- Not a validated co-transcriptional folding simulator for long mt-mRNAs.
  `cofold` and `kinetic` are exploratory parameter sweeps, not predictions.
- Not a tertiary-structure prediction tool. Pseudoknots and tertiary
  contacts are outside standard dot-bracket secondary-structure models.

## Interpretation cautions

1. DMS-derived dot-bracket structures are model-derived, not direct
   observations.
2. DMS mainly informs A/C accessibility/reactivity; interpretation
   differs across bases, chemistries, and modifications.
3. m1A can alter DMS signal; mask or annotate known m1A sites if you
   have a modification table. (Modification-aware analysis is a planned
   stage; the bundled `data/mt_modifications.tsv` is not yet wired into
   the pipeline.)
4. Full-length MFE is a coarse descriptor for long mRNAs; prefer
   `local-probability` for region-level interpretation.
5. Synonymous-null models (`substitution`) test whether observed coding
   sequences have unusual structural properties relative to codon-aware
   nulls. They do **not**, on their own, prove selection on RNA
   structure.
6. Yeast and human mitochondrial transcript architectures differ
   markedly (UTR lengths, GC content, strand origin). Aggregate
   yeast-vs-human comparisons should report strand and region rather
   than pooling silently.
7. Human ND6 is L-strand encoded and has very different composition
   from the H-strand transcripts; do not pool it without annotation.
8. CoFold / kinetic outputs are exploratory unless validated
   experimentally for the transcript and condition of interest.
9. The `.db` header `kcal/mol` value is the upstream caller's reported
   energy under its constrained fold and is not authoritative for
   thermodynamic comparisons; `mtrnafeat` recomputes ΔG via Vienna
   `eval_structure` where it matters.

## Examples

[examples/](examples/) contains shell wrappers for the common run patterns:

| Script | Purpose |
|--------|---------|
| [01_smoke_mini.sh](examples/01_smoke_mini.sh) | Sequential smoke run on `test_data/` fixtures (<2 min). |
| [02_smoke_parallel.sh](examples/02_smoke_parallel.sh) | Same smoke run, with `--parallel`. |
| [03_real_all.sh](examples/03_real_all.sh) | Full pipeline on real data, both species, parallel. |
| [04_substitution_run.sh](examples/04_substitution_run.sh) | `substitution` stage only, with tunable `N` and `MAX_NT`. |
| [05_single_step.sh](examples/05_single_step.sh) | Run one subcommand by name (`STEP=tis ./05_single_step.sh`). |

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for dev setup, tests, lint, and
conventions for adding config fields or stages.

```bash
pip install -e '.[dev]'
pytest                                   # smoke + unit
ruff check src tests                     # lint
```

## Documentation

- [docs/STAGES.md](docs/STAGES.md) — every subcommand: purpose, I/O, flags, notes.
- [docs/CONFIG.md](docs/CONFIG.md) — every YAML config field explained.
- [docs/FIGURES.md](docs/FIGURES.md) — every figure output: what is plotted, how to read each axis, source CSV.
- [CONTRIBUTING.md](CONTRIBUTING.md) — dev guide and conventions.
- [CHANGELOG.md](CHANGELOG.md) — version history.

## Citation

If you use mtrnafeat in academic work, please cite:

```
@software{mtrnafeat,
  author  = {Ahn, Ahram},
  title   = {mtrnafeat: Mitochondrial mRNA structural-feature analysis},
  year    = {2026},
  url     = {https://github.com/Ahram-Ahn/mtRNAFeat},
  version = {0.3.0}
}
```

A machine-readable [`CITATION.cff`](CITATION.cff) is also provided.

## License

MIT — see [LICENSE](LICENSE).
