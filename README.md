# mtrnafeat

Mitochondrial mRNA structural-feature analysis. A Python CLI that compares
DMS-MaPseq–derived RNA secondary structures against thermodynamic predictions
(Vienna MFE and CoFold), runs synonymous-recoding permutation tests against
the wild-type ΔG, sweeps CoFold parameters, scans transcripts with sliding
windows, zooms in on translation initiation sites, and emits a comparative
yeast↔human COX1 panel — all from a single config file.

The tool exists because the structural questions you ask of a
DMS-MaPseq-probed mitochondrial mRNA aren't all answered by one analysis.
You want a per-transcript summary, an MFE-vs-experimental landscape, an
element decomposition, a sliding-window picture, statistical significance
under shuffled nulls, a TIS zoom that handles missing 5'UTRs honestly, a
synonymous-recoding permutation test, a CoFold parameter sweep, and a
cross-species comparative — and you want them all reproducible from one
YAML config and one command.

## Features

The pipeline ships 12 stages (see [docs/STAGES.md](docs/STAGES.md) for the
deep dive on each):

- **stats** — per-transcript length, MFE, foldedness, GC%, paired-pair composition.
- **landscape** — simulated GC-gradient (empirical + symmetric-GC nulls) against the experimental DMS overlay.
- **features** — element decomposition, region-stratified heatmaps, phase-space contour.
- **window** — sliding-window DMS-vs-Vienna ΔG trace per gene.
- **significance** — dinucleotide-shuffle z-scores per gene (and optionally per window).
- **tis** — TIS −50/+50 nt zoom; honors 5'UTR when present, clamps at 5'-end when not.
- **substitution** — synonymous-recoding ΔG permutation test against flat-GC / positional-GC / synonymous null pools.
- **cofold** — CoFold (α, τ) parameter sweep against the DMS ΔG. CoFold is the co-transcriptional folding model of [Proctor & Meyer (2013, *NAR*)](https://academic.oup.com/nar/article/41/19/9090/2411166), which adds a soft penalty `f(d) = α · (1 − exp(−d/τ))` to every candidate base pair of sequence distance `d`. **α** (alpha, kcal/mol) is the asymptotic penalty strength — larger values discourage long-range pairs more strongly. **τ** (tau, nt) is the decay constant — the distance at which the penalty reaches `α·(1 − 1/e) ≈ 0.63·α`; small `τ` makes the penalty bite at short range, large `τ` lets short loops form freely and only penalizes truly long-range contacts. CoFold's published defaults (α = 0.5, τ = 640 nt) correspond to a transcription speed of ~50 nt/s with a ~12.8 s pairing window. The sweep grids both axes and picks the `(α, τ)` that best matches the experimental DMS ΔG per gene.
- **compare** — yeast↔human COX1 codon-aligned comparative (substitution table, ΔG track, directional flux).
- **gene-panel** — per-gene composition + paired-pair + local-foldedness panel.
- **kinetic** — DrTransformer kinetic folding (opt-in only).
- **plot** — re-render plots from cached CSVs.

## Requirements

- Python ≥ 3.11
- ViennaRNA Python bindings (`viennarna ≥ 2.5`)
- (Optional) DrTransformer for the `kinetic` stage

## Installation

```bash
git clone <repo-url>
cd mtrnafeat
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

## Quickstart

### 1. Smoke test (under 2 minutes, no real data)

Runs the full pipeline (minus `kinetic`) on the bundled `test_data/` fixtures:

```bash
bash examples/01_smoke_mini.sh
```

### 2. One stage on real data

```bash
mtrnafeat tis --config configs/all.yaml --outdir runs/tis
```

Outputs land in `runs/tis/tis/` (per-stage CSV + plot) and
`runs/tis/tables/tis_dms_vs_mfe.csv` (centralized copy).

### 3. Full pipeline (all stages, both species)

```bash
mtrnafeat run-all --parallel --config configs/all.yaml --outdir runs/all
```

Wall time: 30–60 min in parallel mode on a recent laptop.

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

## Stages at a glance

| Stage | What it does | Key outputs |
|-------|--------------|-------------|
| [stats](docs/STAGES.md#stats) | per-transcript summary | `stats/per_transcript_statistics.csv` |
| [landscape](docs/STAGES.md#landscape) | GC-gradient sim + experimental overlay | `landscape/landscape_overlay.svg`, `gc_gradient.csv` |
| [features](docs/STAGES.md#features) | element decomposition + heatmaps | `features/raw_motifs.csv`, `heatmap_size_ratios.svg` |
| [window](docs/STAGES.md#window) | sliding-window DMS-vs-Vienna | `window/window_{species}_{gene}.svg` |
| [significance](docs/STAGES.md#significance) | dinuc-shuffle z-scores | `significance/z_per_gene.csv` |
| [tis](docs/STAGES.md#tis) | −50/+50 nt TIS zoom | `tis/tis_dms_vs_mfe.csv`, `tis_zoom_grid.svg` |
| [substitution](docs/STAGES.md#substitution) | synonymous-recoding ΔG perm test | `substitution/substitution_thermo_summary.csv` |
| [cofold](docs/STAGES.md#cofold) | CoFold (α, τ) sweep | `cofold/cofold_best_per_gene.csv` |
| [compare](docs/STAGES.md#compare) | yeast↔human COX1 comparative | `compare/cox1_alignment_table.csv` |
| [gene-panel](docs/STAGES.md#gene-panel) | per-gene composition panel | `gene_panels/panel_{species}_{gene}.svg` |
| [kinetic](docs/STAGES.md#kinetic) | DrTransformer kinetic folding (opt-in) | `kinetic/kinetic_summary.csv` |
| [plot](docs/STAGES.md#plot) | re-render plots from cached CSVs | (stage-specific) |
| [run-all](docs/STAGES.md#run-all) | orchestrate the full pipeline | every output above |

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
├── features/            motifs, spans, heatmaps, phase-space contour
├── window/              per-gene sliding-window traces
├── significance/        z_per_gene.csv (and z_per_window.csv with --per-window)
├── tis/                 TIS −50/+50 zoom CSV + plot
├── substitution/        permutation distribution + summary + KDE panels
├── cofold/              (α, τ) sweep grid + heatmaps
├── compare/             COX1 alignment table, ΔG track, substitution heatmap
├── gene_panels/         one panel figure per (species, gene)
├── kinetic/             DrTransformer trajectories (only if you ran `kinetic`)
└── tables/              centralized copies of per-stage summary CSVs
```

Plots default to `.svg` (editable in Inkscape/Illustrator). Switch to
`png` or `pdf` by setting `plot_format` in the YAML.

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
- [CONTRIBUTING.md](CONTRIBUTING.md) — dev guide and conventions.
- [CHANGELOG.md](CHANGELOG.md) — version history.

## Citation

If you use mtrnafeat in academic work, please cite:

```
@software{mtrnafeat,
  author  = {Ahn, Ahram},
  title   = {mtrnafeat: Mitochondrial mRNA structural-feature analysis},
  year    = {2026},
  url     = {<repo-url>},
  version = {0.1.0}
}
```

## License

MIT — see [LICENSE](LICENSE).
