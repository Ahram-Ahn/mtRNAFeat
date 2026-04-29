# Stage reference

One section per `mtrnafeat` subcommand, in pipeline order. Every stage
takes `--config PATH`, `--outdir PATH`, and `--seed INT` as global flags
(documented once below). Subcommand-specific flags are passed after a `--`
separator on the CLI.

## Common flags (all subcommands)

| Flag | Type | Default | Meaning |
|------|------|---------|---------|
| `--config` | path | none | YAML config file (see [docs/CONFIG.md](CONFIG.md)). |
| `--outdir` | path | from config (`runs/default`) | Output root; per-stage subdirs created underneath. |
| `--seed` | int | from config (`42`) | Override the RNG seed for this run. |

Subcommand-specific args go **after** a literal `--`:

```bash
mtrnafeat substitution --config configs/all.yaml --outdir runs/x -- --n 1000 --max-nt 600
```

---

## Index

- [`stats`](#stats) — per-transcript statistics + boxplot
- [`landscape`](#landscape) — GC-gradient simulation + experimental overlay
- [`features`](#features) — element decomposition, heatmaps, phase-space
- [`window`](#window) — sliding-window DMS-vs-Vienna scan
- [`significance`](#significance) — dinucleotide-shuffle z-scores
- [`tis`](#tis) — −50/+50 nt zoom around start codon
- [`substitution`](#substitution) — synonymous-recoding ΔG permutation
- [`cofold`](#cofold) — CoFold (α, τ) parameter sweep
- [`compare`](#compare) — yeast↔human COX1 comparative
- [`gene-panel`](#gene-panel) — per-gene composition / paired-pair / foldedness
- [`kinetic`](#kinetic) — DrTransformer kinetic folding (opt-in)
- [`plot`](#plot) — re-render plots from cached CSVs
- [`run-all`](#run-all) — orchestrate the full pipeline

---

## stats

**Purpose**. Per-transcript summary statistics (length, MFE, foldedness,
GC%, paired-pair composition) for every (species, gene) in the `.db` files.

**Reads**: `cfg.db_files`, `cfg.data_dir`.
**Writes**:
- `stats/per_transcript_statistics.csv`
- `stats/stats_summary.{svg|png}` (per-species boxplot)
- `tables/per_transcript_statistics.csv` (centralized copy)

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--centroid` | also compute centroid-vs-MFE base-pair distance metrics. |

---

## landscape

**Purpose**. Compares DMS-experimental foldedness against simulated nulls
across a GC% gradient. Three families of nulls:

1. Symmetric-GC clouds at the conditions in `sim_gc_conditions`
   (default: 46% / 30% / 7% GC).
2. Empirical-frequency clouds derived from the `.db` files
   (controlled by `sim_freqs_per_species`; empty → derive from data).
3. A continuous gradient curve at every step in `gradient_steps`.

**Reads**: `cfg.sim_seq_length`, `cfg.sim_num_sequences`, `cfg.gradient_steps`,
`cfg.gradient_seqs_per_step`, `cfg.sim_gc_conditions`, `cfg.sim_freqs_per_species`,
plus all `db_files` for the experimental overlay.
**Writes**:
- `landscape/specific_conditions.csv`
- `landscape/gc_gradient.csv`
- `landscape/experimental_overlay.csv`
- `landscape/landscape_overlay.{svg|png}`
- `landscape/gc_gradient_curves.{svg|png}`
- `landscape/pairing_bias_{GC,AU,GU}.{svg|png}`

**Flags**: none.

---

## features

**Purpose**. Decomposes every (species, gene) structure into elements
(stems, hairpins, internal loops, multibranch loops) and motifs.
Stratifies by region (5'UTR, CDS, 3'UTR) and emits size-ratio heatmaps,
a phase-space contour, and span boxplots.

**Reads**: `cfg.db_files`, `cfg.max_loop_artifact_size`, `cfg.max_heatmap_size`,
plus simulated counterparts using `cfg.sim_*` for null comparison.
**Writes**:
- `features/raw_motifs.csv`
- `features/raw_spans.csv`
- `features/region_stratified_elements.csv`
- `features/heatmap_size_ratios.{svg|png}`
- `features/phase_space_contour.{svg|png}`
- `features/span_boxplot.{svg|png}`

**Flags**: none.

---

## window

**Purpose**. Slides a window of length `window_nt` along every transcript
at step `step_nt` and folds each window under `max_bp_span`. Produces a
DMS-vs-Vienna ΔG trace per gene.

**Reads**: `cfg.window_nt`, `cfg.step_nt`, `cfg.max_bp_span`, `cfg.db_files`,
`cfg.target_genes`.
**Writes**:
- `window/window_scan_metrics.csv`
- `window/window_scan_summary.csv`
- `window/window_{SPECIES}_{GENE}.{svg|png}` (one per gene)

**Flags (after `--`)**:
| Flag | Type | Effect |
|------|------|--------|
| `--window N` | int | override `window_nt` for this run. |
| `--step N` | int | override `step_nt`. |
| `--span N` | int | override `max_bp_span`. |

---

## significance

**Purpose**. Per-gene z-scores from dinucleotide-preserving shuffle nulls.
Optionally also per-window.

**Reads**: `cfg.n_shuffles`, `cfg.db_files`, `cfg.target_genes`,
plus `cfg.window_nt`/`cfg.step_nt`/`cfg.max_bp_span` if `--per-window`.
**Writes**:
- `significance/z_per_gene.csv`
- `significance/z_per_window.csv` (only with `--per-window`)

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--per-window` | also compute z-scores at every window position. |

---

## tis

**Purpose**. TIS zoom: a −50/+50 nt window centered on the annotated start
codon, comparing DMS-derived ΔG to Vienna MFE on the same window.

**5'UTR handling** (commit
[9362371](../README.md): "TIS: honor 5UTR when present, clamp at 5-end when not"):
- When a 5'UTR is annotated, the upstream window is `min(L_5UTR, 50)` —
  the actual UTR length, capped at 50.
- When no 5'UTR is annotated, the window is clamped at the 5' end of the
  transcript (no synthetic upstream).
- Each row carries `L_5UTR_total`, `L_5UTR_in_window`, `L_CDS_in_window`,
  `Has_5UTR`, and `Has_Full_5UTR_Context` so downstream consumers can tell
  whether a row represents a full ±50 context or a truncated one.
- Plots hatch the bars for genes whose upstream context is truncated.

**Reads**: `cfg.db_files`, `cfg.target_genes`, plus per-gene 5'UTR
annotations from `mtrnafeat.io.annotations`.
**Writes**:
- `tis/tis_dms_vs_mfe.csv`
- `tis/tis_zoom_grid.{svg|png}`
- `tables/tis_dms_vs_mfe.csv` (centralized copy)

**Flags**: none.

**Code**: [src/mtrnafeat/analysis/tis.py](../src/mtrnafeat/analysis/tis.py),
[src/mtrnafeat/commands/tis.py](../src/mtrnafeat/commands/tis.py).

---

## substitution

**Purpose**. The most novel piece of the package. For each (species, gene),
compares the wild-type CDS ΔG against three null pools generated by
synonymous recoding:

1. **Flat-GC** — uniform GC composition.
2. **Positional-GC** — per-codon-position GC matched.
3. **Synonymous** — each codon replaced by a synonymous alternative.

Every pool member is folded under CoFold; the wild-type ΔG is reported as
a z-score and empirical p-value against each pool.

**Reads**: `cfg.substitution_n_simulations`, `cfg.substitution_max_nt`,
`cfg.cofold_alpha`, `cfg.cofold_tau`, `cfg.db_files`, `cfg.target_genes`,
`cfg.n_workers`.
**Writes**:
- `substitution/substitution_thermo_distribution.csv` (long-format raw)
- `substitution/substitution_thermo_summary.csv` (per-gene Z, p)
- `substitution/kde_panels_*.{svg|png}` (small-multiples KDE)
- `substitution/z_heatmap.{svg|png}`
- `tables/substitution_thermo_summary.csv` (centralized copy)

**Flags (after `--`)**:
| Flag | Type | Effect |
|------|------|--------|
| `--n N` | int | override `substitution_n_simulations` (default 200). |
| `--max-nt N` | int | override `substitution_max_nt` (default 300). |

---

## cofold

**Purpose**. Two-axis parameter sweep over CoFold's `(alpha, tau)`. For
each (species, gene), folds the CDS at every cell of
`cofold_alpha_sweep × cofold_tau_sweep` and reports how close each cell
gets to the DMS-evaluated ΔG. Optionally also computes a per-window
correlation between CoFold ΔG and DMS ΔG.

**Background — what is CoFold, and what do α and τ mean?** CoFold
(Proctor & Meyer, 2013, *NAR*,
[doi:10.1093/nar/gkt600](https://academic.oup.com/nar/article/41/19/9090/2411166))
is a co-transcriptional folding model that augments Vienna's nearest-
neighbour energy with a soft penalty on every candidate base pair `(i, j)`
of sequence distance `d = j − i`:

```
f(d) = α · (1 − exp(−d / τ))      [kcal/mol]
```

The penalty rises monotonically with `d`, so long-range contacts pay more
than short-range ones — the kinetic intuition being that bases transcribed
far apart in time have less chance to find each other before downstream
sequence appears. The two parameters control the shape of that penalty:

- **α** (`alpha`, kcal/mol) — the **asymptotic penalty** as `d → ∞`. Larger
  α means stronger discrimination against long-range pairs. At α = 0 the
  penalty vanishes and CoFold collapses to plain Vienna MFE; at α = 1 a
  fully long-range pair pays an extra +1 kcal/mol.
- **τ** (`tau`, nt) — the **decay constant**. The penalty reaches
  `α · (1 − 1/e) ≈ 0.63·α` at distance `d = τ`. Small τ makes the penalty
  rise quickly (short loops still pay nearly the full α); large τ keeps
  short-range pairing free and only penalizes truly long-range contacts.

CoFold's published defaults are **α = 0.5 kcal/mol** and **τ = 640 nt**,
which the authors derive from a transcription speed of ~50 nt/s and a
~12.8 s "reach" for downstream pairing. The `cofold` stage tests whether
those defaults — or some other gene-specific setting — make Vienna's
folding match the experimental DMS structure better than plain MFE.

**Reads**: `cfg.cofold_alpha_sweep`, `cfg.cofold_tau_sweep`,
`cfg.window_nt`, `cfg.step_nt`, `cfg.db_files`, `cfg.target_genes`.
**Writes**:
- `cofold/cofold_grid.csv`
- `cofold/cofold_best_per_gene.csv`
- `cofold/cofold_per_window_corr.csv` (unless `--no-window-corr`)
- `cofold/gap_heatmap_*.{svg|png}` (one per gene)
- `cofold/per_window_corr_*.{svg|png}` (unless skipped)
- `tables/cofold_best_per_gene.csv` (centralized copy)

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--no-window-corr` | skip the per-window correlation pass (faster). |

---

## compare

**Purpose**. Yeast↔human COX1 comparative analysis on a codon-aligned
PAL2NL alignment: per-codon substitution table, transition/transversion
summary, directional substitution flux. (The per-codon local-ΔG track
was dropped — local-ΔG along the transcript is already covered by the
`window` stage at higher resolution and with both folding engines.)

**Reads**: `cfg.alignment_file`, `cfg.data_dir`, COX1 records from both
species's `.db` files.
**Writes**:
- `compare/cox1_alignment_table.csv`
- `compare/cox1_substitution_summary.csv`
- `compare/cox1_directional_flux.csv` (when computed)
- `compare/cox1_transition_transversion.csv` (when computed)
- `compare/cox1_substitution_heatmap.{svg|png}`
- `compare/cox1_directional_flux_heatmap.{svg|png}`
- `tables/cox1_*.csv` (centralized copies)

**Flags**: none.

**Notes**: stage skips silently when the alignment file or the COX1
record in either species is missing — keeps the pipeline working on test
fixtures that don't include COX1.

---

## gene-panel

**Purpose**. Generalized version of the legacy ND6-only panel — composition
+ paired-pair + local-foldedness — emitted once per (species, gene).

**Reads**: `cfg.db_files`, `cfg.target_genes`.
**Writes**:
- `gene_panels/panel_{SPECIES}_{GENE}.{svg|png}` (one figure per gene)

**Flags**: none.

---

## kinetic

**Purpose**. Co-transcriptional kinetic folding via DrTransformer. Opt-in
(never run by `run-all`); requires the `drtransformer` binary on PATH
(installed via `pip install -e '.[kinetic]'` or system-wide).

**Reads**: `cfg.db_files`. By default genes = `COX1, ND6, ATP9`.
**Writes**:
- `kinetic/kinetic_summary.csv`
- `kinetic/kinetic_trajectory.csv`
- `kinetic/kinetic_{SPECIES}_{GENE}.{svg|png}`

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--genes COX1,ND6,ATP9` | comma-separated gene list (default shown). |

**Notes**: returns exit code `2` if DrTransformer is not found.

---

## plot

**Purpose**. Re-render plots from a previous run's cached CSVs without
re-running the analysis.

**Usage**:
```bash
mtrnafeat plot landscape --from runs/all
mtrnafeat plot features  --from runs/all
```

**Args (positional)**:
| Arg | Values |
|-----|--------|
| `kind` | `landscape` \| `features` |

**Flags**:
| Flag | Effect |
|------|--------|
| `--from PATH` | source run directory (defaults to `cfg.outdir`). |

---

## run-all

**Purpose**. Orchestrates every independent stage in one invocation. The
pipeline has no cross-stage dependencies after the round-3 cleanup, so
`--parallel` simply fires all stages concurrently as subprocesses.

**Stages run** (in order, sequential mode):
`stats`, `landscape`, `features`, `window`, `significance`, `tis`,
`compare`, `substitution`, `cofold`, `gene_panel`.

**Stages NOT run**: `kinetic` (opt-in only), `plot` (re-render utility).

**Reads**: every config field consumed by the included stages.
**Writes**: every output of the included stages, under `cfg.outdir`.

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--parallel` | fire stages concurrently in subprocesses (pool sized by `cpu_count()`). |
| `--skip a,b,c` | comma-separated stage names to skip. |

**Examples**:
```bash
mtrnafeat run-all --config configs/all.yaml --outdir runs/all
mtrnafeat run-all --parallel --config configs/all.yaml --outdir runs/all
mtrnafeat run-all --parallel --skip significance,cofold --config configs/all.yaml --outdir runs/x
```
