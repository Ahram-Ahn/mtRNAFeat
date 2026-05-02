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

- [`doctor`](#doctor) — environment diagnostics
- [`validate-inputs`](#validate-inputs) — pre-flight config + data checks
- [`stats`](#stats) — per-transcript statistics + boxplot
- [`landscape`](#landscape) — GC-gradient simulation + experimental overlay
- [`features`](#features) — element decomposition, heatmaps, phase-space, span ECDFs
- [`window`](#window) — whole-transcript fold-and-compare trace
- [`local-probability`](#local-probability) — RNAplfold per-position pair probabilities
- [`significance`](#significance) — dinucleotide-shuffle z-scores (+ optional cotrans scan)
- [`tis`](#tis) — −50/+50 nt zoom around start codon
- [`substitution`](#substitution) — synonymous-recoding ΔG permutation
- [`cofold`](#cofold) — CoFold (α, τ) parameter sweep
- [`compare`](#compare) — yeast↔human COX1 comparative
- [`gene-panel`](#gene-panel) — per-gene composition / paired-pair / foldedness
- [`kinetic`](#kinetic) — DrTransformer kinetic folding (opt-in)
- [`plot`](#plot) — re-render plots from cached CSVs
- [`run-all`](#run-all) — orchestrate the full pipeline

---

## doctor

**Purpose**. One-shot environment diagnostics: confirms Python ≥ 3.11,
the ViennaRNA Python binding, RNAplfold availability, RNAstructure
availability and `DATAPATH` (when `fold_engine` requires it), and
whether DrTransformer is on `PATH` for the optional `kinetic` stage.
Prints a fixed-width summary table.

**Reads**: nothing on disk other than the optional `--config` to surface
which config was loaded. The default `Config` is used otherwise.
**Writes**: stdout summary; optionally a JSON list of issues with
`--json PATH`.

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--json PATH` | Also write the issue list as JSON to `PATH`. |

**Exit status**: `0` if no `ERROR` rows; `1` otherwise. `WARN` and
`OPTIONAL` rows never fail the command.

---

## validate-inputs

**Purpose**. Pre-flight check of the supplied config and the data files
it references. Catches mismatched sequence/structure lengths, unbalanced
brackets, non-ACGU characters, mis-sized UTR/CDS annotations, and
missing optional inputs (alignment, modifications table) before an
expensive analysis stage runs.

**Reads**: `cfg`, `cfg.data_dir`, every `.db` file in `cfg.db_files`,
`cfg.alignment_file` (optional), `data/mt_modifications.tsv` (optional).
**Writes**: stdout summary; optionally a JSON list of issues with
`--json PATH`. No output directory is created.

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--json PATH` | Also write the issue list as JSON to `PATH`. |

**Exit status**: `0` if no `ERROR` rows; `1` otherwise. Missing optional
files (alignment, modifications) produce `WARN` and do not fail the
command.

**Note**. Requires `--config`. With no config, the command prints a
hint and exits non-zero.

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

**Purpose**. Whole-transcript **fold-and-compare** for every (species, gene).
For each transcript, folds the sequence three ways and emits per-position
and summary CSVs. The figure shows two of the three traces (DMS-derived
vs the engine-span fold); the unconstrained Vienna-full trace remains in
the CSVs but is not plotted, since for long mRNAs the configured
``max_bp_span`` already captures every realistic contact and the third
line was cluttering the plot.

1. **DMS-derived**: take the dot-bracket from the `.db` file and
   recompute its ΔG with ViennaRNA `eval_structure` (the `.db` header
   MFE is not used directly). **Plotted.**
2. **Vienna full**: pure ViennaRNA MFE on the whole transcript with no
   max-bp-span cap. **Computed and stored in the CSVs but not
   plotted.**
3. **Engine span**: refold the transcript with `max_bp_span = N`
   (default 300 nt) under the chosen engine — `rnastructure` (default,
   matches how the `.db` DMS-derived dot-brackets were produced
   upstream) or `vienna`. **Plotted.**

The plot stacks the two rolling paired-fraction traces (legend outside
the data axis, titled "DMS vs <engine>") and a transcript-architecture
strip (5'UTR / CDS / 3'UTR) below.

**Reads**: `cfg.max_bp_span`, `cfg.rolling_window`, `cfg.fold_engine`,
`cfg.db_files`, `cfg.target_genes`.
**Writes**:
- `window/window_per_position.csv` (one row per (species, gene, position) with the three rolling paired-fraction tracks)
- `window/window_summary.csv` (one row per (species, gene) with the three ΔGs and overall paired fractions)
- `window/window_failed.csv` (only if any genes failed — e.g. RNAstructure DATAPATH unset)
- `window/window_{SPECIES}_{GENE}.{svg|png}` (one per gene)

**Flags (after `--`)**:
| Flag | Type | Effect |
|------|------|--------|
| `--span N` | int | override `cfg.max_bp_span` (default 300). |
| `--rolling-window N` | int | rolling-mean width on the per-position paired vector (default `cfg.rolling_window` = 25). |
| `--engine NAME` | str | `rnastructure` (default) or `vienna`. |

---

## local-probability

**Purpose**. ViennaRNA's RNAplfold per-position pair-probability scan over
each transcript. For every (species, gene), reports the local pair
probability `p(i,j)` along the sequence — useful for spotting positions
that are nearly always paired (high probability across the local window)
versus positions that fluctuate or stay open. Complements the binary
"paired vs. unpaired" picture from the `window` stage.

**Reads**: `cfg.rnaplfold_window`, `cfg.rnaplfold_max_bp_span`,
`cfg.rnaplfold_cutoff`, `cfg.rolling_window` (smoothing for the plot
overlay), `cfg.db_files`, `cfg.target_genes`.
**Writes**:
- `local_probability/local_probability_per_position.csv` (one row per (species, gene, position) with raw and smoothed probability)
- `local_probability/local_probability_{SPECIES}_{GENE}.{svg|png}` (one per gene)

**Flags (after `--`)**:
| Flag | Type | Effect |
|------|------|--------|
| `--window N` | int | RNAplfold window size (default `cfg.rnaplfold_window` = 80). |
| `--max-bp-span N` | int | RNAplfold max base-pair span (default `cfg.rnaplfold_max_bp_span` = 50). |
| `--cutoff F` | float | RNAplfold pair-probability cutoff (default `cfg.rnaplfold_cutoff` = 0.001). |
| `--smooth N` | int | rolling-window width for the smoothed plot overlay (default `cfg.rolling_window`). |

These defaults match RNAplfold's own published defaults (`-W 80 -L 50 -c 0.001`).

---

## significance

**Purpose**. Per-gene structural-significance with two layers, both
restricted to `cfg.target_genes`:

1. **Sequence-level z-score** (always emitted, fast). Folds each whole
   transcript and a dinucleotide-shuffled null pool (Workman & Krogh
   1999 null), reports `Z_MFE = (observed − null_mean) / null_std`.
2. **Per-gene cotranscriptional / sliding scan** (with `--scan`).
   Walks the transcript in the chosen mode and tracks Vienna's MFE,
   ensemble diversity, and paired fraction. Z-scores the first-difference
   signals so structural rearrangement events stand out, and detects
   peaks above `--z-threshold`.

**Reads**: `cfg.n_shuffles`, `cfg.db_files`, `cfg.target_genes`. The
`--scan` pass uses its own `--window`/`--step` (not the `cfg.window_nt`
defaults).
**Writes**:
- `significance/z_per_gene.csv` (always)
- `significance/cotrans_per_window.csv` (only with `--scan`)
- `significance/cotrans_{SPECIES}_{GENE}.{svg|png}` (only with `--scan`, one per gene)

**Flags (after `--`)**:
| Flag | Type | Effect |
|------|------|--------|
| `--scan` | flag | enable the per-gene cotranscriptional scan. |
| `--mode prefix\|sliding` | str | scan mode (default `sliding`; `prefix` grows the prefix length 0..N to mimic nascent RNA). |
| `--window N` | int | sliding-mode window size (default 120). |
| `--step N` | int | prefix-mode growth step or sliding stride (default 30). |
| `--z-threshold T` | float | peak threshold for the smoothed Δ signals (default 2.0). |
| `--per-window` | flag | deprecated alias for `--scan` (preserved so `run-all` keeps working). |

`run-all` invokes `significance` with `--per-window`, so the cotrans
scan runs by default in the full pipeline.

---

## tis

**Purpose**. TIS zoom: a −50/+50 nt window centered on the annotated start
codon, comparing the DMS-derived ΔG (Vienna `eval_structure` on the
projected `.db` dot-bracket) to Vienna MFE on the same window.

**5'UTR handling**:
- When a 5'UTR is annotated, the upstream window is `min(L_5UTR, 50)` —
  the actual UTR length, capped at 50.
- When no 5'UTR is annotated, the window is clamped at the 5' end of the
  transcript (no synthetic upstream).
- Each row carries `L_5UTR_total`, `L_5UTR_in_window`, `L_CDS_in_window`,
  `Has_5UTR`, and `Has_Full_5UTR_Context` so downstream consumers can tell
  whether a row represents a full ±50 context or a truncated one.

**Plot conventions** (the figure is designed to be read on its own):
- Bars whose 5'UTR is shorter than 50 nt are **hatched** — upstream
  context is truncated and the value is not directly comparable to
  full-context bars.
- A `ΔG = 0` value is rendered as a small open-diamond marker plus the
  literal label `0` so it cannot be confused with missing data (the
  projected DMS structure has zero pairs in the window — a real
  measurement, not a NaN).
- A truly missing value (NaN) is annotated `n/a` instead of being
  drawn as a zero-height bar.

**Reads**: `cfg.db_files`, `cfg.target_genes`, plus per-gene 5'UTR
annotations from [`src/mtrnafeat/io/annotations.py`](../src/mtrnafeat/io/annotations.py).
**Writes**:
- `tis/tis_dms_vs_mfe.csv`
- `tis/tis_zoom_grid.{svg|png}`
- `tables/tis_dms_vs_mfe.csv` (centralized copy)

**Flags**: none.

**Code**: [src/mtrnafeat/analysis/tis.py](../src/mtrnafeat/analysis/tis.py),
[src/mtrnafeat/commands/tis.py](../src/mtrnafeat/commands/tis.py),
[src/mtrnafeat/viz/tis_plot.py](../src/mtrnafeat/viz/tis_plot.py).

---

## substitution

**Purpose**. The most novel piece of the package. For each (species, gene),
compares the wild-type CDS ΔG against five null pools generated by
codon-aware permutation:

1. **flat_gc** — uniform GC composition matching the gene's overall GC%
   (G and C are pooled, A and T are pooled).
2. **flat_acgu** — IID nucleotides preserving the gene's full A/C/G/T
   frequency vector independently. Reduces to `flat_gc` only when
   P(G)=P(C) and P(A)=P(T). Use when the transcript has marked G/C or
   A/U asymmetry — e.g. human ND6 (L-strand) has G ≈ 191 vs C ≈ 37, an
   asymmetry that overall-GC nulls hide entirely.
3. **positional_gc** — per-codon-position GC matched (1st / 2nd /
   wobble GC fractions preserved); G and C still pooled within each
   codon position.
4. **positional_acgu** — per-codon-position A/C/G/T frequencies
   preserved independently (12 parameters: 4 nucleotides × 3 codon
   positions). The strictest compositional null short of the
   synonymous pool; isolates structural signal from compositional
   asymmetry that varies by codon position.
5. **synonymous** — each codon replaced by a synonymous alternative,
   weighted by codon-usage frequency × positional GC.

Every pool member is folded with **plain ViennaRNA MFE** (the same
`thermo.fold_mfe` call as the wild-type reference), so the Z-score
compares apples to apples. CoFold is not used here — the dedicated
[`cofold`](#cofold) sweep handles the soft-penalty exploration.

**ΔG provenance — the `.db` header MFE is intentionally bypassed.** Both
wild-type references are recomputed by ViennaRNA on the `.db` records:

- `WT_MFE` — Vienna MFE of the wild-type sequence
  (`thermo.fold_mfe(seq)`); compared against the pool.
- `WT_DMS_Eval` — Vienna `eval_structure` on the `.db` dot-bracket
  truncated to the chunk length (`substitution_max_nt`).
- `WT_DMS_Eval_FullLength` — same but for the full transcript-length
  DMS structure (a sanity check that the chunk-truncated value tracks
  the full-length one).

The summary CSV carries a `DMS_dG_Source` column that spells this out.

**Reads**: `cfg.substitution_n_simulations`, `cfg.substitution_max_nt`,
`cfg.max_bp_span`, `cfg.db_files`, `cfg.target_genes`, `cfg.n_workers`.
**Writes**:
- `substitution/substitution_thermo_distribution.csv` (long-format raw — one row per (species, gene, pool, simulation) with that variant's ΔG)
- `substitution/substitution_kde_panels_{human,yeast}.{svg|png}` (per-species small-multiples KDE)
- `substitution/substitution_z_heatmap_{human,yeast}.{svg|png}` (per-species heatmap of Z(WT_MFE − Pool))
- `tables/substitution_thermo_summary.csv` (per-gene Z and p; the only summary CSV — there is no copy under `substitution/`)

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
`cfg.window_nt`, `cfg.step_nt`, `cfg.max_bp_span`, `cfg.db_files`,
`cfg.target_genes`, `cfg.n_workers`.
**Writes**:
- `cofold/cofold_grid.csv` (long-format: one row per (species, gene, α, τ) with the gap to DMS-eval ΔG)
- `cofold/cofold_best_per_gene.csv` (the (α, τ) cell that minimizes `|CoFold − DMS|` for each gene)
- `cofold/cofold_per_window_corr.csv` (unless `--no-window-corr`; per-window Pearson r between CoFold and DMS-projected ΔG)
- `cofold/cofold_gap_strip_{human,yeast}.{svg|png}` (per-species strip plot of |CoFold − DMS| across the (α, τ) sweep; one dot per cell, hue by τ, red ring on the best dot per gene, black square on the default α=0.5 / τ=640)
- `cofold/cofold_per_window_corr_{human,yeast}.{svg|png}` (unless `--no-window-corr`; one curve per τ, x-axis is α)
- `tables/cofold_best_per_gene.csv` (centralized copy)

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--no-window-corr` | skip the per-window correlation pass (faster; drops the per-window corr CSV and curves). |

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

**Purpose**. Per-(species, gene) 3-panel composition figure:

1. Base counts (A / U / G / C bar chart).
2. Paired-pair composition (G-C / A-U / G-U wobble as % of paired columns).
3. Rolling local paired-fraction trace along the transcript, with a
   dedicated 5'UTR / CDS / 3'UTR architecture strip below (same
   consistent label treatment as the `window` plot).

The panel title reports the overall sequence GC% so cross-gene
comparisons are easy. When 5'UTR/3'UTR annotations from
[`io/annotations.py`](../src/mtrnafeat/io/annotations.py) are missing for
a gene, the architecture strip is omitted (the trace fills the third
panel by itself).

**Reads**: `cfg.db_files`, `cfg.target_genes`, plus per-gene annotations
from `mtrnafeat.io.annotations`.
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

**Purpose**. Orchestrates every independent analysis stage in one
invocation. Stages have no cross-stage data dependencies, so
`--parallel` simply fires all of them concurrently as subprocesses
(pool size = `min(len(stages) + 1, cpu_count())`).

**Stages run** (in order, sequential mode; matches `INDEPENDENT` in
[src/mtrnafeat/commands/pipeline.py](../src/mtrnafeat/commands/pipeline.py)):

`stats`, `landscape`, `features`, `window`, `local_probability`,
`significance`, `tis`, `compare`, `substitution`, `cofold`, `gene_panel`.

`run-all` also passes `--per-window` to `significance` automatically, so
the cotrans scan runs by default.

**Stages NOT run**: `kinetic` (opt-in only — requires DrTransformer on
PATH), `plot` (re-render utility).

**Reads**: every config field consumed by the included stages.
**Writes**: every output of the included stages, under `cfg.outdir`.

**Flags (after `--`)**:
| Flag | Effect |
|------|--------|
| `--parallel` | fire stages concurrently in subprocesses (pool sized by `cpu_count()`). |
| `--skip a,b,c` | comma-separated stage names to skip (use the underscore form for `local_probability` and `gene_panel`, e.g. `--skip kinetic,local_probability`). |

**Examples**:
```bash
mtrnafeat run-all --config configs/all.yaml --outdir runs/all
mtrnafeat run-all --parallel --config configs/all.yaml --outdir runs/all
mtrnafeat run-all --parallel --skip significance,cofold --config configs/all.yaml --outdir runs/x

# Smoke run on the bundled fixtures (Vienna engine, --skip kinetic just
# matches what examples/01_smoke_mini.sh does):
mtrnafeat run-all --config test_data/mini.config.yaml --outdir runs/smoke -- --skip kinetic
```
