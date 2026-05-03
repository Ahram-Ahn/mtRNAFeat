# Configuration reference

Every field of the `Config` dataclass at [src/mtrnafeat/config.py](../src/mtrnafeat/config.py),
in dataclass order. Each section gives **type · default · what it controls ·
when to change**, plus the stage(s) that consume it.

The annotated YAML template at [configs/template.yaml](../configs/template.yaml)
matches this list 1:1. Unknown keys raise `ValueError` at load time.

CLI overrides: `--outdir` and `--seed` on any subcommand override the YAML
value for that run only.

## Index

- [I/O](#io) — `data_dir`, `outdir`, `db_files`, `alignment_file`
- [Determinism](#determinism) — `seed`
- [Simulation](#simulation) — `sim_seq_length`, `sim_num_sequences`,
  `gradient_steps`, `gradient_seqs_per_step`, `sim_gc_conditions`,
  `sim_freqs_per_species`, `n_shuffles`
- [Parallelism](#parallelism) — `n_workers`
- [Genes](#genes) — `target_genes`
- [Window & local-probability](#window--local-probability) — `rolling_window`,
  `rnaplfold_window`, `rnaplfold_max_bp_span`, `rnaplfold_cutoff`,
  `local_probability_scan_window_nt`, `local_probability_scan_step_nt`,
  `tis_upstream_nt`, `tis_downstream_nt`, `tis_n_circular_shifts`,
  `tis_window_sweep_pairs`,
  `structure_deviation_threshold`,
  `structure_deviation_min_region_length`,
  `structure_deviation_merge_gap`,
  `structure_deviation_high_threshold`,
  `structure_deviation_low_threshold`,
  `structure_deviation_tis_upstream`, `structure_deviation_tis_downstream`,
  `structure_deviation_stop_upstream`, `structure_deviation_stop_downstream`,
  `structure_deviation_early_cds_nt`, `structure_deviation_top_labels`,
  `structure_deviation_null_model`, `structure_deviation_n_null`,
  `window_nt`, `step_nt`, `max_bp_span`, `window_size_sweep`,
  `fold_engine`
- [CoFold](#cofold) — `cofold_alpha`, `cofold_tau`, `cofold_alpha_sweep`,
  `cofold_tau_sweep`
- [Features](#features) — `max_loop_artifact_size`, `max_heatmap_size`
- [Substitution-thermo](#substitution-thermo) — `substitution_n_simulations`,
  `substitution_max_nt`
- [Plot](#plot) — `dpi`, `plot_format`

---

## I/O

### `data_dir`
- **Type**: path · **Default**: `data`
- **Controls**: root directory for `db_files` and `alignment_file` lookup.
- **Used by**: every stage that reads sequence/structure data.
- **When to change**: if your `.db` files live outside the repo (e.g. on a
  shared dataset volume).

### `outdir`
- **Type**: path · **Default**: `runs/default`
- **Controls**: root output directory. Per-stage subdirectories
  (`stats/`, `landscape/`, …) and a centralized `tables/` are created underneath.
- **Used by**: every stage.
- **When to change**: every run. Override on the CLI with `--outdir runs/my-run`.

### `db_files`
- **Type**: dict[str, str] · **Default**:
  `{Human: human_mt-mRNA_all.db, Yeast: yeast_mt-mRNA_all.db}`
- **Controls**: species → `.db` filename. Path resolved relative to `data_dir`.
- **Used by**: every stage iterating over species.
- **When to change**: adding a new species, or pointing at an alternative
  reactivity dataset.

### `alignment_file`
- **Type**: str · **Default**: `PAL2NL_aa-dna_alignment_yeast_human.txt`
- **Controls**: codon-aligned PAL2NL alignment file used by the [`compare`](STAGES.md#compare) stage.
- **Used by**: `compare` only. Missing file → stage skips silently.
- **When to change**: replacing the alignment, or adding a new species pair.

---

## Determinism

### `seed`
- **Type**: int · **Default**: `42`
- **Controls**: every RNG used by the pipeline (numpy, dinucleotide shuffle,
  null-pool generation).
- **Used by**: `landscape`, `significance`, `substitution`, anywhere with
  randomness.
- **When to change**: rarely. Override on the CLI with `--seed N` for
  per-run reproducibility checks.

---

## Simulation

### `sim_seq_length`
- **Type**: int · **Default**: `150`
- **Controls**: length (nt) of each simulated random sequence in the
  landscape clouds.
- **Used by**: [`landscape`](STAGES.md#landscape).

### `sim_num_sequences`
- **Type**: int · **Default**: `1500`
- **Controls**: number of simulated sequences per condition for the
  landscape scatter clouds.
- **Used by**: `landscape`. Higher → tighter cloud, slower run.

### `gradient_steps`
- **Type**: int · **Default**: `21`
- **Controls**: number of GC% steps in the landscape gradient curve.
  21 → 0%, 5%, 10%, …, 100%.
- **Used by**: `landscape`.

### `gradient_seqs_per_step`
- **Type**: int · **Default**: `200`
- **Controls**: simulated sequences sampled per gradient step.
- **Used by**: `landscape`.

### `sim_gc_conditions`
- **Type**: dict[str, float] · **Default**: 3 conditions
  (`Sim Human (46% GC)`, `Sim Yeast CDS (30% GC)`, `Sim Yeast 5' UTR (7% GC)`)
- **Controls**: symmetric-GC null conditions for landscape overlay. Key is
  the legend label; value is the target GC fraction (0–1).
  A=U=(1−GC)/2, G=C=GC/2.
- **Used by**: `landscape`.
- **When to change**: adding a new species or testing a different GC null.

### `sim_freqs_per_species`
- **Type**: dict[str, dict[str, float]] · **Default**: `{}` (empty)
- **Controls**: per-species (A, U, G, C) frequencies for simulated nulls.
  Empty → empirical frequencies are derived from the `.db` files at runtime
  (default; honors actual H-strand enrichment in human).
- **Used by**: `landscape`.
- **When to change**: modeling H-strand-only / L-strand-only composition,
  or comparing against published whole-genome frequencies.

### `n_shuffles`
- **Type**: int · **Default**: `200`
- **Controls**: dinucleotide-shuffle iterations per gene for the
  significance stage z-scores.
- **Used by**: [`significance`](STAGES.md#significance).
- **When to change**: tighter z-score CIs need more shuffles (e.g. 1000)
  at the cost of runtime.

---

## Parallelism

### `n_workers`
- **Type**: int · **Default**: `4`
- **Controls**: process pool size for stages that fan out per-gene
  (window, substitution, cofold, kinetic).
- **Used by**: most heavy stages.
- **When to change**: tune to physical cores. Note: `run-all --parallel`
  uses a separate orchestrator pool sized by `cpu_count()`.

---

## Genes

### `target_genes`
- **Type**: tuple[str, ...] · **Default**: 14 mt-mRNAs
  (COX1, COX2, COX3, CYTB, COB, ATP86, ATP9, ND1, ND2, ND3, ND4L4, ND5, ND6, VAR1)
- **Controls**: gene whitelist. Pipeline iterates (species × `target_genes`);
  a pair is skipped if the gene is not in that species's `.db`.
- **Used by**: every stage that loops over genes.
- **When to change**: scoping a quick run to one gene
  (e.g. `target_genes: [COX1]`).

---

## Window & local-probability

The `window` and `local-probability` stages each produce a per-position
track per (species, gene). They share two ideas — a rolling-mean smoother
and a max-bp-span cap — but they answer different questions: `window`
re-folds the transcript and reports paired/unpaired booleans; `local-
probability` runs RNAplfold and reports a continuous pair probability.

### `rolling_window`
- **Type**: int · **Default**: `25`
- **Controls**: rolling-window width (nt) used by the `window` stage's
  paired-fraction smoother and by the `local-probability` plot's
  smoothed overlay. Centered window, 1-nt step.
- **Used by**: [`window`](STAGES.md#window),
  [`local-probability`](STAGES.md#local-probability).
- **CLI override**: `mtrnafeat window … -- --rolling-window N` /
  `mtrnafeat local-probability … -- --smooth N`.

### `rnaplfold_window`
- **Type**: int · **Default**: `80`
- **Controls**: ViennaRNA RNAplfold local-window size `W` (nt). Matches
  RNAplfold's published default `-W 80`.
- **Used by**: [`local-probability`](STAGES.md#local-probability).
- **CLI override**: `mtrnafeat local-probability … -- --window N`.

### `rnaplfold_max_bp_span`
- **Type**: int · **Default**: `50`
- **Controls**: RNAplfold max base-pair span `L` (nt). Matches
  `-L 50`.
- **Used by**: `local-probability`.
- **CLI override**: `mtrnafeat local-probability … -- --max-bp-span N`.

### `rnaplfold_cutoff`
- **Type**: float · **Default**: `0.001`
- **Controls**: RNAplfold pair-probability cutoff. Pairs with
  probability below this threshold are not reported. Matches `-c 0.001`.
- **Used by**: `local-probability`.
- **CLI override**: `mtrnafeat local-probability … -- --cutoff F`.

### `local_probability_scan_window_nt`
- **Type**: int · **Default**: `120`
- **Controls**: sliding-window size (nt) for the per-window
  RNAplfold-vs-DMS agreement table emitted by `local-probability`.
- **Used by**: [`local-probability`](STAGES.md#local-probability).
- **CLI override**: `mtrnafeat local-probability … -- --scan-window N`.

### `local_probability_scan_step_nt`
- **Type**: int · **Default**: `10`
- **Controls**: stride (nt) for the same sliding-window agreement table.
- **Used by**: `local-probability`.
- **CLI override**: `mtrnafeat local-probability … -- --scan-step N`.

### `tis_upstream_nt`
- **Type**: int · **Default**: `50`
- **Controls**: requested upstream (5') context (nt) around the start
  codon for the TIS summary. Truncated to whatever 5'UTR is actually
  annotated; the summary's `Has_Full_Upstream_Context` flag reports
  whether truncation occurred.
- **Used by**: `local-probability` (TIS summary).
- **CLI override**: `mtrnafeat local-probability … -- --tis-upstream N`.

### `tis_downstream_nt`
- **Type**: int · **Default**: `50`
- **Controls**: downstream (3') context (nt) around the start codon
  for the TIS summary.
- **Used by**: `local-probability` (TIS summary).
- **CLI override**: `mtrnafeat local-probability … -- --tis-downstream N`.

### `tis_window_sweep_pairs`
- **Type**: tuple[tuple[int, int], ...] · **Default**: `((30, 30), (50, 50), (100, 100), (200, 200), (500, 500))`
- **Controls**: a sensitivity sweep over TIS context widths. The same
  TIS-summary computation is run at every `(upstream, downstream)` pair
  here and emitted to `local_probability_TIS_sensitivity.csv` next to
  the primary `local_probability_TIS_summary.csv`. A single fixed
  `-50/+50` window hides signal in long-5'UTR yeast genes (Yeast COX1
  is unusually open at `-30/+30` but reads trivial at `-100/+100`) and
  is structurally over-wide for human mt-mRNAs whose 5'UTRs are 0-3 nt.
- **Used by**: `local-probability` (sensitivity table).

### `tis_n_circular_shifts`
- **Type**: int · **Default**: `1000`
- **Controls**: number of circular-shift draws used to compute the TIS
  empirical p-values (`Empirical_P_TIS_Low_*`). Circular shifts
  preserve the autocorrelation of the per-position track; the p-value
  is one-sided (`P(shifted_mean ≤ observed_mean)`), so smaller values
  indicate the TIS region is unusually open relative to its own
  transcript.
- **Used by**: `local-probability` (TIS summary).
- **CLI override**: `mtrnafeat local-probability … -- --tis-n-shuffles N`.

### `structure_deviation_threshold`
- **Type**: float · **Default**: `0.25`
- **Controls**: minimum smoothed `|P_model − P_DMS|` magnitude required
  to call a region in the `structure-deviation` stage.
- **Used by**: [`structure-deviation`](STAGES.md#structure-deviation).
- **CLI override**: `mtrnafeat structure-deviation … -- --threshold F`.

### `structure_deviation_min_region_length`
- **Type**: int · **Default**: `25`
- **Controls**: minimum region length in nt; intervals shorter than
  this are dropped after merging.
- **Used by**: `structure-deviation`.
- **CLI override**: `--min-region-length N`.

### `structure_deviation_merge_gap`
- **Type**: int · **Default**: `10`
- **Controls**: same-sign regions separated by ≤ this gap are merged
  before length filtering.
- **Used by**: `structure-deviation`.
- **CLI override**: `--merge-gap N`.

### `structure_deviation_high_threshold`
- **Type**: float · **Default**: `0.50`
- **Controls**: paired-fraction "high" cutoff used by the
  region-class assignment.
- **Used by**: `structure-deviation`.
- **CLI override**: `--high-threshold F`.

### `structure_deviation_low_threshold`
- **Type**: float · **Default**: `0.30`
- **Controls**: paired-fraction "low" cutoff used by the
  region-class assignment.
- **Used by**: `structure-deviation`.
- **CLI override**: `--low-threshold F`.

### `structure_deviation_tis_upstream` / `structure_deviation_tis_downstream`
- **Type**: int · **Defaults**: `30` / `60`
- **Controls**: TIS context window (nt) used to flag
  `Overlaps_TIS_Window` on each called region and to compute the
  TIS-aggregated metrics in `structure_deviation_gene_summary.csv`.
- **Used by**: `structure-deviation`.

### `structure_deviation_stop_upstream` / `structure_deviation_stop_downstream`
- **Type**: int · **Defaults**: `60` / `30`
- **Controls**: stop-codon context window (nt). Same role as the TIS
  window but at the 3′ end of the CDS.
- **Used by**: `structure-deviation`.

### `structure_deviation_early_cds_nt`
- **Type**: int · **Default**: `300`
- **Controls**: 3′ extent of the `early_CDS` heatmap bin (positions
  `tis_downstream` to `early_cds_nt` from the start codon).
- **Used by**: `structure-deviation`.

### `structure_deviation_top_labels`
- **Type**: int · **Default**: `5`
- **Controls**: reserved — formerly the count of top-|Δ| region labels
  drawn inline on the deviation panel; current per-gene plots use the
  bottom legend + colored region rectangles only.
- **Used by**: `structure-deviation`.
- **CLI override**: `--top-labels N`.

### `structure_deviation_null_model`
- **Type**: str · **Default**: `"none"` · **Values**: `"none"`, `"dinuc"` (planned)
- **Controls**: optional region-level null model. ``"none"`` (default)
  emits effect sizes only.
- **Used by**: `structure-deviation` (Phase 3 / planned).

### `structure_deviation_n_null`
- **Type**: int · **Default**: `0`
- **Controls**: number of null replicates when
  `structure_deviation_null_model != "none"`.
- **Used by**: `structure-deviation` (Phase 3 / planned).

### `window_nt`
- **Type**: int · **Default**: `120`
- **Controls**: sliding-window length (nt). The current `window`
  command does not slide (it folds the whole transcript), but the
  `cofold` per-window correlation pass and the `significance --scan`
  cotrans scan still use this field.
- **Used by**: [`cofold`](STAGES.md#cofold),
  [`significance --scan`](STAGES.md#significance).

### `step_nt`
- **Type**: int · **Default**: `10`
- **Controls**: slide step between consecutive windows (nt) for the
  same per-window passes that read `window_nt`.
- **Used by**: `cofold`, `significance --scan`.

### `max_bp_span`
- **Type**: int · **Default**: `300`
- **Controls**: maximum base-pair span passed to the folder (nt). 300 nt
  is the publication-relevant default; larger values approach
  unrestricted MFE.
- **Used by**: `window`, `substitution`, `cofold`, anywhere a structure
  is folded with a hard span cap.
- **CLI override**: `mtrnafeat window … -- --span N`.

### `window_size_sweep`
- **Type**: tuple[int, ...] · **Default**: `(60, 120, 240)`
- **Controls**: candidate window sizes for ad-hoc sensitivity
  exploration.
- **Used by**: currently unused by the pipeline; retained for
  hand-rolled scripts.

### `fold_engine`
- **Type**: str · **Default**: `rnastructure`
- **Controls**: folding engine for the `window` command. One of
  `rnastructure` or `vienna`.
  - `rnastructure` (default) — matches how the upstream `.db`
    DMS-derived dot-brackets were produced (`Fold -md <max_bp_span>`).
    Requires the `RNAstructure` binary on PATH and the `DATAPATH`
    environment variable pointing at its `data_tables/` directory.
  - `vienna` — `RNA.fold_compound(seq, md)` with `md.max_bp_span = N`.
    Pure-Python; no external binary needed.
- **Used by**: `window`. All other commands fold via ViennaRNA
  unconditionally (`window` is the only stage with a configurable
  engine).
- **CLI override**: `mtrnafeat window … -- --engine vienna|rnastructure`.

---

## CoFold

The CoFold-inspired soft long-range penalty (Proctor & Meyer, *NAR*
2013, Vol 41 Iss 9, e102,
[doi:10.1093/nar/gkt174](https://doi.org/10.1093/nar/gkt174))
augments the ViennaRNA energy model with `f(d) = alpha * (1 - exp(-d / tau))`,
in kcal/mol, applied to every candidate base pair of sequence-distance `d`.
The penalty is monotonically increasing in `d`, so long-range contacts
pay more than short-range ones — the kinetic intuition being that bases
transcribed far apart in time have less chance to find each other before
downstream sequence appears. (Note: this kcal/mol parameterization and
the mtrnafeat default `alpha = 0.5`, `tau = 640` are this implementation's
own choices; the published CoFold uses a dimensionless multiplicative
scaling with different numerical defaults — treat sweep outputs as
exploratory rather than as a literal CoFold reproduction.)

**The two parameters control the penalty's shape:**

- **α** (`alpha`, kcal/mol) is the *asymptotic* penalty as `d → ∞`. Larger
  α more strongly discourages long-range pairs; α = 0 collapses CoFold to
  plain Vienna MFE.
- **τ** (`tau`, nt) is the *decay constant*. The penalty reaches
  `α · (1 − 1/e) ≈ 0.63 · α` at distance `d = τ`. Small τ makes even
  short loops costly; large τ keeps short pairing free and only penalizes
  truly long-range contacts.

CoFold's published defaults (α = 0.5, τ = 640) correspond to a
transcription speed of ~50 nt/s and a ~12.8 s pairing window.

The [`cofold`](STAGES.md#cofold) parameter-sweep stage exercises the
(α, τ) grid directly; other stages (`window`, `tis`, `substitution`) use
plain Vienna unless explicitly invoked otherwise.

### `cofold_alpha`
- **Type**: float · **Default**: `0.5`
- **Controls**: penalty strength α (kcal/mol) — the asymptotic CoFold
  long-range penalty applied to every candidate base pair. Larger α →
  stronger discrimination against long-range pairs (α = 0 → plain
  Vienna).
- **When to change**: rarely directly; the [`cofold`](STAGES.md#cofold) sweep
  finds the best fit per gene.

### `cofold_tau`
- **Type**: float · **Default**: `640.0`
- **Controls**: decay constant τ (nt) — distance at which the CoFold
  penalty reaches `≈ 0.63 · α`. Roughly corresponds to a transcription
  speed of ~50 nt/s with a ~12.8 s pairing window.
- **When to change**: as above.

### `cofold_alpha_sweep`
- **Type**: tuple[float, ...] · **Default**: `(0.0, 0.25, 0.5, 0.75, 1.0)`
- **Controls**: alpha grid for the [`cofold`](STAGES.md#cofold) sweep.
  α = 0 collapses the cell to plain Vienna MFE — keeping it in the grid
  acts as a built-in control.
- **Used by**: `cofold` only.

### `cofold_tau_sweep`
- **Type**: tuple[float, ...] · **Default**:
  `(160.0, 320.0, 640.0, 1280.0, 2560.0)`
- **Controls**: tau grid for the [`cofold`](STAGES.md#cofold) sweep.
- **Used by**: `cofold` only.

---

## Features

### `max_loop_artifact_size`
- **Type**: int · **Default**: `50`
- **Controls**: loop sizes above this threshold (nt) are dropped from the
  element decomposition as folding artifacts.
- **Used by**: [`features`](STAGES.md#features).

### `max_heatmap_size`
- **Type**: int · **Default**: `15`
- **Controls**: axis cap (in element size) for the features-stage
  size-ratio heatmaps.
- **Used by**: `features`.

---

## Substitution-thermo

### `substitution_n_simulations`
- **Type**: int · **Default**: `200`
- **Controls**: number of synonymous-recoding null variants generated
  per pool per (species, gene). The stage uses three pools (flat-GC,
  positional-GC, synonymous), so total folds per gene are
  `3 × substitution_n_simulations + 2` (the `+2` are the two wild-type
  references). Every variant is folded under **plain ViennaRNA MFE**
  with `cfg.max_bp_span` — the wild-type comparison stays apples-to-
  apples.
- **Used by**: [`substitution`](STAGES.md#substitution).
- **When to change**: more iterations → tighter empirical p-values,
  slower run. CLI override: `mtrnafeat substitution -- --n 1000`.

### `substitution_max_nt`
- **Type**: int · **Default**: `300`
- **Controls**: codon-truncation cap (nt). Sequences longer than this are
  truncated to the first `substitution_max_nt` nt before recoding.
- **Used by**: `substitution`. CLI override: `-- --max-nt 600`.
- **When to change**: COX1 etc. are >1.5 kb; truncating keeps folding
  tractable. Increase if you have CPU budget.

---

## Plot

### `dpi`
- **Type**: int · **Default**: `300`
- **Controls**: figure DPI for raster output (PNG) and rasterized layers
  in mixed SVG.
- **Used by**: every plot.

### `plot_format`
- **Type**: str · **Default**: `svg`
- **Controls**: figure file format. `svg` is editable in
  Inkscape/Illustrator; `png` is convenient for quick viewing; `pdf` is
  good for LaTeX.
- **Used by**: every plot.
