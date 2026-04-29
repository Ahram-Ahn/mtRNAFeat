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
- [Window scan](#window-scan) — `window_nt`, `step_nt`, `max_bp_span`,
  `window_size_sweep`
- [CoFold](#cofold) — `cofold_alpha`, `cofold_tau`, `cofold_alpha_sweep`,
  `cofold_tau_sweep`, `fold_engine`
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

## Window scan

### `window_nt`
- **Type**: int · **Default**: `120`
- **Controls**: sliding-window length (nt) for the window stage.
- **Used by**: [`window`](STAGES.md#window).

### `step_nt`
- **Type**: int · **Default**: `10`
- **Controls**: slide step between consecutive windows (nt).
- **Used by**: `window`.

### `max_bp_span`
- **Type**: int · **Default**: `300`
- **Controls**: maximum base-pair span passed to ViennaRNA fold (nt).
  300 nt is the publication-relevant default; larger values approach
  unrestricted MFE.
- **Used by**: `window`, `substitution`, anywhere structure is folded.

### `window_size_sweep`
- **Type**: tuple[int, ...] · **Default**: `(60, 120, 240)`
- **Controls**: candidate window sizes for ad-hoc sensitivity exploration.
- **Used by**: currently unused by the pipeline; retained for hand-rolled
  scripts.

---

## CoFold

The CoFold soft long-range penalty (Proctor & Meyer, *NAR* 2013,
[doi:10.1093/nar/gkt600](https://academic.oup.com/nar/article/41/19/9090/2411166))
augments the ViennaRNA energy model with `f(d) = alpha * (1 - exp(-d / tau))`,
in kcal/mol, applied to every candidate base pair of sequence-distance `d`.
The penalty is monotonically increasing in `d`, so long-range contacts
pay more than short-range ones — the kinetic intuition being that bases
transcribed far apart in time have less chance to find each other before
downstream sequence appears.

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
- **Used by**: `cofold` only.

### `cofold_tau_sweep`
- **Type**: tuple[float, ...] · **Default**:
  `(160.0, 320.0, 640.0, 1280.0, 2560.0)`
- **Controls**: tau grid for the [`cofold`](STAGES.md#cofold) sweep.
- **Used by**: `cofold` only.

### `fold_engine`
- **Type**: str · **Default**: `vienna`
- **Controls**: folding engine selector. Only `"vienna"` is currently
  wired through; the `cofold` stage always uses CoFold regardless.

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
- **Controls**: number of synonymous-recoding null pools generated per
  (species, gene). Each pool is folded under CoFold and compared to the
  wild-type ΔG.
- **Used by**: [`substitution`](STAGES.md#substitution).
- **When to change**: more iterations → tighter empirical p-values, slower
  run. CLI override: `mtrnafeat substitution -- --n 1000`.

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
