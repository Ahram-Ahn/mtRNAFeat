# mtrnafeat

Mitochondrial mRNA structural-feature analysis. DMS-MaPseq vs Vienna MFE, codon-substitution thermodynamic permutation testing, CoFold parameter sweep against the experimental ΔG, sliding-window scans, and per-gene composition / paired-pair / local-foldedness panels.

## Install

```bash
pip install -e .                       # core install
pip install -e '.[kinetic]'            # + DrTransformer (opt-in only)
pip install -e '.[labels]'             # + adjustText for cleaner figure labels
pip install -e '.[dev]'                # + pytest, ruff
```

ViennaRNA Python bindings (`viennarna`) are required. If pip install fails, use bioconda: `conda install -c bioconda viennarna`.

## Layout

```
src/mtrnafeat/    package source
configs/          YAML configs (default + canonical all-species)
tests/            pytest suite — smoke + unit
examples/         shell scripts for common runs
```

## Run

The canonical single-run config covers all 14 unique mt-mRNAs across both Human and Yeast .db files in one pass:

```bash
mtrnafeat run-all --parallel --config configs/all.yaml --outdir runs/all
```

Single-stage invocations:

```
mtrnafeat stats          --config configs/all.yaml --outdir runs/x
mtrnafeat landscape      --config configs/all.yaml --outdir runs/x
mtrnafeat features       --config configs/all.yaml --outdir runs/x
mtrnafeat window         --config configs/all.yaml --outdir runs/x
mtrnafeat significance   --config configs/all.yaml --outdir runs/x   # dinuc-shuffle z
mtrnafeat tis            --config configs/all.yaml --outdir runs/x   # first-50-nt zoom
mtrnafeat substitution   --config configs/all.yaml --outdir runs/x   # synonymous-pool ΔG
mtrnafeat cofold         --config configs/all.yaml --outdir runs/x   # α/τ sweep vs DMS
mtrnafeat compare        --config configs/all.yaml --outdir runs/x   # COX1 yeast↔human
mtrnafeat gene-panel     --config configs/all.yaml --outdir runs/x   # per-gene panels
mtrnafeat kinetic        --config configs/all.yaml --outdir runs/x   # DrTransformer (opt-in)
mtrnafeat plot features  --from runs/x                                # regen plots from CSVs
```

## Development

```
pytest                  # smoke tests (skips ViennaRNA-dependent tests if unavailable)
ruff check src tests
```

## Output layout

Per-run output directory contains one subdirectory per stage plus a `tables/` directory that centralizes the per-gene summary CSVs (per-transcript stats, TIS zoom, substitution Z, CoFold best-fit, comparative substitution summary). Figures default to `.svg` for editability — set `plot_format: png` in the YAML to revert.
