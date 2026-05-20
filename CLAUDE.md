# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**mtrnafeat** is a downstream analysis toolkit for mitochondrial mRNA structural models. It consumes DMS-MaPseq–derived dot-bracket structures and pairs them with thermodynamic predictions (ViennaRNA, RNAstructure, RNAplfold) to analyze structural features, compare against predictions, and perform synonymous-codon recoding tests.

## Setup & Install

Python ≥ 3.11 required.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .           # core (viennarna, numpy, pandas, scipy, matplotlib, seaborn, biopython, tqdm)
pip install -e '.[dev]'    # + pytest, ruff
pip install -e '.[kinetic]' # + DrTransformer (optional)
pip install -e '.[labels]'  # + adjustText (optional)
```

## Commands

```bash
# Run full test suite
pytest

# Run a single test file
pytest tests/test_stats.py -v

# Run a single test function
pytest tests/test_stats.py::test_name

# Lint / format
ruff check src tests
ruff format src tests

# CLI entry point
mtrnafeat <subcommand> --config <config.yaml> --outdir <outdir>

# Smoke test (under 2 min)
bash examples/01_smoke_mini.sh

# Full pipeline
mtrnafeat run-all --config configs/all.yaml --outdir runs/all
```

## Architecture

```
src/mtrnafeat/
├── cli.py               # Entry point; SUBCOMMANDS dict dispatches to commands/
├── config.py            # Config dataclass (~150 fields) + YAML loader
├── core/                # Foundational logic: thermo.py, structure.py, stats.py,
│                        #   cofold.py, shuffle.py, stacking.py, manifest.py
├── analysis/            # Per-stage analysis implementations (~14 modules)
├── commands/            # CLI subcommand implementations (~18 modules, one per stage)
├── engines/             # Folding engine wrappers (RNAstructure, ViennaRNA)
├── io/                  # Parsing/writing: db_parser.py, alignment.py, writers.py, etc.
└── viz/                 # Visualization (~14 plot modules) + style.py for theming
```

**Data flow:** YAML config → `.db` files (dot-bracket) → each command loads config → calls analysis module → writes CSVs + SVGs to `outdir/`.

**All subcommands** follow the signature: `def run(cfg: Config, args: list[str] | None = None) -> int:`

**`Config`** is the central parameter hub — all stages read from it; `config.py:load_config()` merges YAML + CLI overrides.

**`pipeline.py`** orchestrates independent stages in parallel or sequential order for `run-all`.

**CSV output** must go through `io.writers.canonical_csv()` (atomic writes, `%.10g` float formatting for byte-stable reproducibility). Never write CSVs directly.

**Plots** use `viz.style.plot_path()`, `apply_theme()`, and `style_axis()` for consistent publication-quality output.

## Key Conventions

- **Adding a config field**: update `config.py` dataclass, `configs/template.yaml`, and `docs/CONFIG.md`.
- **Adding a stage**: follow the 8-step process in `CONTRIBUTING.md` (analysis module → command module → cli.py → pipeline.py → docs + tests).
- **Tests with ViennaRNA**: mark with `@pytest.mark.needs_rna` (defined in `conftest.py`); they skip gracefully when ViennaRNA is unavailable.
- **Test fixtures**: `test_data/mini.config.yaml`, `test_data/mini_human.db`, `test_data/mini_yeast.db` (ND6 + ATP9 subset).
- Ruff ignores E501 (line length); target is 110 chars with Python 3.11 semantics.

## Key Files

| Path | Purpose |
|------|---------|
| `configs/template.yaml` | Annotated config template — copy to customize |
| `data/human_mt-mRNA_all.db` | Human mitochondrial .db structures |
| `data/yeast_mt-mRNA_all.db` | Yeast mitochondrial .db structures |
| `docs/STAGES.md` | Per-subcommand deep-dive (purpose, reads, writes, flags) |
| `docs/CONFIG.md` | Every config field documented |
| `docs/FIGURES.md` | Every output figure explained |
| `scripts/regenerate_fixtures.py` | Rebuild test fixture CSVs after algorithm changes |
| `scripts/parity_check.py` | Cross-run consistency verification |
