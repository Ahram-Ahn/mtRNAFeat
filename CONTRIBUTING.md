# Contributing to mtrnafeat

Thanks for taking the time to contribute. This document covers the dev
loop, conventions, and how to extend the pipeline cleanly.

## Setup

Editable install with the `dev` extras pulls in pytest and ruff:

```bash
git clone <repo-url>
cd mtrnafeat
python -m venv .venv && source .venv/bin/activate
pip install -e '.[dev]'
```

ViennaRNA Python bindings are required for the analysis tests. If
`pip install viennarna` fails on your platform, fall back to bioconda:

```bash
conda install -c bioconda viennarna
```

For the opt-in kinetic stage, also install DrTransformer:

```bash
pip install -e '.[kinetic]'
```

## Tests

```bash
pytest                           # full suite (skips ViennaRNA-deps if missing)
pytest tests/test_io.py          # one file
pytest -k determinism            # by keyword
pytest --cov=mtrnafeat           # with coverage
```

Tests that require ViennaRNA are gated by the `needs_rna` marker
(see [tests/conftest.py](tests/conftest.py)) and skip gracefully when
the bindings aren't available — useful in CI environments without compiled
RNA libs.

The deterministic-output tests (`tests/test_determinism.py`) verify that
CSV outputs are byte-stable across runs given a fixed seed; if you change
analysis math, expect to update the recorded SHAs in that file.

## Lint

```bash
ruff check src tests             # lint
ruff format src tests            # format (in-place)
```

Ruff config lives in [pyproject.toml](pyproject.toml) — line length 110,
target Python 3.11, rules `E F W I B UP` with `E501` (line length) ignored.

## Smoke run

The shortest path to verify a fresh checkout:

```bash
bash examples/01_smoke_mini.sh
```

Runs the full pipeline (minus `kinetic`) on the `test_data/` fixtures.
Wall time: under 2 minutes.

## Adding a config field

A config field touches three places. Update them together to keep the
template, docs, and code in sync:

1. **Dataclass** — add the field to `Config` in
   [src/mtrnafeat/config.py](src/mtrnafeat/config.py) with a default value.
   Tuples need to be coerced from list in `_apply()` if you want YAML
   list → tuple conversion.
2. **Annotated YAML** — add the field at its default to
   [configs/template.yaml](configs/template.yaml) under the appropriate
   section divider, with a one-line comment.
3. **Docs** — add a section to [docs/CONFIG.md](docs/CONFIG.md) following
   the existing pattern (Type · Default · Controls · Used by · When to change).

Verify by loading the template:

```bash
python -c "from mtrnafeat.config import load_config; load_config('configs/template.yaml')"
```

## Adding a stage

A new stage `foo` requires:

1. **Command module** — `src/mtrnafeat/commands/foo.py` with
   `def run(cfg: Config, args: list[str] | None = None) -> int:`.
2. **CLI registration** — add `"foo": "mtrnafeat.commands.foo"` to
   `SUBCOMMANDS` in [src/mtrnafeat/cli.py](src/mtrnafeat/cli.py).
3. **Pipeline registration** (if it should run inside `run-all`) — add the
   module name to the `INDEPENDENT` tuple in
   [src/mtrnafeat/commands/pipeline.py](src/mtrnafeat/commands/pipeline.py).
   Stages that need external binaries (like `kinetic`) should NOT be added
   here — keep them opt-in.
4. **Output convention** — write CSVs via `canonical_csv()` and figures
   via `plot_path()` (both honor `cfg.plot_format` / `cfg.dpi`). For
   summary tables that belong in `tables/`, also call `tables_csv()`.
5. **Docs** — add a section to [docs/STAGES.md](docs/STAGES.md) with
   Purpose / Reads / Writes / Flags / Notes.
6. **Tests** — at minimum, add the stage to `tests/test_cli_smoke.py`
   so a fresh checkout exercises it.

## Commit style

Terse imperative subject, optional body. Prefix with stage name when the
change is scoped to one stage. Examples from `git log`:

```
TIS: honor 5UTR when present, clamp at 5-end when not
Initial commit: mtrnafeat — mitochondrial mRNA structural-feature analysis
```

Body (if present) explains *why*, not *what* — the diff already says what.

## Pull requests

- Run `pytest && ruff check src tests` before opening.
- One PR = one logical change. Split unrelated cleanups into their own PRs.
- If you change CSV output shape, update the determinism SHAs in
  [tests/test_determinism.py](tests/test_determinism.py) in the same PR
  and explain the shape change in the PR body.
