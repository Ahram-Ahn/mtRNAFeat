# Example run scripts

Run these from the repository root, or execute them directly from
`examples/`. Most scripts accept environment variables so you can reuse
them without editing the file.

| Script | Situation | Typical command |
|---|---|---|
| `01_smoke_mini.sh` | Fast bundled test data run | `examples/01_smoke_mini.sh` |
| `02_smoke_parallel.sh` | Same smoke run, parallel orchestration | `examples/02_smoke_parallel.sh` |
| `03_real_all.sh` | Full built-in Human + Yeast analysis | `examples/03_real_all.sh` |
| `04_substitution_run.sh` | Only synonymous-codon substitution thermodynamics | `examples/04_substitution_run.sh` |
| `05_single_step.sh` | One subcommand only | `STEP=features CONFIG=configs/all.yaml examples/05_single_step.sh` |
| `06_preflight_only.sh` | Check environment and inputs before a long run | `CONFIG=configs/all.yaml examples/06_preflight_only.sh` |
| `07_real_all_with_comparison.sh` | Full Human + Yeast run, explicitly including comparison stages | `examples/07_real_all_with_comparison.sh` |
| `08_multisample_independent.sh` | Independent analysis of many `.db` samples | `CONFIG=configs/my-experiment.yaml examples/08_multisample_independent.sh` |
| `09_fast_overview.sh` | Quick overview: stats, landscape, features, TIS, gene panels | `examples/09_fast_overview.sh` |
| `10_landscape_features_only.sh` | Re-run only changed landscape/features analyses | `examples/10_landscape_features_only.sh` |
| `11_one_gene_debug.sh` | One-gene debug run from an existing config | `GENES=COX1 examples/11_one_gene_debug.sh` |
| `12_db_folder_to_multisample.sh` | Build a config from every `.db` file in one folder, then run | `DB_DIR=data/my_samples examples/12_db_folder_to_multisample.sh` |

Common environment variables:

- `CONFIG`: YAML config to use. Default depends on the script.
- `OUTDIR`: output directory. Defaults to a `runs/...` path.
- `GENES`: comma-separated gene list for `11_one_gene_debug.sh`.
- `DB_DIR`: folder containing `.db` files for
  `12_db_folder_to_multisample.sh`.
- `RUN=0`: for `12_db_folder_to_multisample.sh`, generate the config but
  do not start the analysis.

Use the same virtual environment setup described in the README first.
The scripts run `python -m mtrnafeat.cli` with `PYTHONPATH=src`, using the
`python` currently on your `PATH`. Activate your environment first, or set
`PYTHON=/path/to/python` when calling a script.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
```
