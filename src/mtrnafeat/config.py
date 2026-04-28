"""Configuration dataclass and YAML loader. CLI overrides merge on top."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Config:
    # I/O
    data_dir: Path = Path("data")
    outdir: Path = Path("runs/default")
    db_files: dict[str, str] = field(
        default_factory=lambda: {"Human": "human_mt-mRNA_all.db", "Yeast": "yeast_mt-mRNA_all.db"}
    )
    alignment_file: str = "PAL2NL_aa-dna_alignment_yeast_human.txt"

    # Determinism
    seed: int = 42

    # Simulation defaults (landscape / significance)
    sim_seq_length: int = 150
    sim_num_sequences: int = 1500
    gradient_steps: int = 21
    gradient_seqs_per_step: int = 200
    sim_gc_conditions: dict[str, float] = field(
        default_factory=lambda: {
            "Sim Human (46% GC)": 0.46,
            "Sim Yeast CDS (30% GC)": 0.30,
            "Sim Yeast 5' UTR (7% GC)": 0.07,
        }
    )
    # Per-species (A, U, G, C) frequencies used to seed the simulated null
    # clouds in `analysis.landscape`. Empty dict → empirical frequencies are
    # computed from the .db files at runtime. Override in YAML to model e.g.
    # H-strand-only or L-strand-only nucleotide composition.
    sim_freqs_per_species: dict[str, dict[str, float]] = field(default_factory=dict)
    n_shuffles: int = 200

    n_workers: int = 4                       # multiprocessing pool size
    target_genes: tuple[str, ...] = (
        "COX1", "COX2", "COX3",
        "CYTB", "COB",
        "ATP86", "ATP9",
        "ND1", "ND2", "ND3", "ND4L4", "ND5", "ND6",
        "VAR1",
    )

    # `window` command: whole-transcript fold-and-compare with a rolling-mean
    # smooth on the per-position paired vector.
    rolling_window: int = 25
    # `local-probability` command: ViennaRNA RNAplfold defaults. window = local
    # window size; max_bp_span = hard distance cap; cutoff = min pair-prob to
    # report. Same defaults as RNAplfold's own ``-W 80 -L 50 -c 0.001``.
    rnaplfold_window: int = 80
    rnaplfold_max_bp_span: int = 50
    rnaplfold_cutoff: float = 0.001
    # window_nt / step_nt are kept for `significance` and `cofold` which still
    # do a sliding-window scan over each transcript.
    window_nt: int = 120
    step_nt: int = 10
    # Single canonical max_bp_span. The previous sweep across (50,100,150,300,inf)
    # produced too many figures for too little marginal insight; 300 nt is the
    # publication-relevant default. The cofold/ stage handles the soft-penalty
    # exploration separately.
    max_bp_span: int = 300
    window_size_sweep: tuple[int, ...] = (60, 120, 240)

    # CoFold soft long-range penalty (Proctor & Meyer 2013, NAR).
    # alpha is the penalty strength (kcal/mol), tau the decay constant (nt).
    # CoFold defaults: alpha=0.5, tau=640 (transcription ~50 nt/s).
    # Used by the dedicated `cofold` parameter-sweep stage. window/, tis/,
    # substitution/ all use plain Vienna; only `cofold` exercises CoFold itself.
    cofold_alpha: float = 0.5
    cofold_tau: float = 640.0
    cofold_alpha_sweep: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)
    cofold_tau_sweep: tuple[float, ...] = (160.0, 320.0, 640.0, 1280.0, 2560.0)
    # Per-window folder used by the `window` command. The .db ground-truth
    # dot-brackets were produced upstream with RNAstructure (Moran et al.,
    # `Fold -md 350`), so RNAstructure is the canonical match. ViennaRNA is
    # available via `--engine vienna` for the legacy path. Other commands
    # (stats, landscape, kinetic, …) still use ViennaRNA via core.thermo.
    fold_engine: str = "rnastructure"   # "rnastructure" (default) | "vienna"

    # Features
    max_loop_artifact_size: int = 50
    max_heatmap_size: int = 15

    # Substitution-thermo permutation test
    substitution_n_simulations: int = 200
    substitution_max_nt: int = 300   # codon-truncated to this max length

    # Plot
    dpi: int = 300
    plot_format: str = "svg"   # "svg" (editable) | "png" | "pdf"

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["data_dir"] = str(self.data_dir)
        d["outdir"] = str(self.outdir)
        d["target_genes"] = list(self.target_genes)
        d["window_size_sweep"] = list(self.window_size_sweep)
        d["cofold_alpha_sweep"] = list(self.cofold_alpha_sweep)
        d["cofold_tau_sweep"] = list(self.cofold_tau_sweep)
        return d


def load_config(path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> Config:
    """Load a Config from YAML; apply CLI overrides on top."""
    cfg = Config()
    if path is not None:
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}
        cfg = _apply(cfg, raw)
    if overrides:
        cfg = _apply(cfg, overrides)
    return cfg


def _apply(cfg: Config, raw: dict[str, Any]) -> Config:
    fields = {f for f in cfg.__dataclass_fields__}
    new_kwargs = cfg.to_dict()
    for k, v in raw.items():
        if k not in fields:
            raise ValueError(f"Unknown config key: {k}")
        new_kwargs[k] = v
    new_kwargs["data_dir"] = Path(new_kwargs["data_dir"])
    new_kwargs["outdir"] = Path(new_kwargs["outdir"])
    new_kwargs["target_genes"] = tuple(new_kwargs["target_genes"])
    new_kwargs["window_size_sweep"] = tuple(int(x) for x in new_kwargs["window_size_sweep"])
    new_kwargs["cofold_alpha_sweep"] = tuple(float(x) for x in new_kwargs["cofold_alpha_sweep"])
    new_kwargs["cofold_tau_sweep"] = tuple(float(x) for x in new_kwargs["cofold_tau_sweep"])
    return Config(**new_kwargs)
