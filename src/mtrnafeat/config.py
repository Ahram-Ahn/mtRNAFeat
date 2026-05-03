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
    # local-probability DMS-agreement outputs. Sliding window over the
    # per-position track to compute mean RNAplfold P(paired) vs DMS-derived
    # paired fraction; TIS summary uses circular-shift nulls (preserves
    # autocorrelation, no refolding).
    local_probability_scan_window_nt: int = 120
    local_probability_scan_step_nt: int = 10
    tis_upstream_nt: int = 50
    tis_downstream_nt: int = 50
    tis_n_circular_shifts: int = 1000
    # Sensitivity sweep over TIS context widths so a single fixed window
    # doesn't hide signal in long-5'UTR yeast transcripts (e.g. yeast COX1
    # has a 460-nt 5'UTR; the −50/+50 default is far too narrow). Each
    # (upstream, downstream) pair is summarized in
    # ``local_probability_TIS_sensitivity.csv`` next to the primary TIS
    # summary. Compute is cheap (slice + circular-shift means).
    tis_window_sweep_pairs: tuple[tuple[int, int], ...] = (
        (30, 30), (50, 50), (100, 100), (200, 200), (500, 500),
    )

    # ─────────── structure-deviation stage ───────────
    # Region-discovery pass that compares RNAplfold local pairing
    # probability (sequence-only thermodynamic prior) to the DMS-derived
    # paired fraction (experimental measurement) and calls regions where
    # the two diverge. Replaces the role of the old significance scan
    # with a biologically interpretable, region-centered output.
    structure_deviation_threshold: float = 0.25
    structure_deviation_min_region_length: int = 25
    structure_deviation_merge_gap: int = 10
    structure_deviation_high_threshold: float = 0.50
    structure_deviation_low_threshold: float = 0.30
    structure_deviation_tis_upstream: int = 30
    structure_deviation_tis_downstream: int = 60
    structure_deviation_stop_upstream: int = 60
    structure_deviation_stop_downstream: int = 30
    structure_deviation_early_cds_nt: int = 300
    # Mid-CDS bin = middle ``[lo, hi]`` fraction of the CDS, used by the
    # cross-gene heatmap so each gene contributes a comparable interior
    # window. Defaults span the central 40% of the CDS, leaving 30% at
    # each end for the early-CDS / late-CDS bins.
    structure_deviation_mid_cds_lo: float = 0.30
    structure_deviation_mid_cds_hi: float = 0.70
    # Late-CDS bin width (nt before the stop_upstream window). 300 nt is
    # roughly the largest mt-mRNA "interior" length that fits before the
    # stop-codon proximal window in human/yeast transcripts.
    structure_deviation_late_cds_window_nt: int = 300
    structure_deviation_top_labels: int = 5
    # Optional region-level null model. ``"none"`` (default) emits effect
    # sizes only; ``"dinuc"`` runs a dinucleotide-shuffle null and
    # max-statistic correction across called regions per gene.
    structure_deviation_null_model: str = "none"
    structure_deviation_n_null: int = 0
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
    # Per-window folder used by the `window` command. The .db DMS-derived
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

    def __post_init__(self) -> None:
        """Validate numeric ranges and enum choices.

        Failing here means a YAML config or CLI override has a value that
        would silently produce NaN tracks, flipped windows, or invalid
        upstream-tool calls. We'd rather hard-fail at load time than emit
        a junk run.
        """
        def _ge(name: str, value: float, lo: float) -> None:
            if not (value >= lo):
                raise ValueError(f"{name} must be >= {lo}, got {value!r}")

        def _gt(name: str, value: float, lo: float) -> None:
            if not (value > lo):
                raise ValueError(f"{name} must be > {lo}, got {value!r}")

        def _in_unit(name: str, value: float) -> None:
            if not (0.0 <= value <= 1.0):
                raise ValueError(f"{name} must be in [0, 1], got {value!r}")

        def _choice(name: str, value: Any, allowed: tuple[str, ...]) -> None:
            if value not in allowed:
                raise ValueError(
                    f"{name} must be one of {allowed}, got {value!r}"
                )

        # Determinism / parallelism
        _ge("seed", self.seed, 0)
        _gt("n_workers", self.n_workers, 0)

        # Simulation
        _gt("sim_seq_length", self.sim_seq_length, 0)
        _gt("sim_num_sequences", self.sim_num_sequences, 0)
        _gt("gradient_steps", self.gradient_steps, 0)
        _gt("gradient_seqs_per_step", self.gradient_seqs_per_step, 0)
        _gt("n_shuffles", self.n_shuffles, 0)

        # Smoothing / window scans
        _gt("rolling_window", self.rolling_window, 0)
        _gt("rnaplfold_window", self.rnaplfold_window, 0)
        _gt("rnaplfold_max_bp_span", self.rnaplfold_max_bp_span, 0)
        if self.rnaplfold_max_bp_span > self.rnaplfold_window:
            raise ValueError(
                f"rnaplfold_max_bp_span ({self.rnaplfold_max_bp_span}) must be "
                f"<= rnaplfold_window ({self.rnaplfold_window}); RNAplfold "
                "rejects L > W."
            )
        _in_unit("rnaplfold_cutoff", self.rnaplfold_cutoff)
        _gt("local_probability_scan_window_nt", self.local_probability_scan_window_nt, 0)
        _gt("local_probability_scan_step_nt", self.local_probability_scan_step_nt, 0)
        _ge("tis_upstream_nt", self.tis_upstream_nt, 0)
        _ge("tis_downstream_nt", self.tis_downstream_nt, 0)
        _gt("tis_n_circular_shifts", self.tis_n_circular_shifts, 0)

        # structure-deviation
        _in_unit("structure_deviation_threshold", self.structure_deviation_threshold)
        _in_unit("structure_deviation_high_threshold", self.structure_deviation_high_threshold)
        _in_unit("structure_deviation_low_threshold", self.structure_deviation_low_threshold)
        if self.structure_deviation_low_threshold > self.structure_deviation_high_threshold:
            raise ValueError(
                f"structure_deviation_low_threshold ({self.structure_deviation_low_threshold}) "
                f"must be <= structure_deviation_high_threshold "
                f"({self.structure_deviation_high_threshold})."
            )
        _gt("structure_deviation_min_region_length", self.structure_deviation_min_region_length, 0)
        _ge("structure_deviation_merge_gap", self.structure_deviation_merge_gap, 0)
        _ge("structure_deviation_tis_upstream", self.structure_deviation_tis_upstream, 0)
        _ge("structure_deviation_tis_downstream", self.structure_deviation_tis_downstream, 0)
        _ge("structure_deviation_stop_upstream", self.structure_deviation_stop_upstream, 0)
        _ge("structure_deviation_stop_downstream", self.structure_deviation_stop_downstream, 0)
        _gt("structure_deviation_early_cds_nt", self.structure_deviation_early_cds_nt, 0)
        _in_unit("structure_deviation_mid_cds_lo", self.structure_deviation_mid_cds_lo)
        _in_unit("structure_deviation_mid_cds_hi", self.structure_deviation_mid_cds_hi)
        if self.structure_deviation_mid_cds_lo > self.structure_deviation_mid_cds_hi:
            raise ValueError(
                f"structure_deviation_mid_cds_lo ({self.structure_deviation_mid_cds_lo}) "
                f"must be <= structure_deviation_mid_cds_hi "
                f"({self.structure_deviation_mid_cds_hi})."
            )
        _gt("structure_deviation_late_cds_window_nt",
            self.structure_deviation_late_cds_window_nt, 0)
        _ge("structure_deviation_top_labels", self.structure_deviation_top_labels, 0)
        _choice(
            "structure_deviation_null_model",
            self.structure_deviation_null_model,
            ("none", "dinuc"),
        )
        _ge("structure_deviation_n_null", self.structure_deviation_n_null, 0)
        if self.structure_deviation_null_model != "none" and self.structure_deviation_n_null == 0:
            raise ValueError(
                "structure_deviation_null_model is "
                f"{self.structure_deviation_null_model!r} but "
                "structure_deviation_n_null is 0; set n_null > 0 to draw "
                "any null samples."
            )

        # window / max_bp_span / cofold
        _gt("window_nt", self.window_nt, 0)
        _gt("step_nt", self.step_nt, 0)
        _gt("max_bp_span", self.max_bp_span, 0)
        _ge("cofold_alpha", self.cofold_alpha, 0.0)
        _gt("cofold_tau", self.cofold_tau, 0.0)

        # Engines
        _choice("fold_engine", self.fold_engine, ("rnastructure", "vienna"))

        # Substitution
        _gt("substitution_n_simulations", self.substitution_n_simulations, 0)
        _gt("substitution_max_nt", self.substitution_max_nt, 0)

        # Plot
        _gt("dpi", self.dpi, 0)
        _choice("plot_format", self.plot_format, ("svg", "png", "pdf"))

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["data_dir"] = str(self.data_dir)
        d["outdir"] = str(self.outdir)
        d["target_genes"] = list(self.target_genes)
        d["window_size_sweep"] = list(self.window_size_sweep)
        d["cofold_alpha_sweep"] = list(self.cofold_alpha_sweep)
        d["cofold_tau_sweep"] = list(self.cofold_tau_sweep)
        d["tis_window_sweep_pairs"] = [list(p) for p in self.tis_window_sweep_pairs]
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
    new_kwargs["tis_window_sweep_pairs"] = tuple(
        (int(p[0]), int(p[1])) for p in new_kwargs["tis_window_sweep_pairs"]
    )
    return Config(**new_kwargs)
