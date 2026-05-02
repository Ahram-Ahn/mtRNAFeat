"""Reusable input validation primitives.

Used by `mtrnafeat doctor` and `mtrnafeat validate-inputs` to surface
config / data problems before expensive analysis stages run. Every
checker returns a list of `ValidationIssue` records; nothing prints,
nothing exits. The CLI commands compose these checkers and decide how
to render the result.
"""
from __future__ import annotations

import os
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from mtrnafeat.config import Config, load_config
from mtrnafeat.io.annotations import annotation_df

Level = Literal["OK", "WARN", "ERROR", "OPTIONAL"]
_ALPHABET = set("ACGU")


@dataclass(frozen=True)
class ValidationIssue:
    level: Level
    item: str
    message: str

    def to_dict(self) -> dict[str, str]:
        return asdict(self)


def has_error(issues: Iterable[ValidationIssue]) -> bool:
    return any(i.level == "ERROR" for i in issues)


# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------

def check_config(path: str | Path | None) -> tuple[list[ValidationIssue], Config | None]:
    """Verify config file exists and loads. Returns issues and the loaded Config (or None)."""
    issues: list[ValidationIssue] = []
    if path is None:
        issues.append(ValidationIssue("ERROR", "config", "No --config path supplied"))
        return issues, None

    p = Path(path)
    if not p.exists():
        issues.append(ValidationIssue("ERROR", "config", f"File not found: {p}"))
        return issues, None

    try:
        cfg = load_config(p)
    except ValueError as exc:
        issues.append(ValidationIssue("ERROR", "config", f"{p}: {exc}"))
        return issues, None
    except Exception as exc:  # pragma: no cover — defensive; YAML parse errors etc.
        issues.append(ValidationIssue("ERROR", "config", f"{p}: {type(exc).__name__}: {exc}"))
        return issues, None

    issues.append(ValidationIssue("OK", "config", str(p)))
    return issues, cfg


def check_data_dir(cfg: Config) -> list[ValidationIssue]:
    p = Path(cfg.data_dir)
    if not p.exists():
        return [ValidationIssue("ERROR", "data_dir", f"Not found: {p}")]
    if not p.is_dir():
        return [ValidationIssue("ERROR", "data_dir", f"Not a directory: {p}")]
    return [ValidationIssue("OK", "data_dir", str(p))]


def check_target_genes(cfg: Config) -> list[ValidationIssue]:
    if not cfg.target_genes:
        return [ValidationIssue("WARN", "target_genes", "Empty target_genes list")]
    bad = [g for g in cfg.target_genes if not isinstance(g, str) or not g.strip()]
    if bad:
        return [ValidationIssue("ERROR", "target_genes", f"Invalid entries: {bad}")]
    return [ValidationIssue("OK", "target_genes", f"{len(cfg.target_genes)} genes")]


def check_fold_engine(cfg: Config) -> list[ValidationIssue]:
    allowed = {"vienna", "rnastructure"}
    if cfg.fold_engine not in allowed:
        return [
            ValidationIssue(
                "ERROR",
                "fold_engine",
                f"{cfg.fold_engine!r} not in {sorted(allowed)}",
            )
        ]
    return [ValidationIssue("OK", "fold_engine", cfg.fold_engine)]


# --------------------------------------------------------------------------
# .db files
# --------------------------------------------------------------------------

def check_db_file(path: str | Path) -> list[ValidationIssue]:
    """Validate a single .db file. Reuses parse_db (which checks length parity
    and bracket balance via DbRecord.__post_init__) and adds an alphabet check
    on top.
    """
    p = Path(path)
    label = p.name

    if not p.exists():
        return [ValidationIssue("ERROR", label, f"File not found: {p}")]

    # Imported lazily so this module stays import-light even when other I/O
    # helpers are unavailable (e.g. early in `doctor`).
    from mtrnafeat.io.db_parser import parse_db

    try:
        records = parse_db(p)
    except ValueError as exc:
        return [ValidationIssue("ERROR", label, str(exc))]
    except Exception as exc:  # pragma: no cover — defensive
        return [ValidationIssue("ERROR", label, f"{type(exc).__name__}: {exc}")]

    if not records:
        return [ValidationIssue("WARN", label, "No records parsed")]

    issues: list[ValidationIssue] = []
    for rec in records:
        bad = sorted({c for c in rec.sequence if c not in _ALPHABET})
        if bad:
            issues.append(
                ValidationIssue(
                    "ERROR",
                    f"{label}:{rec.raw_gene}",
                    f"Sequence contains non-ACGU characters: {bad}",
                )
            )

    if not has_error(issues):
        issues.insert(0, ValidationIssue("OK", label, f"{len(records)} records"))
    return issues


# --------------------------------------------------------------------------
# Annotations
# --------------------------------------------------------------------------

def check_annotations(cfg: Config) -> list[ValidationIssue]:
    """Cross-check db-file transcript lengths against the bundled annotation
    tables. UTR/CDS coordinates that don't fit within the transcript are
    errors; CDS lengths not divisible by 3 are warnings (some mt-genes have
    incomplete stop codons that get filled by polyadenylation).
    """
    from mtrnafeat.io.db_parser import parse_db

    issues: list[ValidationIssue] = []
    for species, fname in cfg.db_files.items():
        path = Path(cfg.data_dir) / fname
        if not path.exists():
            # already reported by check_db_file
            continue
        try:
            ann = annotation_df(species)
        except ValueError:
            issues.append(
                ValidationIssue(
                    "WARN",
                    f"annotation:{species}",
                    f"No bundled annotation for species {species!r}",
                )
            )
            continue
        try:
            records = parse_db(path)
        except Exception:  # pragma: no cover — already reported elsewhere
            continue
        ann_by_gene = {row["transcript"]: row for _, row in ann.iterrows()}
        for rec in records:
            row = ann_by_gene.get(rec.gene)
            if row is None:
                issues.append(
                    ValidationIssue(
                        "WARN",
                        f"annotation:{species}:{rec.gene}",
                        "No annotation row (UTR/CDS coords unavailable)",
                    )
                )
                continue
            l_tr = int(row["l_tr"])
            l_utr5 = int(row["l_utr5"])
            l_utr3 = int(row["l_utr3"])
            l_cds = int(row["l_cds"])
            if l_tr != len(rec.sequence):
                issues.append(
                    ValidationIssue(
                        "WARN",
                        f"annotation:{species}:{rec.gene}",
                        f"annotation l_tr={l_tr} but sequence length={len(rec.sequence)}",
                    )
                )
            if l_utr5 + l_cds + l_utr3 != l_tr:
                issues.append(
                    ValidationIssue(
                        "ERROR",
                        f"annotation:{species}:{rec.gene}",
                        f"l_utr5 + l_cds + l_utr3 ({l_utr5 + l_cds + l_utr3}) != l_tr ({l_tr})",
                    )
                )
            if l_cds > 0 and l_cds % 3 != 0:
                issues.append(
                    ValidationIssue(
                        "WARN",
                        f"annotation:{species}:{rec.gene}",
                        f"l_cds={l_cds} not divisible by 3 (mt-mRNAs sometimes have incomplete stops)",
                    )
                )
    if not issues:
        issues.append(ValidationIssue("OK", "annotations", "all transcripts match bundled annotation"))
    elif not has_error(issues):
        # Promote the first OK signal so output isn't confusing.
        n_warns = sum(1 for i in issues if i.level == "WARN")
        issues.insert(0, ValidationIssue("OK", "annotations", f"{n_warns} warning(s)"))
    return issues


# --------------------------------------------------------------------------
# Optional inputs
# --------------------------------------------------------------------------

def check_alignment(cfg: Config) -> list[ValidationIssue]:
    """The codon alignment used by `compare`. Optional — missing is a warning."""
    if not cfg.alignment_file:
        return [ValidationIssue("WARN", "alignment", "no alignment_file configured")]
    p = Path(cfg.data_dir) / cfg.alignment_file
    if not p.exists():
        return [
            ValidationIssue(
                "WARN",
                "alignment",
                f"file missing ({p}); compare stage will be skipped",
            )
        ]
    return [ValidationIssue("OK", "alignment", str(p))]


def check_modifications(cfg: Config) -> list[ValidationIssue]:
    """Modifications input is reserved for the Phase 4 `modifications` stage.
    For now: just report whether the canonical file exists at the expected
    location.
    """
    p = Path(cfg.data_dir) / "mt_modifications.tsv"
    if p.exists():
        return [ValidationIssue("OK", "modifications", str(p))]
    return [
        ValidationIssue(
            "WARN",
            "modifications",
            "no modification table at data/mt_modifications.tsv (modification-aware analysis disabled)",
        )
    ]


# --------------------------------------------------------------------------
# Aggregate
# --------------------------------------------------------------------------

def run_all_checks(config_path: str | Path) -> list[ValidationIssue]:
    """Composite check used by `validate-inputs`. Order matters: config first,
    then anything that depends on it.
    """
    issues, cfg = check_config(config_path)
    if cfg is None:
        return issues
    issues.extend(check_data_dir(cfg))
    issues.extend(check_target_genes(cfg))
    issues.extend(check_fold_engine(cfg))
    for fname in cfg.db_files.values():
        issues.extend(check_db_file(Path(cfg.data_dir) / fname))
    issues.extend(check_annotations(cfg))
    issues.extend(check_modifications(cfg))
    issues.extend(check_alignment(cfg))
    return issues


# --------------------------------------------------------------------------
# Rendering helpers (used by both commands)
# --------------------------------------------------------------------------

_LEVEL_WIDTH = 8


def format_table(rows: list[tuple[str, str, str]]) -> str:
    """Three-column fixed-width table: item | level | message."""
    if not rows:
        return ""
    item_w = max(len(r[0]) for r in rows)
    item_w = max(item_w, 24)
    out_lines = []
    for item, level, message in rows:
        out_lines.append(f"{item.ljust(item_w)}  {level.ljust(_LEVEL_WIDTH)}  {message}")
    return "\n".join(out_lines)


def issues_to_rows(issues: Iterable[ValidationIssue]) -> list[tuple[str, str, str]]:
    return [(i.item, i.level, i.message) for i in issues]


def issues_writable(outdir: str | Path) -> ValidationIssue:
    """Probe whether a directory is writable. Used by `doctor --outdir`."""
    p = Path(outdir)
    try:
        p.mkdir(parents=True, exist_ok=True)
        probe = p / ".mtrnafeat_writable_probe"
        probe.write_text("ok")
        probe.unlink()
        return ValidationIssue("OK", "outdir writable", str(p))
    except OSError as exc:
        return ValidationIssue("ERROR", "outdir writable", f"{p}: {exc}")


def serialize_issues(issues: Iterable[ValidationIssue]) -> list[dict[str, Any]]:
    return [i.to_dict() for i in issues]


def env_pathvar(name: str) -> str | None:
    return os.environ.get(name)
