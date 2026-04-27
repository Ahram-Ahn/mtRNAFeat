"""True kinetic folding via DrTransformer subprocess.

This is the experiment that decides whether the user's "nascent ≈ DMS"
hypothesis survives — DrTransformer is a deterministic heuristic
cotranscriptional folder. We run it, parse the trajectory, and compare the
end-state pair-table to (a) the DMS structure and (b) the global Vienna MFE.

Reference: Badelt et al. 2023, Bioinformatics 39(1):btad034.
https://academic.oup.com/bioinformatics/article/39/1/btad034/6992659

DrTransformer is an optional dependency (`pip install drtransformer`).
This module degrades to a clear ImportError if it isn't available.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mtrnafeat.config import Config
from mtrnafeat.constants import KIND_TRUE_KINETIC, canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.io.db_parser import parse_db

KIND = KIND_TRUE_KINETIC


@dataclass
class KineticResult:
    species: str
    gene: str
    sequence: str
    end_structures: list[tuple[float, str]]   # (occupancy, dot-bracket) at end
    trajectory: pd.DataFrame                  # parsed time-courses


def has_drtransformer() -> bool:
    return shutil.which("DrTransformer") is not None


def run_drtransformer(seq: str, outdir: Path, transcription_rate: float = 30.0,
                      logml_t: float = -2.0, name: str = "mtrnafeat") -> tuple[Path, Path]:
    """Run DrTransformer on `seq`. Returns (drf_path, log_path)."""
    if not has_drtransformer():
        raise RuntimeError(
            "DrTransformer CLI not found on PATH. Install with `pip install drtransformer`."
        )
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "DrTransformer",
        "--name", name,
        "--t-end", "60",
        "--t-lin", "10",
        "--t-log", "5",
        "--logfile",
    ]
    proc = subprocess.run(cmd, input=f">{name}\n{seq}\n", text=True,
                           capture_output=True, cwd=str(outdir))
    if proc.returncode != 0:
        raise RuntimeError(f"DrTransformer failed (exit {proc.returncode}): {proc.stderr[:500]}")
    drf = outdir / f"{name}.drf"
    log = outdir / f"{name}.log"
    if not drf.exists():
        raise RuntimeError(f"Expected DrTransformer output {drf} not found")
    return drf, log


def parse_drf(drf_path: Path) -> pd.DataFrame:
    """Parse DrTransformer .drf trajectory file → long-form DataFrame.

    Format (one struct per row):
        id  transcript_len  occupancy  energy  structure
    """
    rows = []
    with open(drf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                rows.append({
                    "Struct_ID": parts[0],
                    "Transcript_Length": int(parts[1]),
                    "Occupancy": float(parts[2]),
                    "Energy": float(parts[3]),
                    "Structure": parts[4],
                })
            except ValueError:
                continue
    return pd.DataFrame(rows)


def end_state_population(traj: pd.DataFrame) -> list[tuple[float, str]]:
    """Return (occupancy, structure) pairs at the longest transcript length,
    sorted by descending occupancy."""
    if traj.empty:
        return []
    L = traj["Transcript_Length"].max()
    final = traj[traj["Transcript_Length"] == L].sort_values("Occupancy", ascending=False)
    return [(float(o), s) for o, s in zip(final["Occupancy"], final["Structure"])]


def compare_kinetic_to_dms(end_pop: list[tuple[float, str]], dms_struct: str, mfe_struct: str) -> dict:
    """Most-populated kinetic end-state vs DMS vs MFE — bp_distance and Hamming."""
    if not end_pop:
        return {}
    occ, kin = end_pop[0]
    return {
        "Top_Occupancy": occ,
        "Kinetic_End_Structure": kin,
        "BPdist_Kinetic_vs_DMS": float(thermo.bp_distance(kin, dms_struct)),
        "BPdist_Kinetic_vs_MFE": float(thermo.bp_distance(kin, mfe_struct)),
        "BPdist_DMS_vs_MFE": float(thermo.bp_distance(dms_struct, mfe_struct)),
    }


def run_kinetic_for_genes(cfg: Config, genes: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not has_drtransformer():
        raise RuntimeError("DrTransformer not available; install with `pip install drtransformer`.")
    if genes is None:
        genes = ["COX1", "ND6", "ATP9"]
    summaries = []
    trajectories = []
    for species, fname in cfg.db_files.items():
        path = cfg.data_dir / fname
        recs = {r.gene: r for r in parse_db(path)}
        for gene in genes:
            target = canonical_gene(gene)
            if target not in recs:
                continue
            rec = recs[target]
            with tempfile.TemporaryDirectory() as tmp:
                drf, _log = run_drtransformer(rec.sequence, Path(tmp), name=f"{species}_{target}")
                traj = parse_drf(drf)
            traj["Species"] = species
            traj["Gene"] = target
            trajectories.append(traj)
            mfe_struct, mfe = thermo.fold_mfe(rec.sequence)
            cmp = compare_kinetic_to_dms(end_state_population(traj), rec.structure, mfe_struct)
            summaries.append({"Species": species, "Gene": target, "Length": len(rec.sequence),
                              "MFE": mfe, **cmp})
    summary_df = pd.DataFrame(summaries)
    traj_df = pd.concat(trajectories, ignore_index=True) if trajectories else pd.DataFrame()
    return summary_df, traj_df
