"""CoFold parameter sweep: which (alpha, tau) makes Vienna MFE track DMS ΔG?

For each (species, gene), we fold the wild-type CDS at every (alpha, tau)
combination on the configured sweep grid (`cfg.cofold_alpha_sweep`
× `cfg.cofold_tau_sweep`). For each combination we compare to:
  - DMS-evaluated ΔG (the energy of the experimental structure under the
    standard Vienna model).

Outputs:
  - per-gene "best fit" table: which (alpha, tau) minimizes
    |CoFold_MFE − DMS_eval| for that gene.
  - per-gene heatmap data (alpha × tau → gap).
  - optional per-window correlation: for each (alpha, tau) and each gene,
    Pearson r between the CoFold per-window ΔG track and the DMS-projected
    per-window ΔG track. This answers "does the soft penalty improve the
    *shape* of the energy profile, not just its magnitude?"

The whole point of this stage is the user's hypothesis: "by penalizing
long-range pairs, Vienna ΔG mimics DMS ΔG more closely." This stage
quantifies that hypothesis at multiple parameter settings.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

import numpy as np
import pandas as pd

from mtrnafeat.analysis.window import sliding_intervals
from mtrnafeat.config import Config
from mtrnafeat.constants import canonical_gene
from mtrnafeat.core import thermo
from mtrnafeat.core.cofold import cofold_dG
from mtrnafeat.core.projection import project_structure_to_window
from mtrnafeat.core.structure import filter_max_bp_span
from mtrnafeat.io.db_parser import parse_db
from mtrnafeat.progress import progress, step


@dataclass(frozen=True)
class _GeneJob:
    species: str
    gene: str
    sequence: str
    structure: str
    max_bp_span: int
    alphas: tuple[float, ...]
    taus: tuple[float, ...]
    window_nt: int
    window_step: int
    do_window_corr: bool


def _rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def _gene_full_grid(job: _GeneJob) -> pd.DataFrame:
    rna = _rna(job.sequence)
    dms_struct = job.structure or ""
    dms_eval = float("nan")
    if dms_struct:
        try:
            dms_eval = thermo.eval_structure(rna, dms_struct, max_bp_span=job.max_bp_span)
        except Exception:
            dms_eval = float("nan")

    rows = []
    for alpha in job.alphas:
        for tau in job.taus:
            try:
                _, mfe = cofold_dG(rna, alpha=alpha, tau=tau, max_bp_span=job.max_bp_span)
            except Exception:
                mfe = float("nan")
            gap = float(mfe - dms_eval) if (np.isfinite(mfe) and np.isfinite(dms_eval)) else float("nan")
            rows.append({
                "Species": job.species,
                "Gene": canonical_gene(job.gene),
                "alpha": float(alpha),
                "tau": float(tau),
                "CoFold_MFE": float(mfe),
                "DMS_Eval_dG": float(dms_eval),
                "Gap_CoFold_minus_DMS": gap,
                "Abs_Gap": float(abs(gap)) if np.isfinite(gap) else float("nan"),
            })
    return pd.DataFrame(rows)


def _gene_window_corr(job: _GeneJob) -> pd.DataFrame:
    """Per (alpha, tau): Pearson r between per-window CoFold and DMS-projected ΔG."""
    rna = _rna(job.sequence)
    intervals = sliding_intervals(len(rna), job.window_nt, job.window_step)
    if len(intervals) < 3:
        return pd.DataFrame()

    # Pre-compute the DMS-projected per-window ΔG track (same for every alpha, tau).
    dms_track = []
    valid = []
    for s, e in intervals:
        seq_w = rna[s:e]
        dms_proj = project_structure_to_window(job.structure or "", s, e)
        dms_proj_span, _ = filter_max_bp_span(dms_proj, job.max_bp_span)
        try:
            dms_e = thermo.eval_structure(seq_w, dms_proj_span, max_bp_span=job.max_bp_span)
            dms_track.append(float(dms_e))
            valid.append(True)
        except Exception:
            dms_track.append(float("nan"))
            valid.append(False)
    dms_arr = np.array(dms_track)

    rows = []
    for alpha in job.alphas:
        for tau in job.taus:
            track = []
            for s, e in intervals:
                seq_w = rna[s:e]
                try:
                    _, mfe = cofold_dG(seq_w, alpha=alpha, tau=tau, max_bp_span=job.max_bp_span)
                    track.append(float(mfe))
                except Exception:
                    track.append(float("nan"))
            arr = np.array(track)
            mask = np.isfinite(arr) & np.isfinite(dms_arr)
            if mask.sum() >= 3:
                # Pearson r.
                a = arr[mask]; b = dms_arr[mask]
                if a.std() > 0 and b.std() > 0:
                    r = float(np.corrcoef(a, b)[0, 1])
                else:
                    r = float("nan")
                rmse = float(np.sqrt(np.mean((a - b) ** 2)))
                mean_gap = float(np.mean(a - b))
            else:
                r = rmse = mean_gap = float("nan")
            rows.append({
                "Species": job.species,
                "Gene": canonical_gene(job.gene),
                "alpha": float(alpha),
                "tau": float(tau),
                "Pearson_r_CoFold_vs_DMS": r,
                "Per_Window_RMSE": rmse,
                "Per_Window_Mean_Gap": mean_gap,
                "N_Windows_Valid": int(mask.sum()),
            })
    return pd.DataFrame(rows)


def _run_one(job: _GeneJob) -> tuple[pd.DataFrame, pd.DataFrame]:
    full = _gene_full_grid(job)
    win = _gene_window_corr(job) if job.do_window_corr else pd.DataFrame()
    return full, win


def run_cofold_sweep(cfg: Config, do_window_corr: bool = True
                      ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (full_grid_df, window_corr_df)."""
    step(f"running cofold-sweep ({len(cfg.cofold_alpha_sweep)}×{len(cfg.cofold_tau_sweep)} grid"
         + ", per-window corr ON" if do_window_corr else ", per-window corr OFF" + ")")
    jobs: list[_GeneJob] = []
    for species, fname in cfg.db_files.items():
        rec_by_gene = {r.gene: r for r in parse_db(cfg.data_dir / fname)}
        for gene in cfg.target_genes:
            target = canonical_gene(gene)
            if target not in rec_by_gene:
                continue
            rec = rec_by_gene[target]
            jobs.append(_GeneJob(
                species=species, gene=target,
                sequence=rec.sequence, structure=rec.structure or "",
                max_bp_span=int(cfg.max_bp_span),
                alphas=tuple(cfg.cofold_alpha_sweep),
                taus=tuple(cfg.cofold_tau_sweep),
                window_nt=int(cfg.window_nt),
                window_step=int(cfg.step_nt),
                do_window_corr=do_window_corr,
            ))
    if not jobs:
        return pd.DataFrame(), pd.DataFrame()

    workers = max(1, int(cfg.n_workers))
    full_frames: list[pd.DataFrame] = []
    win_frames: list[pd.DataFrame] = []
    if workers == 1 or len(jobs) == 1:
        for j in progress(jobs, desc="cofold-sweep (genes)", unit="gene"):
            f, w = _run_one(j)
            full_frames.append(f)
            if not w.empty:
                win_frames.append(w)
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(_run_one, j): j for j in jobs}
            for fut in progress(as_completed(futures), desc="cofold-sweep (genes)",
                                  total=len(futures), unit="gene"):
                f, w = fut.result()
                full_frames.append(f)
                if not w.empty:
                    win_frames.append(w)

    full = pd.concat(full_frames, ignore_index=True) if full_frames else pd.DataFrame()
    win = pd.concat(win_frames, ignore_index=True) if win_frames else pd.DataFrame()
    return full, win


def best_per_gene(full: pd.DataFrame) -> pd.DataFrame:
    """Per-gene row with the (alpha, tau) that minimizes |CoFold_MFE − DMS_eval|."""
    if full.empty:
        return pd.DataFrame()
    rows = []
    for (species, gene), g in full.groupby(["Species", "Gene"]):
        valid = g.dropna(subset=["Abs_Gap"])
        if valid.empty:
            continue
        best = valid.sort_values("Abs_Gap").iloc[0]
        rows.append({
            "Species": species,
            "Gene": gene,
            "Best_alpha": float(best["alpha"]),
            "Best_tau": float(best["tau"]),
            "Best_CoFold_MFE": float(best["CoFold_MFE"]),
            "DMS_Eval_dG": float(best["DMS_Eval_dG"]),
            "Best_Abs_Gap": float(best["Abs_Gap"]),
            "Plain_Vienna_at_alpha0": float(g[(g["alpha"] == 0.0)].iloc[0]["CoFold_MFE"])
                                          if (g["alpha"] == 0.0).any() else float("nan"),
        })
    return pd.DataFrame(rows)
