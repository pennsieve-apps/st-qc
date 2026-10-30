"""
st-qc: quality control and normalization of spatial transcriptomics
samples, the step after st-prep in every spatial workflow.

Reads st-prep's <sample_id>.h5ad files (port `spatial`) and samples.csv (port
`metadata`), one sample at a time, and writes

  spatial/<sample_id>.h5ad   filtered, X log-normalized, layers["counts"], HVGs, PCA
  qc/qc_summary.csv          one row per sample: counts before/after, medians, thresholds
  qc/thresholds.csv          one row per sample x metric: how each cut was chosen
  qc/removed.csv             sample_id, obs_id, reason for every removed observation
  qc/<sample_id>_qc.png, qc/<sample_id>_spatial_counts.png
  qc/summary.json, provenance.json (in both output folders)

See processor/qc.py for the method.
"""
from __future__ import annotations

import csv
import json
import logging
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("st-qc")
for noisy in ("matplotlib", "PIL", "fontTools", "numba", "h5py"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-qc"
REPO_ROOT = Path(__file__).resolve().parent.parent

DEFAULTS = {
    "FILTER_MODE": "mad",
    "MADS": "5",
    "MIN_COUNTS": "0",
    "MIN_GENES": "0",
    "MAX_PCT_MT": "20",
    "MAX_CONTROL_FRACTION": "0.1",
    "MIN_CELLS_PER_GENE": "3",
    "TARGET_SUM": "0",
    "N_HVG": "2000",
    "N_PCS": "50",
    "SEED": "42",
    "THREADS": "0",
}
PARAMETERS = tuple(DEFAULTS)
FILTER_MODES = ("mad", "fixed", "none")


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (env.get(k) if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    p: dict = {}
    p["FILTER_MODE"] = raw["FILTER_MODE"].strip().lower()
    if p["FILTER_MODE"] not in FILTER_MODES:
        raise ParameterError(f"FILTER_MODE must be one of {', '.join(FILTER_MODES)}, got {raw['FILTER_MODE']!r}")
    try:
        for k in ("MADS", "MIN_COUNTS", "MIN_GENES", "MAX_PCT_MT", "MAX_CONTROL_FRACTION", "TARGET_SUM"):
            p[k] = float(raw[k])
        for k in ("MIN_CELLS_PER_GENE", "N_HVG", "N_PCS", "SEED", "THREADS"):
            p[k] = int(float(raw[k]))
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if p["MADS"] <= 0:
        raise ParameterError("MADS must be > 0")
    if min(p["MIN_COUNTS"], p["MIN_GENES"], p["TARGET_SUM"], p["MIN_CELLS_PER_GENE"], p["N_HVG"], p["THREADS"]) < 0:
        raise ParameterError("MIN_COUNTS, MIN_GENES, TARGET_SUM, MIN_CELLS_PER_GENE, N_HVG and THREADS must be >= 0")
    if not 0 <= p["MAX_PCT_MT"] <= 100:
        raise ParameterError("MAX_PCT_MT is a percentage, 0-100")
    if not 0 <= p["MAX_CONTROL_FRACTION"] <= 1:
        raise ParameterError("MAX_CONTROL_FRACTION is a fraction, 0-1")
    if p["N_PCS"] < 1:
        raise ParameterError("N_PCS must be >= 1")
    return p


def _rel(path: Path | None, root: Path) -> str | None:
    if path is None:
        return None
    try:
        return str(Path(path).relative_to(root))
    except ValueError:
        return str(path)


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return ""
    if isinstance(v, float):
        return "" if math.isnan(v) else f"{v:.6g}"
    return v


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r.get(c, "")) for c in columns])


def _clean_json(obj):
    import numpy as np

    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        v = float(obj)
        return None if math.isnan(v) or math.isinf(v) else v
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, np.ndarray)):
        return [_clean_json(v) for v in obj]
    return obj


def write_json(path: Path, data) -> None:
    with open(path, "w") as fh:
        json.dump(_clean_json(data), fh, indent=2)
        fh.write("\n")


def tool_versions() -> dict:
    from importlib.metadata import version

    out = {"python": platform.python_version()}
    for name in ("scanpy", "anndata", "numpy", "scipy", "pandas", "h5py", "scikit-learn", "scikit-misc", "matplotlib"):
        try:
            out[name] = version(name)
        except Exception:   # noqa: BLE001
            out[name] = None
    return out


def peak_rss_mb() -> float | None:
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def cgroup_peak_mb() -> float | None:
    try:
        return round(int(open("/sys/fs/cgroup/memory.peak").read()) / 1024 / 1024, 1)
    except (OSError, ValueError):
        return None


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def process_sample(sample, params: dict, out_spatial: Path, out_qc: Path) -> dict:
    """One sample through QC; returns the qc_summary row plus the threshold and removed rows."""
    import anndata as ad
    import numpy as np

    from processor import figures as F
    from processor import qc as Q

    t0 = time.monotonic()
    adata = ad.read_h5ad(sample.path)
    sid = sample.sample_id
    if "sample_id" not in adata.obs:
        adata.obs["sample_id"] = sid
    kind = Q.kind_of(adata, sample.row)
    meta = adata.uns.get("pennsieve") if isinstance(adata.uns.get("pennsieve"), dict) else {}
    n_obs_before, n_genes_before = adata.n_obs, adata.n_vars
    if n_obs_before == 0 or n_genes_before == 0:
        raise ValueError(f"empty matrix ({n_obs_before} observations x {n_genes_before} genes)")
    if "spatial" not in adata.obsm:
        raise ValueError("no obsm['spatial']")
    t_read = time.monotonic() - t0

    Q.compute_qc_metrics(adata)
    thr = Q.choose_thresholds(adata, params, kind)
    flags = Q.flag_observations(adata, thr)
    removed = flags.any(axis=1).to_numpy()
    reasons = Q.reasons_of(flags)
    n_by_reason = {r: int(flags[r].sum()) for r in Q.REASONS}
    obs_before = adata.obs[[c for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "control_fraction") if c in adata.obs]].copy()
    xy_before = np.asarray(adata.obsm["spatial"])
    uns_spatial = adata.uns.get("spatial")
    n_after = int((~removed).sum())
    if n_after < 2:
        raise ValueError(f"{n_after} of {n_obs_before} observations would be kept; check the thresholds")
    t1 = time.monotonic()
    Q.filter_observations(adata, ~removed)
    n_genes_dropped = Q.filter_genes(adata, params["MIN_CELLS_PER_GENE"])
    if adata.n_vars < 2:
        raise ValueError(f"{adata.n_vars} genes left after MIN_CELLS_PER_GENE={params['MIN_CELLS_PER_GENE']}")
    if n_genes_dropped:
        Q.compute_qc_metrics(adata)   # obs metrics describe the matrix that is written
    stats = Q.summarise_obs(adata)
    t_filter = time.monotonic() - t1

    t2 = time.monotonic()
    Q.stage_counts(adata)
    hvg_how = Q.highly_variable(adata, params["N_HVG"])   # on the counts, before X is normalized
    target_sum = Q.normalize(adata, params["TARGET_SUM"])
    n_comps = Q.pca(adata, params["N_PCS"], params["SEED"])
    t_norm = time.monotonic() - t2

    meta = dict(meta)
    meta.update({"kind": kind, "normalized": True, "target_sum": float(target_sum), "qc_version": app_version(),
                 "qc_filter_mode": params["FILTER_MODE"]})
    adata.uns["pennsieve"] = meta
    t3 = time.monotonic()
    adata.write_h5ad(out_spatial / f"{sid}.h5ad", compression="gzip")
    t_write = time.monotonic() - t3

    t4 = time.monotonic()
    figs = []
    try:
        figs.append(F.qc_histograms(sid, kind, obs_before, removed, thr, n_obs_before, n_after, out_qc / f"{sid}_qc.png"))
        before = ad.AnnData(obs=obs_before, obsm={"spatial": xy_before})
        if uns_spatial is not None:
            before.uns["spatial"] = uns_spatial
        figs.append(F.spatial_counts(sid, before, removed, reasons, out_qc / f"{sid}_spatial_counts.png"))
    except Exception as exc:   # figures are a convenience; the tables and the h5ad stand without them
        log.exception("%s: figures failed: %s", sid, exc)
    t_fig = time.monotonic() - t4

    row = {"sample_id": sid, "kind": kind, "platform": meta.get("platform") or sample.row.get("platform", ""),
           "n_obs_before": n_obs_before, "n_obs_after": n_after, "n_genes_before": n_genes_before,
           "n_genes_after": adata.n_vars, "pct_removed": 100.0 * (n_obs_before - n_after) / n_obs_before,
           **stats, "filter_mode": thr.mode, "min_counts": thr.min_counts, "min_genes": thr.min_genes,
           "max_pct_mt": thr.max_pct_mt, "max_control_fraction": thr.max_control_fraction,
           "n_removed_low_counts": n_by_reason["low_counts"], "n_removed_low_genes": n_by_reason["low_genes"],
           "n_removed_mt": n_by_reason["high_mt"], "n_removed_control": n_by_reason["high_control"],
           "n_genes_dropped": n_genes_dropped, "n_hvg": int(adata.var["highly_variable"].sum()), "n_pcs": n_comps,
           "target_sum": target_sum, "status": "ok"}
    thr_rows = [{"sample_id": sid, "kind": kind, "mode": thr.mode, **d} for d in thr.details]
    removed_rows = [{"sample_id": sid, "obs_id": oid, "reason": why}
                    for oid, why in zip(reasons.index[removed], reasons.to_numpy()[removed])]
    timings = {"read": t_read, "filter": t_filter, "normalize_hvg_pca": t_norm, "write": t_write, "figures": t_fig,
               "total": time.monotonic() - t0}
    log.info("%s: %s %s, %d -> %d obs (%.1f %% removed: %s), %d -> %d genes, median counts %.0f, HVG %s, %d PCs, "
             "target_sum %.0f, %.1fs (read %.1f, filter %.1f, norm+pca %.1f, write %.1f, figures %.1f)",
             sid, kind, meta.get("platform", ""), n_obs_before, n_after, row["pct_removed"],
             ", ".join(f"{k} {v}" for k, v in n_by_reason.items() if v) or "none", n_genes_before, adata.n_vars,
             stats["median_counts"] or 0, hvg_how, n_comps, target_sum, timings["total"], t_read, t_filter, t_norm,
             t_write, t_fig)
    return {"row": row, "thresholds": thr_rows, "removed": removed_rows, "figures": [f.name for f in figs],
            "timings": {k: round(v, 2) for k, v in timings.items()}, "hvg": hvg_how, "file": sample.path}


SUMMARY_COLUMNS = ["sample_id", "kind", "platform", "n_obs_before", "n_obs_after", "pct_removed", "n_genes_before",
                   "n_genes_after", "median_counts", "median_genes", "median_pct_mt", "median_control_fraction",
                   "filter_mode", "min_counts", "min_genes", "max_pct_mt", "max_control_fraction",
                   "n_removed_low_counts", "n_removed_low_genes", "n_removed_mt", "n_removed_control",
                   "n_genes_dropped", "n_hvg", "n_pcs", "target_sum", "status"]
THRESHOLD_COLUMNS = ["sample_id", "kind", "mode", "metric", "direction", "threshold", "median", "mad", "mads", "source"]


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2
    threads = os.environ.get("ST_QC_THREADS") or os.environ.get("OMP_NUM_THREADS") or "unset"
    log.info("parameters: %s; threads %s (%s)", {k: v for k, v in params.items()}, threads,
             os.environ.get("ST_QC_THREADS_SOURCE", "environment"))

    from processor.inputs import discover

    try:
        found = discover(input_dir)
    except (OSError, ValueError) as exc:
        log.error("cannot read the inputs: %s", exc)
        return 1
    for name, reason in found.skipped.items():
        log.warning("%s: skipped, %s", name, reason)
    if not found.samples:
        log.error("no <sample_id>.h5ad found under %s", input_dir)
        return 1
    out_spatial, out_qc = output_dir / "spatial", output_dir / "qc"
    out_spatial.mkdir(parents=True, exist_ok=True)
    out_qc.mkdir(parents=True, exist_ok=True)
    log.info("%d samples (samples.csv: %s)", len(found.samples), _rel(found.samples_csv, input_dir))

    rows, thr_rows, removed_rows, statuses, figures = [], [], [], [], []
    for s in found.samples:
        t = time.monotonic()
        status = {"sample_id": s.sample_id, "file": _rel(s.path, input_dir)}
        try:
            result = process_sample(s, params, out_spatial, out_qc)
            rows.append(result["row"])
            thr_rows.extend(result["thresholds"])
            removed_rows.extend(result["removed"])
            figures.extend(result["figures"])
            r = result["row"]
            status.update(status="ok", n_obs_before=r["n_obs_before"], n_obs_after=r["n_obs_after"],
                          n_genes_before=r["n_genes_before"], n_genes_after=r["n_genes_after"], hvg=result["hvg"],
                          timings_seconds=result["timings"], peak_rss_mb=peak_rss_mb())
        except Exception as exc:   # one bad sample is skipped, not fatal
            status.update(status="failed", error=f"{type(exc).__name__}: {exc}", seconds=round(time.monotonic() - t, 2))
            log.exception("%s: failed: %s", s.sample_id, exc)
            rows.append({"sample_id": s.sample_id, "status": f"failed: {type(exc).__name__}: {exc}"})
        statuses.append(status)
    n_ok = sum(1 for s in statuses if s["status"] == "ok")
    n_failed = len(statuses) - n_ok

    write_csv(out_qc / "qc_summary.csv", rows, SUMMARY_COLUMNS)
    write_csv(out_qc / "thresholds.csv", thr_rows, THRESHOLD_COLUMNS)
    write_csv(out_qc / "removed.csv", removed_rows, ["sample_id", "obs_id", "reason"])
    ok_rows = [r for r in rows if r.get("status") == "ok"]
    summary = {
        "n_samples": len(statuses), "n_ok": n_ok, "n_failed": n_failed,
        "n_obs_before": sum(r["n_obs_before"] for r in ok_rows), "n_obs_after": sum(r["n_obs_after"] for r in ok_rows),
        "n_removed": len(removed_rows),
        "n_removed_by_reason": {k: sum(r[f"n_removed_{k}"] for r in ok_rows)
                                for k in ("low_counts", "low_genes", "mt", "control")},
        "kinds": sorted({r["kind"] for r in ok_rows}),
        "samples": [{k: r.get(k) for k in SUMMARY_COLUMNS} for r in rows],
        "parameters": params, "figures": sorted(figures),
    }
    write_json(out_qc / "summary.json", summary)
    provenance = {
        "app": APP_NAME, "version": app_version(), "tools": tool_versions(), "parameters": params,
        "threads": {"n": threads, "source": os.environ.get("ST_QC_THREADS_SOURCE", "environment")},
        "method": {
            "qc_metrics": "scanpy.pp.calculate_qc_metrics(qc_vars=['mt']); control_fraction = control_counts / (total_counts + control_counts)",
            "filter": {"mad": "outlier when log1p(total_counts) or log1p(n_genes_by_counts) < median - MADS*MAD, or pct_counts_mt > max(median + MADS*MAD, MAX_PCT_MT); control_fraction > MAX_CONTROL_FRACTION when present (MAD unscaled, as in scipy.stats.median_abs_deviation)",
                       "fixed": "total_counts < MIN_COUNTS, n_genes_by_counts < MIN_GENES (kind defaults when 0: spot 500/250, bin 50/20, cell 10/5), pct_counts_mt > MAX_PCT_MT, control_fraction > MAX_CONTROL_FRACTION",
                       "none": "nothing removed"}[params["FILTER_MODE"]],
            "genes": f"dropped when detected in fewer than MIN_CELLS_PER_GENE={params['MIN_CELLS_PER_GENE']} kept observations",
            "normalization": "layers['counts'] = raw counts; X = log1p(normalize_total(target_sum=TARGET_SUM or median total counts)), float32 CSR",
            "hvg": f"scanpy.pp.highly_variable_genes(flavor='seurat_v3', n_top_genes={params['N_HVG']}, layer='counts') when the panel has more genes, else all",
            "pca": f"scanpy.pp.pca(n_comps={params['N_PCS']}, mask_var='highly_variable', zero_center=True, svd_solver='arpack', random_state={params['SEED']}); X not scaled",
        },
        "inputs": {"spatial": [_rel(p, input_dir) for p in found.files_seen], "samples_csv": _rel(found.samples_csv, input_dir),
                   "skipped": found.skipped},
        "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
        "timings_seconds": {"total": round(time.monotonic() - t0, 2)},
        "peak_rss_mb": peak_rss_mb(), "cgroup_memory_peak_mb": cgroup_peak_mb(),
        "samples": statuses, "n_ok": n_ok, "n_failed": n_failed,
        "outputs": {"spatial": sorted(p.name for p in out_spatial.iterdir()) + ["provenance.json"],
                    "qc": sorted(p.name for p in out_qc.iterdir()) + ["provenance.json"]},
    }
    write_json(out_qc / "provenance.json", provenance)
    write_json(out_spatial / "provenance.json", provenance)
    log.info("done: %d samples (%d failed), %d -> %d observations, %.1fs, peak RSS %s MB, cgroup peak %s MB -> %s",
             n_ok, n_failed, summary["n_obs_before"], summary["n_obs_after"], time.monotonic() - t0,
             provenance["peak_rss_mb"], provenance["cgroup_memory_peak_mb"], output_dir)
    if n_ok == 0:
        log.error("every sample failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
