"""
The method, per sample: QC metrics, thresholds (MAD / fixed / none), filtering
of observations and genes, normalization (counts -> layers["counts"], X ->
normalize_total + log1p), highly variable genes (seurat_v3 on counts) and PCA.

Observations are never reordered: filtering is a boolean mask in file order.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

log = logging.getLogger("st-qc")

KINDS = ("spot", "bin", "cell")
REASONS = ("low_counts", "low_genes", "high_mt", "high_control")
# Fixed-mode defaults when MIN_COUNTS / MIN_GENES are 0, by kind of observation.
FIXED_DEFAULTS = {"spot": (500, 250), "bin": (50, 20), "cell": (10, 5)}


@dataclass
class Thresholds:
    mode: str
    kind: str
    min_counts: float | None = None          # observations with total_counts below this go
    min_genes: float | None = None
    max_pct_mt: float | None = None          # pct_counts_mt above this goes (MAD mode: and above the MAD cut)
    max_control_fraction: float | None = None
    details: list[dict] = field(default_factory=list)   # one row per metric for thresholds.csv

    def as_dict(self) -> dict:
        return {"mode": self.mode, "kind": self.kind, "min_counts": self.min_counts, "min_genes": self.min_genes,
                "max_pct_mt": self.max_pct_mt, "max_control_fraction": self.max_control_fraction}


def kind_of(adata: AnnData, row: dict | None = None) -> str:
    """spot / bin / cell: from uns["pennsieve"], else samples.csv, else guessed from obs."""
    meta = adata.uns.get("pennsieve") if isinstance(adata.uns.get("pennsieve"), dict) else {}
    kind = str(meta.get("kind") or (row or {}).get("kind") or "").strip().lower()
    if kind in KINDS:
        return kind
    if "control_counts" in adata.obs or "cell_area" in adata.obs:
        return "cell"
    if "bin_size_um" in meta:
        return "bin"
    return "spot"


def is_integer_matrix(x) -> bool:
    data = x.data if sp.issparse(x) else np.asarray(x)
    if data.size == 0:
        return True
    if np.issubdtype(data.dtype, np.integer):
        return True
    sample = data if data.size <= 2_000_000 else data[:: max(1, data.size // 2_000_000)]
    return bool(np.all(np.mod(sample, 1) == 0)) and bool(np.all(np.mod(data, 1) == 0))


def compute_qc_metrics(adata: AnnData) -> None:
    """scanpy.pp.calculate_qc_metrics with mt genes (var["mt"], else MT-/mt- prefix), plus
    obs["control_fraction"] = control_counts / (total_counts + control_counts) when present."""
    if "mt" not in adata.var:
        adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    adata.var["mt"] = adata.var["mt"].astype(bool)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "total_counts_mt"):
        adata.obs[c] = adata.obs[c].astype(np.float32 if c == "pct_counts_mt" else np.float64)
    if "control_counts" in adata.obs:
        ctrl = pd.to_numeric(adata.obs["control_counts"], errors="coerce").fillna(0.0).to_numpy(dtype=np.float64)
        total = adata.obs["total_counts"].to_numpy(dtype=np.float64) + ctrl
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total > 0, ctrl / total, 0.0)
        adata.obs["control_fraction"] = frac.astype(np.float32)


def _mad(values: np.ndarray) -> tuple[float, float]:
    med = float(np.median(values))
    return med, float(np.median(np.abs(values - med)))


def choose_thresholds(adata: AnnData, params: dict, kind: str) -> Thresholds:
    mode = params["FILTER_MODE"]
    thr = Thresholds(mode=mode, kind=kind)
    has_control = "control_fraction" in adata.obs
    if mode == "none":
        return thr
    obs = adata.obs
    if mode == "fixed":
        d_counts, d_genes = FIXED_DEFAULTS[kind]
        thr.min_counts = float(params["MIN_COUNTS"] or d_counts)
        thr.min_genes = float(params["MIN_GENES"] or d_genes)
        thr.max_pct_mt = float(params["MAX_PCT_MT"])
        thr.details = [
            {"metric": "total_counts", "direction": "below", "threshold": thr.min_counts, "source": "MIN_COUNTS" if params["MIN_COUNTS"] else f"default for {kind}"},
            {"metric": "n_genes_by_counts", "direction": "below", "threshold": thr.min_genes, "source": "MIN_GENES" if params["MIN_GENES"] else f"default for {kind}"},
            {"metric": "pct_counts_mt", "direction": "above", "threshold": thr.max_pct_mt, "source": "MAX_PCT_MT"},
        ]
    else:   # mad
        k = params["MADS"]
        lc = np.log1p(obs["total_counts"].to_numpy(dtype=np.float64))
        lg = np.log1p(obs["n_genes_by_counts"].to_numpy(dtype=np.float64))
        mt = obs["pct_counts_mt"].to_numpy(dtype=np.float64)
        for name, metric, values, direction in (("min_counts", "total_counts", lc, "below"),
                                                ("min_genes", "n_genes_by_counts", lg, "below"),
                                                ("max_pct_mt", "pct_counts_mt", mt, "above")):
            med, mad = _mad(values)
            row = {"metric": metric, "direction": direction, "median": med, "mad": mad, "mads": k}
            if direction == "above":
                # Both conditions: above the MAD cut and above MAX_PCT_MT (also when the MAD is 0,
                # e.g. panels without mitochondrial genes, where the cut is simply MAX_PCT_MT).
                cut = float(max(med + k * mad, params["MAX_PCT_MT"]))
                thr.max_pct_mt = cut
                row.update(threshold=cut, source=f"max(median + {k:g} MAD, MAX_PCT_MT={params['MAX_PCT_MT']:g})")
            elif mad == 0:
                row.update(threshold=None, source="skipped: MAD is 0")
                log.warning("%s: MAD of %s is 0; no %s-%s cut", adata.obs["sample_id"].iloc[0] if "sample_id" in adata.obs else "", metric, direction, metric)
            else:
                cut = float(np.expm1(med - k * mad))
                setattr(thr, name, cut)
                row.update(threshold=cut, source=f"expm1(median - {k:g} MAD) of log1p({metric})")
            thr.details.append(row)
    if has_control:
        thr.max_control_fraction = float(params["MAX_CONTROL_FRACTION"])
        thr.details.append({"metric": "control_fraction", "direction": "above", "threshold": thr.max_control_fraction,
                            "source": "MAX_CONTROL_FRACTION"})
    return thr


def flag_observations(adata: AnnData, thr: Thresholds) -> pd.DataFrame:
    """One boolean column per reason, in obs order."""
    obs = adata.obs
    n = adata.n_obs
    flags = pd.DataFrame({r: np.zeros(n, dtype=bool) for r in REASONS}, index=obs.index)
    if thr.mode == "none":
        return flags
    if thr.min_counts is not None:
        flags["low_counts"] = obs["total_counts"].to_numpy() < thr.min_counts
    if thr.min_genes is not None:
        flags["low_genes"] = obs["n_genes_by_counts"].to_numpy() < thr.min_genes
    if thr.max_pct_mt is not None:
        flags["high_mt"] = obs["pct_counts_mt"].to_numpy() > thr.max_pct_mt
    if thr.max_control_fraction is not None and "control_fraction" in obs:
        flags["high_control"] = obs["control_fraction"].to_numpy() > thr.max_control_fraction
    return flags


def reasons_of(flags: pd.DataFrame) -> pd.Series:
    """'low_counts;high_mt' per flagged observation (empty when kept)."""
    out = pd.Series([""] * len(flags), index=flags.index, dtype=object)
    for r in REASONS:
        hit = flags[r].to_numpy()
        out[hit] = np.where(out[hit] == "", r, out[hit] + ";" + r)
    return out


def filter_observations(adata: AnnData, keep: np.ndarray) -> None:
    if keep.all():
        return
    adata._inplace_subset_obs(np.flatnonzero(keep))


def filter_genes(adata: AnnData, min_cells: int) -> int:
    """Drops genes detected in fewer than min_cells observations; returns how many went."""
    if min_cells <= 0:
        return 0
    before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return before - adata.n_vars


def stage_counts(adata: AnnData) -> None:
    """layers["counts"] = the raw counts (int32 CSR when integer-valued); X becomes a float32 CSR copy to normalize."""
    x = adata.X
    if not sp.issparse(x):
        x = sp.csr_matrix(x)
    elif not isinstance(x, sp.csr_matrix):
        x = sp.csr_matrix(x)
    x.sum_duplicates()
    if is_integer_matrix(x):
        counts = x if x.dtype == np.int32 else x.astype(np.int32)
    else:
        log.warning("X is not integer-valued; layers['counts'] keeps it as float32")
        counts = x if x.dtype == np.float32 else x.astype(np.float32)
    adata.layers["counts"] = counts
    adata.X = x if (x.dtype == np.float32 and x is not counts) else counts.astype(np.float32)


def normalize(adata: AnnData, target_sum: float) -> float:
    """X = log1p(normalize_total(X, target_sum)), float32 CSR. Returns the target sum used
    (the median of the total counts when target_sum is 0, scanpy's default)."""
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    ts = float(target_sum) if target_sum > 0 else float(np.median(totals[totals > 0])) if np.any(totals > 0) else 1.0
    sc.pp.normalize_total(adata, target_sum=ts)
    sc.pp.log1p(adata)
    if adata.X.dtype != np.float32:
        adata.X = adata.X.astype(np.float32)
    return ts


def highly_variable(adata: AnnData, n_hvg: int) -> str:
    """var["highly_variable"]: seurat_v3 on layers["counts"] when the panel has more than n_hvg genes,
    else every gene. Returns how it was chosen."""
    if adata.n_vars <= n_hvg or n_hvg <= 0:
        adata.var["highly_variable"] = True
        return f"all {adata.n_vars} genes (panel has at most N_HVG={n_hvg})"
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_hvg, layer="counts", check_values=False)
    adata.var["highly_variable"] = adata.var["highly_variable"].astype(bool)
    return f"seurat_v3, top {n_hvg} of {adata.n_vars} genes on layers['counts']"


def pca(adata: AnnData, n_pcs: int, seed: int) -> int:
    """PCA on the log-normalized X of the highly variable genes (no scaling), seeded; obsm["X_pca"], uns["pca"].
    Returns the number of components computed (0 when the sample is too small)."""
    n_hvg = int(adata.var["highly_variable"].sum())
    n_comps = min(n_pcs, adata.n_obs - 1, n_hvg - 1)
    if n_comps < 1:
        log.warning("too few observations (%d) or genes (%d) for PCA", adata.n_obs, n_hvg)
        return 0
    sc.pp.pca(adata, n_comps=n_comps, mask_var="highly_variable", zero_center=True, svd_solver="arpack",
              random_state=seed, dtype="float32")
    return n_comps


def summarise_obs(adata: AnnData) -> dict:
    obs = adata.obs
    out = {"median_counts": float(np.median(obs["total_counts"])) if adata.n_obs else None,
           "median_genes": float(np.median(obs["n_genes_by_counts"])) if adata.n_obs else None,
           "median_pct_mt": float(np.median(obs["pct_counts_mt"])) if adata.n_obs else None}
    if "control_fraction" in obs and adata.n_obs:
        out["median_control_fraction"] = float(np.median(obs["control_fraction"]))
    return out
