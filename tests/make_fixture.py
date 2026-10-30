"""
A synthetic st-prep output that follows the workflow contract
(spatial/<sample_id>.h5ad + metadata/{samples.csv, prep.json}), for the unit
tests and `make run`:

  python -m tests.make_fixture data/fixture [--layout linked|flat|both] [--seed 7]

Two samples with planted low-quality observations:

- `visium_a` (kind spot, platform visium): a 20 x 20 hex grid of 400 spots in
  full-resolution pixel coordinates with a lowres "H&E" image and scalefactors
  in uns["spatial"]; 300 genes incl. 13 MT- genes; a planted domain (the left
  half expresses one gene set more). 12 spots have ~1 % of the normal counts
  (low counts and few genes), 8 spots have ~60 % mitochondrial counts, and gene
  `UNDETECTED1` is never seen (dropped by MIN_CELLS_PER_GENE).
- `cells_a` (kind cell, platform xenium): 450 cells in microns, a 60-gene panel
  without MT genes, obs["control_counts"] and obs["cell_area"]. 10 cells have a
  control fraction of ~0.5, 10 cells have 2 counts.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

PLANTED = {
    "visium_a": {"low_counts": 12, "high_mt": 8},
    "cells_a": {"high_control": 10, "low_counts": 10},
}


def _counts(rng: np.random.Generator, n: int, mu: np.ndarray, size_factor: np.ndarray, dispersion: float = 0.3):
    lam = np.outer(size_factor, mu)
    shape = 1.0 / dispersion
    return rng.poisson(rng.gamma(shape, lam / shape)).astype(np.int32)


def make_visium(rng: np.random.Generator, sid: str = "visium_a") -> ad.AnnData:
    rows, cols = 20, 20
    array_row, array_col = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    array_row, array_col = array_row.ravel(), array_col.ravel() * 2 + (array_row.ravel() % 2)
    spacing = 200.0   # full-resolution pixels between spot centres
    x = 1000.0 + array_col * spacing / 2.0
    y = 1000.0 + array_row * spacing * 0.866
    n = rows * cols
    genes = [f"MT-G{i}" for i in range(13)] + [f"GENE{i}" for i in range(286)] + ["UNDETECTED1"]
    g = len(genes)
    mu = rng.gamma(2.0, 2.0, g)
    mu[:13] = 3.0           # mt genes ~5 % of the counts
    mu[-1] = 0.0
    left = (array_col < cols).astype(float)
    domain_genes = np.zeros(g)
    domain_genes[20:40] = 1.0
    mu_obs = mu[None, :] * (1.0 + 2.0 * np.outer(left, domain_genes))
    size = rng.lognormal(0.0, 0.25, n)
    shape = 1.0 / 0.3
    X = rng.poisson(rng.gamma(shape, size[:, None] * mu_obs / shape)).astype(np.int32)
    # planted: 12 spots with ~1 % of the counts, 8 spots with ~60 % mitochondrial counts
    low = rng.choice(n, 12, replace=False)
    X[low] = rng.poisson(X[low] * 0.01)
    rest = np.setdiff1d(np.arange(n), low)
    mt_high = rng.choice(rest, 8, replace=False)
    X[np.ix_(mt_high, np.arange(13))] = rng.poisson(X[mt_high, 13:].sum(1, keepdims=True) * 1.5 / 13).astype(np.int32)
    obs = pd.DataFrame({
        "sample_id": pd.Categorical([sid] * n), "x": x, "y": y, "in_tissue": np.ones(n, dtype=np.int64),
        "array_row": array_row.astype(np.int64), "array_col": array_col.astype(np.int64),
        "condition": pd.Categorical(["kidney"] * n), "planted": pd.Categorical(np.where(np.isin(np.arange(n), low), "low_counts", np.where(np.isin(np.arange(n), mt_high), "high_mt", "")))},
        index=[f"{sid}:{'ACGT'[i % 4]}{i:05d}-1" for i in range(n)])
    var = pd.DataFrame({"gene_ids": [f"ENSG{i:011d}" for i in range(g)], "feature_types": "Gene Expression",
                        "mt": [name.startswith("MT-") for name in genes]}, index=genes)
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.obsm["spatial"] = np.column_stack([x, y]).astype(np.float64)
    lowres_scale = 0.05
    img = np.full((int(4500 * lowres_scale), int(4500 * lowres_scale), 3), 235, dtype=np.uint8)
    yy, xx = np.mgrid[0 : img.shape[0], 0 : img.shape[1]]
    tissue = ((xx - 125) ** 2 / 90.0**2 + (yy - 110) ** 2 / 85.0**2) < 1
    img[tissue] = (225, 170, 190)
    adata.uns["spatial"] = {sid: {"images": {"lowres": img},
                                  "scalefactors": {"tissue_hires_scalef": 0.2, "tissue_lowres_scalef": lowres_scale,
                                                   "spot_diameter_fullres": 130.0, "fiducial_diameter_fullres": 210.0}}}
    adata.uns["pennsieve"] = {"platform": "visium", "kind": "spot", "coordinate_unit": "pixel", "spot_diameter": 130.0,
                              "grid": "hex", "source": "visium_a/", "prep_version": "0.1.0"}
    return adata


def make_cells(rng: np.random.Generator, sid: str = "cells_a") -> ad.AnnData:
    n = 450
    genes = [f"PANEL{i}" for i in range(60)]
    g = len(genes)
    xy = rng.uniform(0, 800, (n, 2))
    mu = rng.gamma(1.5, 1.5, g)
    types = rng.integers(0, 3, n)
    mu_obs = mu[None, :] * (1.0 + 2.0 * (types[:, None] == (np.arange(g) % 3)[None, :]))
    size = rng.lognormal(0.0, 0.3, n)
    shape = 1.0 / 0.4
    X = rng.poisson(rng.gamma(shape, size[:, None] * mu_obs / shape)).astype(np.int32)
    totals = X.sum(1)
    control = rng.poisson(0.02 * totals + 0.5).astype(np.int64)
    ctrl_high = rng.choice(n, 10, replace=False)
    control[ctrl_high] = (totals[ctrl_high] * 1.0).astype(np.int64)   # control fraction ~0.5
    rest = np.setdiff1d(np.arange(n), ctrl_high)
    low = rng.choice(rest, 10, replace=False)
    X[low] = 0
    X[low, rng.integers(0, g, 10)] = 2
    control[low] = 0
    obs = pd.DataFrame({
        "sample_id": pd.Categorical([sid] * n), "x": xy[:, 0], "y": xy[:, 1], "in_tissue": np.ones(n, dtype=np.int64),
        "cell_area": rng.gamma(4.0, 20.0, n), "nucleus_area": rng.gamma(4.0, 6.0, n), "control_counts": control,
        "condition": pd.Categorical(["tumour"] * n),
        "planted": pd.Categorical(np.where(np.isin(np.arange(n), ctrl_high), "high_control", np.where(np.isin(np.arange(n), low), "low_counts", "")))},
        index=[f"{sid}:cell_{i}" for i in range(n)])
    var = pd.DataFrame({"gene_ids": genes, "feature_types": "Gene Expression", "mt": False}, index=genes)
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.obsm["spatial"] = xy.astype(np.float64)
    adata.uns["pennsieve"] = {"platform": "xenium", "kind": "cell", "coordinate_unit": "micron", "spot_diameter": 10.0,
                              "grid": "none", "source": "cells_a/", "prep_version": "0.1.0"}
    return adata


def make_fixture(root: Path, layout: str = "both", seed: int = 7) -> dict[str, Path]:
    rng = np.random.default_rng(seed)
    samples = [make_visium(rng), make_cells(rng)]
    linked = root / "linked"
    (linked / "spatial").mkdir(parents=True, exist_ok=True)
    (linked / "metadata").mkdir(parents=True, exist_ok=True)
    rows, prep = [], {"app": "st-prep", "samples": {}}
    for a in samples:
        sid = str(a.obs["sample_id"].iloc[0])
        a.write_h5ad(linked / "spatial" / f"{sid}.h5ad")
        meta = a.uns["pennsieve"]
        rows.append({"sample_id": sid, "path": meta["source"], "platform": meta["platform"], "kind": meta["kind"],
                     "n_obs": a.n_obs, "n_genes": a.n_vars, "n_obs_in_tissue": a.n_obs, "condition": str(a.obs["condition"].iloc[0])})
        prep["samples"][sid] = {"platform": meta["platform"], "files": [meta["source"]], "features_dropped": {},
                                "n_obs_dropped_out_of_tissue": 0, "coordinate_unit": meta["coordinate_unit"],
                                "image": meta["platform"] == "visium", "warnings": []}
    with open(linked / "metadata" / "samples.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(linked / "metadata" / "prep.json", "w") as fh:
        json.dump(prep, fh, indent=2)
    out = {"linked": linked}
    if layout in ("flat", "both"):
        flat = root / "flat"
        flat.mkdir(parents=True, exist_ok=True)
        for p in (linked / "spatial").iterdir():
            shutil.copy(p, flat / p.name)
        for p in (linked / "metadata").iterdir():
            shutil.copy(p, flat / p.name)
        out["flat"] = flat
    if layout == "flat":
        shutil.rmtree(linked)
        out.pop("linked")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    for name, path in make_fixture(args.out, args.layout, args.seed).items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
