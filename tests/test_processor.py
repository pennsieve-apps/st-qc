"""The synthetic fixture through the processor (both input layouts), the contract on the written h5ad, the modes."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from processor import main as app
from processor import qc as Q
from tests.make_fixture import PLANTED, make_fixture

QC_FILES = {"qc_summary.csv", "thresholds.csv", "removed.csv", "summary.json", "provenance.json"}


@pytest.fixture(scope="session")
def fixture_root(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv("MPLCONFIGDIR", "/tmp/matplotlib")
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, str(value))
    return app.main()


def rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def check_contract(src: ad.AnnData, out: ad.AnnData, n_pcs: int) -> None:
    # Observations: a subset of the input, in input order, never reordered.
    kept = src.obs_names.isin(out.obs_names)
    assert list(out.obs_names) == list(src.obs_names[kept])
    assert list(out.var_names) == list(src.var_names[src.var_names.isin(out.var_names)])
    # layers["counts"]: the input counts, integer, CSR.
    counts = out.layers["counts"]
    assert sp.isspmatrix_csr(counts) and np.issubdtype(counts.dtype, np.integer)
    expected = src[out.obs_names, out.var_names].X
    assert (counts != expected).nnz == 0
    # X: log1p(normalize_total) float32 CSR, totals equal target_sum.
    assert sp.isspmatrix_csr(out.X) and out.X.dtype == np.float32
    ts = out.uns["pennsieve"]["target_sum"]
    totals = np.asarray(np.expm1(out.X).sum(axis=1)).ravel()
    assert np.allclose(totals, ts, rtol=1e-3)
    assert out.uns["pennsieve"]["normalized"] is True or out.uns["pennsieve"]["normalized"] == True  # noqa: E712
    assert out.var["highly_variable"].dtype == bool
    assert out.obsm["X_pca"].shape == (out.n_obs, n_pcs) and "pca" in out.uns
    assert np.array_equal(out.obsm["spatial"], src[out.obs_names].obsm["spatial"])
    for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "sample_id", "x", "y"):
        assert c in out.obs


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_planted_outliers_are_removed(fixture_root, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixture_root[layout], out, monkeypatch) == 0
    assert sorted(p.name for p in out.iterdir()) == ["qc", "spatial"]
    assert QC_FILES <= {p.name for p in (out / "qc").iterdir()}
    for sid in PLANTED:
        assert (out / "qc" / f"{sid}_qc.png").exists() and (out / "qc" / f"{sid}_spatial_counts.png").exists()
    summary = {r["sample_id"]: r for r in rows(out / "qc" / "qc_summary.csv")}
    assert set(summary) == set(PLANTED) and all(r["status"] == "ok" for r in summary.values())
    assert summary["visium_a"]["kind"] == "spot" and summary["cells_a"]["kind"] == "cell"
    removed = rows(out / "qc" / "removed.csv")
    assert list(removed[0]) == ["sample_id", "obs_id", "reason"]

    for sid, planted in PLANTED.items():
        src = ad.read_h5ad(fixture_root[layout] / ("spatial" if layout == "linked" else "") / f"{sid}.h5ad")
        res = ad.read_h5ad(out / "spatial" / f"{sid}.h5ad")
        check_contract(src, res, 50)
        planted_ids = {why: set(src.obs_names[src.obs["planted"] == why]) for why in planted}
        got = {r["obs_id"]: r["reason"].split(";") for r in removed if r["sample_id"] == sid}
        for why, ids in planted_ids.items():
            assert len(ids) == planted[why]
            assert all(oid in got and why in got[oid] for oid in ids), f"{sid}: {why} not all removed for that reason"
        # Nothing else is removed, within a small tolerance of random outliers.
        all_planted = set().union(*planted_ids.values())
        assert len(set(got) - all_planted) <= 3
        assert int(summary[sid]["n_obs_before"]) == src.n_obs and int(summary[sid]["n_obs_after"]) == res.n_obs
        assert res.n_obs == src.n_obs - len(got)
        assert res.uns["pennsieve"]["kind"] == summary[sid]["kind"]
    # Low-count spots fail both the counts and the genes cut; the never-seen gene is dropped.
    visium_removed = {r["obs_id"]: r["reason"] for r in removed if r["sample_id"] == "visium_a"}
    assert any("low_counts;low_genes" == why for why in visium_removed.values())
    assert int(summary["visium_a"]["n_removed_low_counts"]) >= 12 and int(summary["visium_a"]["n_removed_mt"]) >= 8
    assert int(summary["visium_a"]["n_genes_dropped"]) >= 1
    res = ad.read_h5ad(out / "spatial" / "visium_a.h5ad")
    assert "UNDETECTED1" not in res.var_names and "spatial" in res.uns and "visium_a" in res.uns["spatial"]
    # The 60-gene panel: every gene highly variable; control fraction thresholds applied to cells only.
    cells = ad.read_h5ad(out / "spatial" / "cells_a.h5ad")
    assert cells.var["highly_variable"].all() and "control_fraction" in cells.obs
    assert int(summary["cells_a"]["n_removed_control"]) >= 10 and summary["visium_a"]["max_control_fraction"] == ""
    thr = rows(out / "qc" / "thresholds.csv")
    assert {(t["sample_id"], t["metric"]) for t in thr} >= {("visium_a", "total_counts"), ("visium_a", "pct_counts_mt"), ("cells_a", "control_fraction")}
    assert all(t["mode"] == "mad" for t in thr)

    js = json.loads((out / "qc" / "summary.json").read_text())
    assert js["n_ok"] == 2 and js["n_failed"] == 0 and js["n_removed"] == len(removed)
    prov = json.loads((out / "qc" / "provenance.json").read_text())
    assert prov["app"] == "st-qc" and prov["version"] == "0.1.0"
    assert {"scanpy", "anndata", "numpy"} <= set(prov["tools"])
    assert (out / "spatial" / "provenance.json").exists()


def test_filter_mode_none_keeps_everything(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, FILTER_MODE="none", MIN_CELLS_PER_GENE="0", N_PCS="10", N_HVG="100") == 0
    assert rows(out / "qc" / "removed.csv") == []
    for sid in PLANTED:
        src = ad.read_h5ad(fixture_root["linked"] / "spatial" / f"{sid}.h5ad")
        res = ad.read_h5ad(out / "spatial" / f"{sid}.h5ad")
        assert res.n_obs == src.n_obs and res.n_vars == src.n_vars
        check_contract(src, res, 10)
    # seurat_v3 picked 100 of the 300 visium genes; the 60-gene panel keeps all.
    assert int(ad.read_h5ad(out / "spatial" / "visium_a.h5ad").var["highly_variable"].sum()) == 100
    assert ad.read_h5ad(out / "spatial" / "cells_a.h5ad").var["highly_variable"].all()


def test_fixed_mode_uses_kind_defaults(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, FILTER_MODE="fixed", N_PCS="5") == 0
    summary = {r["sample_id"]: r for r in rows(out / "qc" / "qc_summary.csv")}
    assert float(summary["visium_a"]["min_counts"]) == 500 and float(summary["visium_a"]["min_genes"]) == 250
    assert float(summary["cells_a"]["min_counts"]) == 10 and float(summary["cells_a"]["min_genes"]) == 5
    assert float(summary["cells_a"]["max_control_fraction"]) == 0.1
    thr = rows(out / "qc" / "thresholds.csv")
    assert any(t["source"] == "default for spot" for t in thr)
    removed = rows(out / "qc" / "removed.csv")
    assert sum(1 for r in removed if r["sample_id"] == "cells_a" and "low_counts" in r["reason"]) >= 10
    assert run(fixture_root["linked"], tmp_path / "out2", monkeypatch, FILTER_MODE="fixed", MIN_COUNTS="100", N_PCS="5") == 0
    assert float({r["sample_id"]: r for r in rows(tmp_path / "out2" / "qc" / "qc_summary.csv")}["cells_a"]["min_counts"]) == 100


def test_target_sum(fixture_root, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixture_root["linked"], out, monkeypatch, TARGET_SUM="100", N_PCS="5") == 0
    res = ad.read_h5ad(out / "spatial" / "cells_a.h5ad")
    assert res.uns["pennsieve"]["target_sum"] == 100
    assert np.allclose(np.asarray(np.expm1(res.X).sum(axis=1)).ravel(), 100, rtol=1e-3)


@pytest.mark.parametrize("params", [{"FILTER_MODE": "strict"}, {"MADS": "0"}, {"MAX_PCT_MT": "150"}, {"N_PCS": "x"}])
def test_bad_parameters(fixture_root, tmp_path, monkeypatch, params):
    assert run(fixture_root["linked"], tmp_path / "out", monkeypatch, **params) == 2


def test_a_bad_sample_is_skipped(fixture_root, tmp_path, monkeypatch):
    src = tmp_path / "in"
    import shutil

    shutil.copytree(fixture_root["linked"], src)
    (src / "spatial" / "cells_a.h5ad").write_bytes(b"not an h5ad")
    out = tmp_path / "out"
    assert run(src, out, monkeypatch, N_PCS="5") == 0
    prov = json.loads((out / "qc" / "provenance.json").read_text())
    assert prov["n_failed"] == 1 and prov["n_ok"] == 1
    assert not (out / "spatial" / "cells_a.h5ad").exists()
    summary = {r["sample_id"]: r for r in rows(out / "qc" / "qc_summary.csv")}
    assert summary["cells_a"]["status"].startswith("failed")


def test_no_input_fails(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1


def test_mad_thresholds():
    import pandas as pd

    rng = np.random.default_rng(0)
    counts = np.exp(rng.normal(8, 0.2, 1000))
    counts[:5] = 10
    genes = counts / 3
    mt = np.clip(rng.normal(5, 1, 1000), 0, None)
    mt[10:13] = 60
    a = ad.AnnData(obs=pd.DataFrame({"total_counts": counts, "n_genes_by_counts": genes, "pct_counts_mt": mt},
                                    index=[str(i) for i in range(1000)]))
    thr = Q.choose_thresholds(a, {"FILTER_MODE": "mad", "MADS": 5.0, "MAX_PCT_MT": 20.0, "MAX_CONTROL_FRACTION": 0.1}, "spot")
    med = np.median(np.log1p(counts))
    mad = np.median(np.abs(np.log1p(counts) - med))
    assert thr.min_counts == pytest.approx(np.expm1(med - 5 * mad))
    assert thr.max_pct_mt >= 20
    flags = Q.flag_observations(a, thr)
    assert flags["low_counts"][:5].all() and flags["high_mt"][10:13].all()
    assert list(Q.reasons_of(flags)[:1]) == ["low_counts;low_genes"]
    assert flags.any(axis=1).sum() <= 10


def test_mt_cut_applies_when_mad_is_zero():
    """A panel where nearly every observation has 0 % mt (MAD 0) still drops one above MAX_PCT_MT."""
    from processor import qc as Q
    from processor.main import read_parameters

    rng = np.random.default_rng(0)
    counts = rng.poisson(20, size=(200, 30)).astype(np.int32)
    names = [f"G{i}" for i in range(29)] + ["MT-CO1"]
    counts[:, -1] = 0
    counts[0, -1] = 500          # one observation with ~45 % mt
    a = ad.AnnData(sp.csr_matrix(counts), var=pd.DataFrame(index=names))
    Q.compute_qc_metrics(a)
    thr = Q.choose_thresholds(a, read_parameters({}), "spot")
    assert thr.max_pct_mt == 20
    flags = Q.flag_observations(a, thr)
    assert flags["high_mt"].to_numpy().nonzero()[0].tolist() == [0]


def test_empty_observations_do_not_break_the_mt_cut():
    """Cells with no counts (0/0 % mt) must not turn the mitochondrial threshold into NaN."""
    from processor import qc as Q
    from processor.main import read_parameters

    rng = np.random.default_rng(1)
    counts = rng.poisson(20, size=(100, 10)).astype(np.int32)
    counts[:3] = 0
    a = ad.AnnData(sp.csr_matrix(counts), var=pd.DataFrame(index=[f"G{i}" for i in range(10)]))
    Q.compute_qc_metrics(a)
    assert not a.obs["pct_counts_mt"].isna().any()
    thr = Q.choose_thresholds(a, read_parameters({}), "cell")
    assert thr.max_pct_mt == 20
    assert Q.flag_observations(a, thr)["low_counts"].to_numpy()[:3].all()
