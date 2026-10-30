# st-qc

Quality control and normalization of spatial transcriptomics samples, the step after `st-prep`
in all three spatial workflows:

- **A. Spatial Domain Mapping & Condition Comparison:** `st-prep → st-qc → st-domains → {st-svg, st-compare} → st-report`
- **B. Cell-type Deconvolution & Niches:** `st-prep → st-qc → cell2location → st-niches → st-report`
- **C. Single-cell Spatial Profiling:** `st-prep → st-qc → st-celltype → {st-niches, st-ccc} → st-report`

Not every spot or cell in a spatial dataset is worth analysing. On a Visium slide, spots at a
torn or folded edge, over a bubble, or over tissue that was poorly permeabilized capture few
transcripts; spots over damaged tissue are dominated by mitochondrial transcripts. On Xenium,
MERSCOPE and CosMx, segmentation produces fragments with a handful of transcripts, and some
cells collect a large share of negative-control probe or blank-codeword signal, which says
their counts are mostly noise. Left in, such observations form their own clusters ("low
quality" domains or cell types) and distort normalization. st-qc finds them per sample,
removes them, and prepares every sample for the steps after it: it keeps the raw counts,
log-normalizes the expression, picks the highly variable genes and computes a PCA, so
`st-domains`, `st-celltype` and `cell2location` all start from the same, documented state.

For the computational analyst: per sample, streaming one `.h5ad` at a time,
`scanpy.pp.calculate_qc_metrics` (mitochondrial genes = `var["mt"]`, else `MT-`/`mt-`), the
control fraction `control_counts / (total_counts + control_counts)` for imaging panels,
outlier removal by median absolute deviations (MAD, the default, as recommended by the
scverse best-practices review) or fixed cutoffs, `scanpy.pp.filter_genes(min_cells=…)`,
`layers["counts"]` = raw counts, `X` = `log1p(normalize_total(target_sum))` as float32 CSR,
`highly_variable_genes(flavor="seurat_v3")` on the counts, and `pca` on the HVGs of `X`.
Observations are only removed, never reordered.

## Method

1. **QC metrics.** `total_counts`, `n_genes_by_counts`, `pct_counts_mt` are recomputed from
   `X` (st-prep already wrote them; recomputing makes st-qc correct on any contract-conforming
   `.h5ad`). Imaging-based samples (an `obs["control_counts"]` column from st-prep) also get
   `obs["control_fraction"]`.
2. **Observation filter** (`FILTER_MODE`):
   - `mad` (default): thresholds adapt to each sample. An observation is an outlier when
     `log1p(total_counts)` or `log1p(n_genes_by_counts)` is more than `MADS` MADs below the
     sample's median, or when `pct_counts_mt` is more than `MADS` MADs above its median **and**
     above `MAX_PCT_MT`. The log scale makes the cut symmetric for these right-skewed metrics;
     MAD-based outlier detection is what Heumos et al. (2023) recommend (Germain et al. 2020
     compared it with fixed cutoffs), and 5 MADs is the permissive value of their companion
     book (sc-best-practices.org): it removes clear outliers but not biology (a lymph-node germinal centre has fewer counts than
     the T-cell zone, and should stay). The mitochondrial rule needs both conditions so that
     a tissue with uniformly high mitochondrial content (heart, kidney, muscle) is not emptied,
     and so that a panel without mitochondrial genes (MAD 0) falls back to `MAX_PCT_MT`. When the
     MAD of counts or genes is 0 (rare; a constant column) that cut is skipped and logged.
   - `fixed`: `total_counts < MIN_COUNTS`, `n_genes_by_counts < MIN_GENES`,
     `pct_counts_mt > MAX_PCT_MT`. With `MIN_COUNTS` / `MIN_GENES` at 0 the defaults depend on
     the kind of observation (`uns["pennsieve"]["kind"]`): Visium spots 500 / 250, Visium HD
     bins 50 / 20, cells 10 / 5.
   - `none`: nothing is removed; the sample is only normalized.

   In `mad` and `fixed` mode, imaging panels also drop cells with
   `control_fraction > MAX_CONTROL_FRACTION`. Every removed observation is listed with its
   reasons in `qc/removed.csv`.
3. **Gene filter.** Genes detected in fewer than `MIN_CELLS_PER_GENE` kept observations are
   dropped (`0` keeps every gene). The QC metrics are recomputed afterwards so `obs` describes
   the matrix that is written.
4. **Normalization.** The raw counts move to `layers["counts"]` (int32 CSR when integer-valued).
   `X` = `log1p` of the counts scaled to `TARGET_SUM` per observation (`0` = the sample's median
   total counts, scanpy's default), float32 CSR. `uns["pennsieve"]["normalized"] = true`,
   `uns["pennsieve"]["target_sum"]` records the value used.
5. **Highly variable genes.** `scanpy.pp.highly_variable_genes(flavor="seurat_v3",
   n_top_genes=N_HVG, layer="counts")` (Stuart et al. 2019): a variance-stabilizing fit on the raw
   counts. Targeted panels with at most `N_HVG` genes (Xenium 313-gene, MERSCOPE, CosMx 1k)
   mark every gene as highly variable; a panel is already a selection.
6. **PCA.** `scanpy.pp.pca` on the HVGs of the log-normalized `X`, centred but not scaled
   (sparse, ARPACK, `random_state=SEED`), `min(N_PCS, n_obs − 1, n_hvg − 1)` components:
   `obsm["X_pca"]`, `varm["PCs"]`, `uns["pca"]`.

Each sample is independent: thresholds, the target sum and HVGs are per sample, and nothing
is pooled (batch integration is the job of the clustering apps).

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `spatial` | input (subfolder) | `<sample_id>.h5ad` | `st-prep` output: raw counts in `X` (CSR), `obsm["spatial"]`, `uns["pennsieve"]` (`platform`, `kind`), `uns["spatial"]` images for Visium. Found in `INPUT_DIR/spatial/`, else anywhere under `INPUT_DIR`. |
| `metadata` | input (subfolder) | `samples.csv` (+ `prep.json`, unused) | Which samples to process, and `kind` when the `.h5ad` lacks it. Without a `samples.csv` every `.h5ad` is a sample. |
| `spatial` | output | `spatial/<sample_id>.h5ad`, `spatial/provenance.json` | The kept observations (input order) and genes; contract below. |
| `qc` | output | `qc/*.csv`, `qc/<sample_id>_*.png`, `qc/summary.json`, `qc/provenance.json` | Tables and figures, below. |

### `spatial/<sample_id>.h5ad`

Everything st-prep wrote (`obs` columns, `obsm["spatial"]`, `uns["spatial"]`, `uns["pennsieve"]`)
for the kept observations and genes, plus:

| Slot | Content |
|---|---|
| `X` | log-normalized expression, float32 CSR |
| `layers["counts"]` | raw counts, int32 CSR |
| `obs["control_fraction"]` | imaging panels only |
| `var["highly_variable"]` | bool; with seurat_v3 also `highly_variable_rank`, `means`, `variances`, `variances_norm` |
| `var["n_cells"]` | observations in which the gene is detected (after filtering) |
| `obsm["X_pca"]`, `varm["PCs"]`, `uns["pca"]` | PCA (`variance`, `variance_ratio`, `params`) |
| `uns["pennsieve"]` | adds `normalized: true`, `target_sum`, `qc_version`, `qc_filter_mode` |

Written with gzip compression, like st-prep's files.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `FILTER_MODE` | `mad` \| `fixed` \| `none` | `mad` | How observations are filtered (see Method). |
| `MADS` | number | `5` | MADs from the median for `mad` mode. 3 is strict, 5 permissive. |
| `MIN_COUNTS` | integer | `0` | `fixed` mode: minimum total counts; 0 = by kind (spot 500, bin 50, cell 10). |
| `MIN_GENES` | integer | `0` | `fixed` mode: minimum genes detected; 0 = by kind (spot 250, bin 20, cell 5). |
| `MAX_PCT_MT` | number | `20` | Maximum % mitochondrial counts (`fixed`); in `mad` mode a floor under the MAD cut. |
| `MAX_CONTROL_FRACTION` | number | `0.1` | Imaging panels: maximum fraction of control-probe / blank-codeword counts. |
| `MIN_CELLS_PER_GENE` | integer | `3` | Drop genes detected in fewer kept observations; 0 keeps all. |
| `TARGET_SUM` | number | `0` | Counts per observation before `log1p`; 0 = the sample's median total counts. 100 is common for small imaging panels, 10,000 ("CP10k") for whole-transcriptome data when samples must share a scale. |
| `N_HVG` | integer | `2000` | Highly variable genes for PCA; panels with ≤ `N_HVG` genes use all. |
| `N_PCS` | integer | `50` | Principal components (capped by the sample's size). |
| `SEED` | integer | `42` | PCA random state. |
| `THREADS` | integer | `0` | BLAS / numba threads; 0 = from the container's CPU quota (2 when there is none). |

## Outputs

### `qc/qc_summary.csv`

One row per sample.

| Column | Meaning |
|---|---|
| `sample_id`, `kind`, `platform` | From `uns["pennsieve"]` / `samples.csv` |
| `n_obs_before`, `n_obs_after`, `pct_removed` | Observations in and out |
| `n_genes_before`, `n_genes_after`, `n_genes_dropped` | Genes in and out (`MIN_CELLS_PER_GENE`) |
| `median_counts`, `median_genes`, `median_pct_mt`, `median_control_fraction` | Of the kept observations |
| `filter_mode`, `min_counts`, `min_genes`, `max_pct_mt`, `max_control_fraction` | The thresholds applied (counts and genes on the linear scale; empty: no cut) |
| `n_removed_low_counts`, `n_removed_low_genes`, `n_removed_mt`, `n_removed_control` | Observations failing each test (one observation can fail several) |
| `n_hvg`, `n_pcs`, `target_sum` | What was computed |
| `status` | `ok`, or `failed: <error>` (the sample is skipped; no `.h5ad`) |

### `qc/thresholds.csv`

One row per sample × metric: `sample_id, kind, mode, metric, direction` (`below` / `above`),
`threshold` (linear scale), `median`, `mad` (on the log1p scale for counts and genes), `mads`,
`source` (the formula or the parameter it came from).

### `qc/removed.csv`

`sample_id, obs_id, reason`: every removed observation (`obs_id` = its `obs_names` value),
`reason` = `low_counts`, `low_genes`, `high_mt`, `high_control`, joined with `;` when several.

### Figures

- `qc/<sample_id>_qc.png`: histograms of total counts (log scale), genes detected, % mitochondrial
  counts and (imaging) control fraction, kept and removed stacked, with the threshold.
- `qc/<sample_id>_spatial_counts.png`: on the tissue (the lowres H&E for Visium), left the
  log total counts of every input observation, right the kept observations in grey and the
  removed ones coloured by their first reason. Removed spots that cluster at a tissue edge or
  fold are expected; removed spots spread over a whole region deserve a look.

### `qc/summary.json`, `provenance.json`

`summary.json`: totals over samples, removals by reason, and the `qc_summary.csv` rows, for
`st-report`. `provenance.json` (in both output folders): app version, tool versions,
parameters, the method as run, inputs seen, per-sample status and timings, peak memory.

## Results on the demo data

Default parameters (`FILTER_MODE=mad`, `MADS=5`), full data as written by `st-prep` 0.1.0;
outputs in `workflow-spatial-transcriptomics/data/work/<demo>/qc/` for the apps downstream.

| Demo | Sample | Kind | Observations before → after (% removed) | Genes before → after | Thresholds: counts / genes / % mt / control | Removed: low counts / few genes / control |
|---|---|---|---|---|---|---|
| lymph | lymph_node | spot | 4,035 → 3,928 (2.7 %) | 36,601 → 22,397 | 6,453 / 3,384 / 20 / – | 56 / 107 / – |
| tls | KC1 | spot | 4,173 → 4,090 (2.0 %) | 17,845 → 17,206 | 167 / 162 / 20 / – | 77 / 83 / – |
| tls | KC2 | spot | 2,166 → 2,165 (0.05 %) | 17,845 → 17,631 | 116 / 111 / 20 / – | 1 / 1 / – |
| tls | KC3 | spot | 3,217 → 3,189 (0.9 %) | 17,845 → 16,567 | 449 / 396 / 20 / – | 16 / 28 / – |
| tls | LC1 | spot | 4,364 → 4,358 (0.1 %) | 17,845 → 17,610 | 99 / 102 / 20 / – | 5 / 6 / – |
| tls | LC2 | spot | 4,456 → 4,431 (0.6 %) | 17,845 → 17,767 | 88 / 217 / 20 / – | 0 / 25 / – |
| tls | LC3 | spot | 4,178 → 4,173 (0.1 %) | 17,845 → 17,644 | 53 / 66 / 20 / – | 3 / 5 / – |
| tls | LC4 | spot | 3,009 → 3,009 (0 %) | 17,845 → 17,099 | 15 / 18 / 20 / – | 0 / 0 / – |
| tls | LC5 | spot | 2,152 → 2,082 (3.3 %) | 17,845 → 16,393 | 128 / 238 / 20 / – | 11 / 70 / – |
| xenium | breast_rep1 | cell | 167,780 → 160,674 (4.2 %) | 313 → 313 | 13.1 / 18.3 / 20 / 0.1 | 4,543 / 7,105 / 31 |
| xenium | breast_rep2 | cell | 118,752 → 116,801 (1.6 %) | 313 → 313 | 13.9 / 20.0 / 20 / 0.1 | 608 / 1,951 / 13 |

(One observation can fail several tests, so the removal columns overlap.) `xenium-20k`, st-prep's
seeded 20,000-cell subsample of each replicate, gives the same picture: 4.1 % and 1.7 % removed
with count thresholds within 1 count of the full data.

- **Lymph node** (Visium, whole transcriptome): spots with fewer than ~6,500 counts or 3,400
  genes go, 2.7 %. On the H&E they sit in a low-count patch at the left tissue margin, in two
  small clusters inside the node and as single scattered spots. Mitochondrial content is low (median 1 %), so the
  mitochondrial cut (20 %) removes nothing. 14,204 genes are detected in fewer than 3 spots.
- **TLS** (FFPE Visium, probe-based): the human probe set has **no mitochondrial genes**, so
  `pct_counts_mt` is 0 everywhere and only the count and gene cuts act. The samples differ
  more than 10-fold in depth (median 480 to 6,000 counts per spot) and in spread; the MAD thresholds follow
  each sample, from 449 counts (KC3) down to 15 (LC4, whose counts span 30 to 10,000 per spot and
  where nothing is an outlier on the log scale). 0-3 % of spots go per sample. A biologist
  comparing kidney and lung should look at those depth differences in `qc_summary.csv` first.
- **Xenium breast** (313-gene panel): cells with fewer than ~13 transcripts or ~19 genes go
  (4.2 % and 1.6 %), a strip of low-count cells along the top edge of replicate 1 and
  scattered cells (fragments and small cells) everywhere else; 1,417 cells in replicate 1 and 61 in
  replicate 2 had no counts at all. Very few cells exceed a control fraction of 0.1 (31 and 13). No genes are
  dropped, and with 313 genes ≤ `N_HVG` all genes are "highly variable".

## Runtime

`computeTypes: [standard]`, 2 vCPU / 16 GiB, 2 h timeout. Samples are processed one at a time
and nothing is pooled, so memory follows the largest sample, not the dataset.

Measured with the default parameters, amd64 image under emulation on an M-series Mac (Docker
Desktop, `--memory 3g` or `4g`, `--network none`, 2 threads), peak memory = the container
cgroup's `memory.peak`. Wall times include container start (~2 s) and varied with the other
builds running in parallel:

| Input | Observations × genes (non-zeros) | Wall time | Peak memory |
|---|---|---|---|
| lymph (1 Visium sample) | 4,035 × 36,601 (~24 M) | 16-31 s | 1.04-1.15 GB |
| tls (8 Visium samples) | 27,715 spots in all, 17,845 genes | 67 s | 1.52 GB |
| xenium (2 replicates, full) | 167,780 + 118,752 cells × 313 | 26-48 s | 1.10-1.18 GB |
| synthetic Visium HD-like bins | 200,000 × 18,000 (27.6 M) | 22 s | 1.42 GB |
| synthetic Visium HD-like bins | 400,000 × 18,000 (58.5 M) | 39 s | 2.17 GB |

Memory grows with the non-zeros of the largest sample, about 25 bytes per non-zero over a
~0.7 GB base (the matrix is held as raw int32 counts plus float32 log-normalized values, and
anndata and scanpy copy it while subsetting and selecting HVGs). A Visium HD 8 µm sample
(~500 k bins, a few hundred genes per bin: ~150-250 M non-zeros) or a Xenium 5k-panel slide
(~500 k cells, ~300 M non-zeros) should therefore peak at about 5-9 GB; 16 GiB leaves room
for denser samples, up to roughly 600 M non-zeros. 2 µm Visium HD bins (millions of bins) are
out of scope. Time is dominated by the PCA and the gzip-compressed write, about a minute per
100 M non-zeros under emulation; native amd64 is several times faster.

Image: 996 MB unpacked, 236 MB compressed (`python:3.12-slim` + scanpy and its dependencies; no
igraph, leidenalg or torch).

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python       # pytest on a synthetic fixture (no Docker)
make run PYTHON=.venv/bin/python        # build the image, run it over data/fixture/linked -> data/output
make validate PYTHON=.venv/bin/python   # app.yml against the manifest schema, CITATION.cff
```

The fixture (`tests/make_fixture.py`) has a 400-spot Visium-like sample with planted low-count
and high-mitochondrial spots and a 450-cell Xenium-like sample with planted high-control and
near-empty cells; the tests check that exactly those are removed, the AnnData contract of the
output, and both input layouts (linked subfolders, and everything flat). Real input:
`sample-input/README.md`.

```sh
docker build --platform linux/amd64 -t st-qc .
docker run --rm --platform linux/amd64 --memory 3g --network none \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "/path/to/st-prep/output:/data/input:ro" -v "$PWD/data/output:/data/output" st-qc
```

## Limitations

- **Thresholds are per sample.** MAD thresholds describe each sample's own distribution: a
  sample that is poor throughout (low counts everywhere) keeps most of its spots, and its
  medians in `qc_summary.csv` are the place to notice it. Compare `median_counts` across
  samples before trusting a condition comparison.
- **Low counts are not always low quality.** Spots over fat, fibrous stroma or the edge of a
  lumen, and small cells (lymphocytes on Xenium), have fewer counts for biological reasons. With
  5 MADs on the log scale they are rarely removed, but check the spatial figure: a removed
  region that matches a tissue structure rather than an edge or fold suggests `MADS` is too
  low or `fixed` mode is the better fit.
- **No doublet or ambient-RNA correction.** A Visium spot holds several cells by design, so
  doublet detection does not apply to spots; for imaging platforms, segmentation errors (two
  cells merged) are not detected here. Ambient RNA (SoupX, CellBender) is not removed.
- **Mitochondrial content** needs `MT-`/`mt-` genes in the panel. Xenium, MERSCOPE and CosMx
  panels and the probe sets of probe-based (FFPE) Visium usually have none; then only the count
  and gene cuts act, plus the control fraction for imaging panels.
- **Normalization is library-size scaling.** It does not correct for spatially varying
  capture efficiency or for the cell-number differences between spots (that is what
  cell2location models). Samples are not integrated or batch-corrected.
- **Size.** Designed for Visium (≤ 5k spots), Visium HD 8 µm bins (~500k bins × 18k genes) and
  Xenium / MERSCOPE / CosMx slides (~500k cells); 2 µm Visium HD bins (millions of bins) are out
  of scope for this memory size.

## Licences

Apache-2.0 for this app (see `LICENSE`). Bundled components and their licences:
`THIRD_PARTY.md` (all permissive; no GPL code in the image, since st-qc does no clustering and
installs neither igraph nor leidenalg).

## Citations

Best-practice QC (MAD outlier detection on log counts and genes) and the overall preprocessing:

> Heumos L, Schaar AC, Lance C, Litinetskaya A, Drost F, Zappia L, Lücken MD, Strobl DC, Henao J,
> Curion F, Single-cell Best Practices Consortium, Aliee H, Ansari M, Badia-i-Mompel P, Büttner M,
> Dann E, Dimitrov D, Dony L, Frishberg A, He D, Hediyeh-zadeh S, Hetzel L, Ibarra IL, Jones MG,
> Lotfollahi M, Martens LD, Müller CL, Nitzan M, Ostner J, Palla G, Patro R, Piran Z,
> Ramírez-Suástegui C, Saez-Rodriguez J, Sarkar H, Schubert B, Sikkema L, Srivastava A, Tanevski J,
> Virshup I, Weiler P, Schiller HB, Theis FJ. **Best practices for single-cell analysis across
> modalities.** Nat Rev Genet 24(8):550-572 (2023). https://doi.org/10.1038/s41576-023-00586-w,
> PMID 37002403.

MAD-based filtering evaluated against fixed cutoffs:

> Germain PL, Sonrel A, Robinson MD. **pipeComp, a general framework for the evaluation of
> computational pipelines, reveals performant single cell RNA-seq preprocessing tools.** Genome
> Biol 21(1):227 (2020). https://doi.org/10.1186/s13059-020-02136-7, PMID 32873325.

The seurat_v3 highly-variable-gene method:

> Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, Hao Y, Stoeckius M,
> Smibert P, Satija R. **Comprehensive Integration of Single-Cell Data.** Cell
> 177(7):1888-1902.e21 (2019). https://doi.org/10.1016/j.cell.2019.05.031, PMID 31178118.

Scanpy and anndata:

> Wolf FA, Angerer P, Theis FJ. **SCANPY: large-scale single-cell gene expression data
> analysis.** Genome Biol 19(1):15 (2018). https://doi.org/10.1186/s13059-017-1382-0,
> PMID 29409532.

> Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. **anndata: Access and store annotated data
> matrices.** J Open Source Softw 9(101):4371 (2024). https://doi.org/10.21105/joss.04371
> (not in PubMed; checked on Crossref).

Demo data:

> 10x Genomics. *Human Lymph Node* (Visium, Space Ranger 1.1.0) and *FFPE Human Breast Cancer*
> Xenium replicates 1 and 2, CC BY 4.0, from the 10x Genomics dataset pages (no DOI). The
> Xenium data are described in: Janesick A, Shelansky R, Gottscho AD, Wagner F, Williams SR,
> Rouault M, Beliakoff G, Morrison CA, Oliveira MF, Sicherman JT, Kohlway A, Abousoud J,
> Drennon TY, Mohabbat SH, 10x Development Teams, Taylor SEB. **High resolution mapping of the
> tumor microenvironment using integrated single-cell, spatial and in situ analysis.** Nat
> Commun 14(1):8353 (2023). https://doi.org/10.1038/s41467-023-43458-x, PMID 38114474.

> Dawo S, Nonchev K, Silina K. **10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung
> (5) Cancer with Tertiary Lymphoid Structures.** Zenodo (2025).
> https://doi.org/10.5281/zenodo.14620362 (CC BY 4.0; checked on DataCite).

Every DOI above was checked on Crossref (or DataCite) and every PMID on PubMed (E-utilities
`esummary`: title, author list and year match). All are in `CITATION.cff` with their full
author lists.
