# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image it builds contains the following third-party
components, installed from PyPI at the pinned versions in `requirements.txt`; their licences
apply to the image. Nothing here is GPL (this app needs no clustering, so neither igraph nor
leidenalg is installed), so the image as a whole can be redistributed under these permissive
licences (BSD, MIT, PSF, Apache-2.0, MPL-2.0; MPL-2.0 is file-level copyleft and leaves our code
untouched).

| Component | What | Licence | Source |
|---|---|---|---|
| scanpy 1.12.5 | `calculate_qc_metrics`, `filter_genes`, `normalize_total`, `log1p`, `highly_variable_genes(flavor="seurat_v3")`, `pca` | BSD-3-Clause | https://github.com/scverse/scanpy |
| anndata 0.13.4 | Reads and writes the `.h5ad` files. | BSD-3-Clause | https://github.com/scverse/anndata |
| scikit-misc 0.5.3 | The loess fit inside scanpy's seurat_v3 HVG flavour. | BSD-3-Clause | https://github.com/has2k1/scikit-misc |
| scikit-learn 1.9.1 | PCA (arpack solver on sparse input, through scanpy). | BSD-3-Clause | https://scikit-learn.org |
| NumPy 2.5.4, SciPy 1.18.1 | Arrays, CSR sparse matrices, medians and MADs. | BSD-3-Clause (NumPy also bundles 0BSD/MIT/Zlib/CC0 parts) | https://numpy.org, https://scipy.org |
| pandas 3.1.0 | `obs`/`var` tables, samples.csv. | BSD-3-Clause | https://pandas.pydata.org |
| h5py 3.16.0 | HDF5 backend of anndata. | BSD-3-Clause | https://www.h5py.org |
| Matplotlib 3.11.2, Pillow 13.0.0 | PNG figures; image arrays. | Matplotlib licence (PSF-based, BSD-compatible); MIT-CMU | https://matplotlib.org, https://python-pillow.org |
| numba 0.68.1, llvmlite 0.50.0 | Compiled kernels used by scanpy's QC metrics. | BSD-2-Clause (llvmlite also Apache-2.0 with LLVM exception) | https://numba.pydata.org |
| umap-learn 0.5.12, pynndescent 0.7.0 | Dependencies of scanpy (not used by this app). | BSD-3-Clause, BSD-2-Clause | PyPI |
| statsmodels 0.15.0, patsy 1.0.3, seaborn 0.13.2, networkx 3.7, joblib 1.6.0, threadpoolctl 3.7.0, natsort 8.4.0 | Dependencies of scanpy. | BSD-3-Clause / BSD-2-Clause (patsy) / MIT (natsort) | PyPI |
| fast-array-utils 1.5.1, legacy-api-wrap 1.5, session-info2 0.4.2, scverse-misc 0.1.8 | scverse helpers pulled in by scanpy/anndata. | MPL-2.0 (scverse-misc BSD-3-Clause) | https://github.com/scverse |
| zarr 3.4.2, numcodecs 0.17.1, donfig 0.8.1, google-crc32c 1.9.0 | Dependencies of anndata (zarr backend, unused here). | MIT (google-crc32c Apache-2.0) | PyPI |
| array-api-compat, contourpy, cycler, fonttools, kiwisolver, packaging, pyparsing, python-dateutil, six, tqdm, typing_extensions, tzdata, wrapt | Dependencies of the above (versions in `requirements.txt`). | MIT / BSD / Apache-2.0 or BSD / MPL-2.0+MIT (tqdm) / PSF-2.0 | PyPI |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

## Not bundled

- **Demo data** (not in the image or the repository): 10x Genomics *Human Lymph Node* (Visium)
  and *FFPE Human Breast Cancer* Xenium replicates 1 and 2 (CC BY 4.0; Janesick et al. 2023,
  doi:10.1038/s41467-023-43458-x); Dawo S, Nonchev K, Silina K, *10x Visium Spatial
  Transcriptomics Dataset: Kidney (3) and Lung (5) Cancer with Tertiary Lymphoid Structures*,
  Zenodo 2025, doi:10.5281/zenodo.14620362 (CC BY 4.0).
