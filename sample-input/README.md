# sample-input

This app reads the output of `st-prep`: `spatial/<sample_id>.h5ad` (raw counts in `X`,
`obsm["spatial"]`, `uns["pennsieve"]`) and `metadata/samples.csv`. The contract is in
`workflow-spatial-transcriptomics/docs/conventions.md`.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked/ -> data/output/{spatial,qc}/
```

Two samples: `visium_a`, a 20 × 20 hex grid of 400 Visium-like spots with a lowres "H&E" image,
300 genes (13 `MT-`), 12 spots with ~1 % of the normal counts and 8 with ~60 % mitochondrial
counts; `cells_a`, 450 Xenium-like cells on a 60-gene panel with `control_counts`, 10 cells
with a control fraction of ~0.5 and 10 cells with 2 counts.

## Real data

Run `st-prep` on any of the demo datasets of the workflow repository (`scripts/fetch-sample-data.sh`):

- 10x Genomics *Human Lymph Node* (Visium, Space Ranger 1.1.0), 4,035 spots.
- Dawo S, Nonchev K, Silina K. *10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung
  (5) Cancer with Tertiary Lymphoid Structures.* Zenodo 2025, https://doi.org/10.5281/zenodo.14620362
  (CC BY 4.0): 8 FFPE sections.
- 10x Genomics *FFPE Human Breast Cancer* Xenium replicates 1 and 2 (Janesick et al. 2023,
  https://doi.org/10.1038/s41467-023-43458-x; CC BY 4.0): 313-gene panel, ~167 k and ~118 k cells.

Then point this app at st-prep's output:

```sh
docker run --rm --platform linux/amd64 --memory 3g \
  -e INPUT_DIR=/data/input -e OUTPUT_DIR=/data/output \
  -v "/path/to/st-prep/output:/data/input:ro" -v "$PWD/data/output:/data/output" st-qc
```
