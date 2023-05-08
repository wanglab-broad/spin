# SPIN: <ins>sp</ins>atial <ins>in</ins>tegration of spatially-resolved transcriptomics (SRT) data

## Installation

### From GitHub:
```
pip install git+https://github.com/wanglab-broad/spin@main
```
Installation should complete within ~2 mins.

## Requirements:

### Software:
Tested on MacOS (Monterey, Ventura) and Linux (Red Hat Enterprise Linux 7).

For Python package dependencies, see `pyproject.toml`.

### Data:
Requires one or more SRT datasets in `.h5ad` format, each including an expression matrix under `.X` and spatial coordinates under `.obsm`.

## Usage
### In Python:
Consider the marmoset and mouse data we provide as a demo:
```python
adata = sc.read_h5ad(
    'data/demo.h5ad',
     backup_url='https://zenodo.org/record/TEMP/files/demo.h5ad?download=1'
)
```
This dataset contains expression and spatial data from marmoset and mouse brains, corresponding to the labels `marmoset` and `mouse` under the key `.obs['species']`.

To spatially integrate this data, the single dataset can be passed into `spin.integrate` while specifying the batch key along with the unique batch labels:
```python
adata = spin.integrate(
    adata,
    batch_key='species',
    batch_labels=['marmoset', 'mouse'],
    n_nbrs=30,
    n_samples=12,
)
```
In short, this spatially smooths each dataset individually, performs PCA jointly, integrates the resulting PCs using Harmony, and stores the output under `adata.obsm['X_pca_spin']`.

Alternatively, one can provide multiple datasets corresponding to each batch:
```python
adata = spin.integrate(
    [adata_marmoset, adata_mouse],
    batch_key='species',
    batch_labels=['marmoset', 'mouse'],
    n_nbrs=30,
    n_samples=12,
)
```

... or a single dataset for finding regions in just one sample:
```python
adata = spin.integrate(
    adata_marmoset,
    n_nbrs=30,
    n_samples=12,
)
```

After integration, the data can be clustered using `spin.cluster`, which performs UMAP and Leiden on the integrated expression PCs:
```python
adata = spin.cluster(
    adata,
    resolution=0.5
)
```
Region clusters are stored under `.obs['region']` and UMAP coordinates under `.obsm['X_umap_spin']`.

The resulting region clusters can then be visualized using standard Scanpy functions:
```python
# In physical space
sc.set_figure_params(figsize=(8,5))
sc.pl.embedding(adata, basis='spatial', color='region', s=7)

# In UMAP space
sc.set_figure_params(figsize=(5,5))
sc.pl.embedding(adata, basis='X_umap_spin', color='region', s=3)
```

For details on each parameter, run `help(spin)` after importing SPIN into Python.

For details on downstream analysis (e.g. DEG analysis, trajectory inference), please see the tutorial.

### From the shell:
Requires a read path to the relevant dataset(s) as well as a write path for the output dataset. Otherwise, provide the same parameters you would when running in Python (e.g. above):
```python
python spin/src/spin.py \
--adata_paths data/demo.h5ad \
--write_path data/demo_integrated.h5ad \
--batch_key species \
--batch_labels marmoset mouse \
--n_nbrs 30 \
--n_samples 12 \
--resolution "0.5"
```

Just as when running in Python, multiple datasets can be passed in instead:
```python
python spin/src/spin.py \
--adata_paths data/demo_marmoset.h5ad data/demo_mouse.h5ad \
--write_path data/demo_integrated.h5ad \
--batch_key species \
--batch_labels marmoset mouse \
--n_nbrs 30 \
--n_samples 12 \
--resolution "0.5"
```

... or just a single dataset:
```python
python spin/src/spin.py \
--adata_paths data/demo_marmoset.h5ad\
--write_path data/demo_integrated.h5ad \
--n_nbrs 30 \
--n_samples 12 \
--resolution "0.5"
```
