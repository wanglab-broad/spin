# SPIN: <ins>sp</ins>atial <ins>in</ins>tegration of spatially-resolved transcriptomics (SRT) data
[![Biorxiv badge](https://zenodo.org/badge/doi/TEMP)](https://doi.org/TEMP) ⬅️ read the preprint here <br>
[![Zenodo badge](https://zenodo.org/badge/doi/TEMP)](https://doi.org/TEMP) ⬅️ access the data here <br>

This package is an implementation of the approach described in the manuscript linked above. In short, conventional single-cell analysis can identify molecular *cell types* by considering each cell individually, without spatial information.
![alt text]

Arguably the simplest way to incorporate spatial information and identify molecular *tissue regions* is to spatially smooth gene expression features across cells in the tissue. This can be done by setting the features of each cell to the average of its spatial neighborhood.
![alt text]

However, a problem arises when smoothed representations of each cell are compared to one another. Physically adjacent cells will have almost identical neighborhoods and thus identical smoothed representations.
![alt text]

Because conventional methods for downstream anlaysis rely on the nearest neighbors graph in feature space, we run into a problem: nearest neighbors in feature space are just nearest neighbors in physical space. This leads to reconstruction of physical space in latent space rather than identifying the desired large-scale molecular patterns.

Here, we implement an approach in which each cell's spatial neighborhood is randomly subsampled before averaging, allowing the *exact neighborhood* composition to vary while still maintaining the *general molecular* composition.
![alt text]

Ultimately, this approach enables the application of conventional single-cell tools to spatial molecular features, yielding regional analogies for each tool. For more details, please refer to the manuscript.

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
This dataset contains expression and spatial data from marmoset and mouse brains, corresponding to the cell labels `'marmoset'` and `'mouse'` under the key `.obs['species']`.

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
In short, this performs spatial subsampling and smoothing to each dataset individually (details in manuscript linked above), performs PCA jointly, integrates the resulting PCs using Harmony, and stores the output under `adata.obsm['X_pca_spin']`.

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

After integration, the data can then be clustered using `spin.cluster`, which performs UMAP and Leiden on the integrated expression PCs:
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
Downstream analysis (e.g. DEG analysis, trajectory inference) can then be performed using standard Scanpy functions as well.

For details on the parameters of `spin.integrate` and `spin.cluster`, run `help(spin)` after importing SPIN into Python.

For details on downstream analysis, please see the [tutorial](tutorials/tutorial.ipynb) notebook.

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
