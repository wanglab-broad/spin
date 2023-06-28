"""
SPatially INtegrate spatially resolved transcriptomics (SRT) datasets.
"""

from __future__ import annotations

import logging
from typing import Optional, Collection

import scanpy as sc
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from anndata import AnnData


# Create logger
logger = logging.getLogger('SPIN')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


def spin(
    adatas: Optional[Collection[AnnData] | AnnData] = None,
    adata_paths: Optional[Collection[str]] = None,
    write_path: Optional[str] = None,
    batch_key: Optional[str] = None,
    batch_labels: Optional[Collection[str]] = None,
    n_nbrs: Collection[int] | int = 30,
    n_samples: Collection[Optional[int]] | Optional[int] = None,
    spatial_key: str = 'spatial',
    n_pcs: int = 50,
    svd_solver: str = 'randomized',
    pca_key: str = 'X_pca_spin',
    region_key: str = 'region',
    umap_key: str = 'X_umap_spin',
    resolution: float = 0.5,
    verbose: bool = True,
    random_state: int = 0,
):
    """\
    Spatially integrate and cluster SRT data using SPIN from the shell.
    
    Parameters
    ----------
    adatas
        One or more SRT datasets.
        Assumed to be have been normalized prior.
    adata_paths
        Paths to one or more SRT datasets.
    write_path
        Path to write integrated data to.
    batch_key
        The key to batch information within `adata.obs`.
    batch_labels
        Labels corresponding to each batch. Relevant when integrating across multiple
        `adatas`. Will be stored under `adata.obs[<batch_key>]`.
    n_nbrs
        Number of nearest neighbors to find for each cell.
    n_samples
        Number of random neighbor samples used for averaging.
    spatial_key
        The key to spatial coordinates within `adata.obsm`.
    n_pcs
        Number of PCs to calculate for dimension reduction.
    svd_solver
        SVD solver to use.
    pca_key
        The key to store PCA output under within `adata.obsm`.
    region_key
        The key to store region labels under within `adata.obs`.
    umap_key
        The key to store UMAP output under within `adata.obsm`.
    resolution
        Resolution for Leiden clustering
    verbose
        Display updates on function progress.
    random_state
        Random seed used for smoothing, PCA, and Harmony.
    """
    adata_source_types = np.array([type(sc.AnnData()), type(None)])
    n_adata_sources = np.sum(a==type(None))
    assert n_adata_sources == 1, "Requires either a list of paths OR a list of AnnDatas"
    
    # Read data
    if adata_paths:
        adatas = []
        for path in adata_paths:
            if verbose:
                logger.info(f'Reading {path}')
            adatas.append(sc.read_h5ad(path))

    # Integrate spatial features across samples
    adata = _integrate(
        adatas,
        batch_key=batch_key,
        batch_labels=batch_labels,
        n_nbrs=n_nbrs,
        n_samples=n_samples,
        spatial_key=spatial_key,
        n_pcs=n_pcs,    
        svd_solver=svd_solver,
        pca_key=pca_key,
        random_state=random_state,
        verbose=verbose,
    )

    # Cluster integrated samples to find regions
    adata = _cluster(
        adata,
        pca_key=pca_key,
        region_key=region_key,
        umap_key=umap_key,
        resolution=resolution,
        verbose=verbose,
    )

    # Write data
    if write_path:
        if verbose:
            adata.write(write_path)
            logger.info(f'Written to {write_path}')
    else:
        return adata


def _integrate(
    adatas: Collection[AnnData] | AnnData,
    batch_key: Optional[str] = None,
    batch_labels: Optional[Collection[str]] = None,
    n_nbrs: Collection[int] | int = 30,
    n_samples: Collection[Optional[int]] | Optional[int] = None,
    spatial_key: str = 'spatial',
    n_pcs: int = 50,
    svd_solver: str = 'randomized',
    pca_key: str = 'X_pca_spin',
    random_state: int = 0,
    verbose: bool = True,
) -> AnnData:
    """\
    Smooth and integrate SRT datasets.
    
    Parameters
    ----------
    adatas
        One or more SRT datasets.
        Assumed to be have been normalized prior.
    batch_key
        The key to batch information within `adata.obs`. Relevant when integrating across
        multiple samples.
        If only analyzing a single sample, leave as `None`.
    batch_labels
        Labels corresponding to each adata. Will be stored under `adata.obs[<batch_key>]`.
        Required if passing in multiple adatas. Otherwise, leave as `None`.
    n_nbrs
        Number of nearest neighbors to find for each cell.
        Either provide n_nbrs for each dataset or a single n_nbrs for all datasets.
        Needs to be the same length as n_samples.
    n_samples
        Number of random neighbor samples used for averaging.
        Either provide n_samples for each dataset or a single n_samples for all datasets.
        Needs to be the same length as n_nbrs.
    spatial_key
        The key to spatial coordinates within `adata.obsm`.
    n_pcs
        Number of PCs to calculate for dimension reduction.
    svd_solver
        SVD solver to use.
    pca_key
        The key to store PCA output under within `adata.obsm`.
    random_state
        Random seed used for smoothing, PCA, and Harmony.
    verbose
        Display updates on function progress.

    Returns
    -------
    Copy of adata input containing integrated spatial expression PCs.
    """
    # Handle non-Collection input (e.g. single AnnData, n_nbrs, and/or n_samples)
    if type(adatas) == AnnData:
        adatas = [adatas]
    if type(n_nbrs) == int:
        n_nbrs = [n_nbrs]
    if (type(n_samples) == int) or (n_samples == None):
        n_samples = [n_samples]

    # Split single AnnData by batch, if batch_key provided
    n_adatas = len(adatas)
    if (n_adatas == 1):
        if batch_key:
            if verbose:
                logger.info(f'Splitting adata by `{batch_key}`')
            batch_labels = adatas[0].obs[batch_key].unique()
            adatas = [adatas[0][adatas[0].obs[batch_key]==batch] for batch in batch_labels]
            n_adatas = len(adatas)

    # If single n_nbrs/n_samples for multiple batches, clone
    if len(n_nbrs) < n_adatas:
        n_nbrs *= n_adatas
    if len(n_samples) < n_adatas:
        n_samples *= n_adatas

    # Smooth each batch independently
    if verbose:
        logger.info('Smoothing')
    for i in range(n_adatas):
        _get_nbrs(
            adatas[i],
            n_nbrs[i],
            spatial_key,
        )
        _smooth(
            adatas[i],
            n_samples[i],
            random_state,
        )

    # Concatenate batches and add metadata
    if n_adatas > 1:
        adata = sc.concat(adatas, keys=batch_labels, label=batch_key, join='inner')
        adata.uns['n_nbrs'] = n_nbrs
        adata.uns['n_samples'] = n_samples
    else:
        adata = adatas[0]

    # Run PCA
    if verbose:
        logger.info('Performing PCA')
    pca = PCA(n_components=n_pcs, svd_solver=svd_solver, random_state=random_state)
    adata.obsm[pca_key] = pca.fit_transform(adata.layers['smooth'])
    adata.varm[pca_key] = pca.components_.T

    # Integrate PCs
    if batch_labels:
        if verbose:
            logger.info('Integrating')
        sc.external.pp.harmony_integrate(
            adata,
            batch_key,
            basis=pca_key,
            adjusted_basis=pca_key,
            random_state=random_state,
            verbose=verbose,
        )
        if verbose:
            logger.info('Integration complete')

    return adata.copy() # copying necessary for multiple runs on single AnnData


def _cluster(
    adata: AnnData,
    pca_key: str = 'X_pca_spin',
    region_key: str = 'region',
    umap_key: str = 'X_umap_spin',
    resolution: float = 0.5,
    verbose: bool = True,
) -> AnnData:
    """\
    Create nearest neighbors graph in latent space and perform UMAP and Leiden.
    
    Parameters
    --------
    adata
        SRT dataset
    pca_key
        The key to PCA output within `adata.obsm`.
    region_key
        The key to store region labels under within `adata.obs`.
    umap_key
        The key to store UMAP output under within `adata.obsm`.
    resolution
        Resolution for Leiden clustering.
    verbose
        Display updates on function progress.

    Returns
    --------
    Copy of adata input containing region cluster labels and UMAP coordinates
    """
    if verbose:
        logger.info('Finding latent neighbors')
    sc.pp.neighbors(
        adata,
        use_rep=pca_key,
        key_added=region_key,
    )
    if verbose:
        logger.info('Leiden clustering')
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=region_key,
        neighbors_key=region_key,
    )
    if umap_key:
        if verbose:
            logger.info('Performing UMAP')
        umap = sc.tl.umap(
            adata,
            neighbors_key=region_key,
            copy=True,
        ).obsm['X_umap']
        adata.obsm[umap_key] = umap
    if verbose:
        logger.info('Clustering complete')
    
    return adata.copy()


def _get_nbrs(
    adata: AnnData,
    n_nbrs: int,
    spatial_key: str,
) -> np.ndarray:
    """\
    Find spatial nearest neighbors of each cell.
    
    Parameters
    ----------
    adata
        SRT dataset.
    n_nbrs
        Number of nearest neighbors to find for each cell.
    spatial_key
        The key to spatial coordinates within `adata.obsm`.

    Returns
    -------
    Updates `adata` with the following fields:

    `.obsm['nbr_idxs']`:
        Matrix of neighbor indices of shape `n_obs` × `n_nbrs`
    `.uns['n_nbrs']`:
        Number of neighbors used
    `.uns['spatial_key']`:
        The key to `.obsm` for the domain used to find neighbors  
    """
    # Find spatial nearest neighbors
    coordinates = adata.obsm[spatial_key]
    nbrs = NearestNeighbors(n_neighbors=n_nbrs)
    nbrs.fit(coordinates)
    _, nbr_idxs = nbrs.kneighbors(coordinates)

    # Save metadata
    adata.obsm['nbr_idxs'] = nbr_idxs
    adata.uns['n_nbrs'] = n_nbrs
    adata.uns['spatial_key'] = spatial_key


def _smooth(
    adata: AnnData,
    n_samples: Optional[int],
    random_state: int,
) -> np.ndarray:
    """\
    Set each cell's representation to the average of its randomly-subsampled neighborhood.
    
    Parameters
    ----------
    adata
        SRT dataset.
    nbr_idxs
        Matrix of neighbor indices of shape `n_obs` × `n_nbrs`. Rows correspond to cells
        and columns to the adata indices of their nearest neighbors.
    n_samples
        Number of random neighbor samples used for averaging.
        If None, sample 1/3 of n_nbrs, which is inferred from `nbr_idxs`.
    random_state
        Random seed for randomly subsampling neighborhoods.

    Returns
    -------
    Updates `adata.layers` with the following fields:

    `smooth`:
        Smoothed data matrix of shape `n_obs` × `n_vars`
    `n_samples`:
        Number of random neighbor samples used for averaging
    """
    # Randomly subsample each cell's neighborhood
    if not n_samples:
        n_samples = adata.obsm['nbr_idxs'].shape[1] // 3
    np.random.seed(random_state)
    nbr_idxs_sampled = np.array([
        np.random.choice(idxs, size=n_samples, replace=False)
        for idxs in adata.obsm['nbr_idxs']
    ])

    # Set each cell's representation to the average of its subsampled neighborhood
    X_smooth = np.zeros(adata.X.shape)
    for nth_nbrs in np.array(nbr_idxs_sampled).T:
        X_smooth += adata.X[nth_nbrs] / n_samples

    # Save metadata
    adata.obsm['nbr_idxs_sampled'] = nbr_idxs_sampled
    adata.layers['smooth'] = np.array(X_smooth) # in case adata.X is sparse
    adata.uns['n_samples'] = n_samples
