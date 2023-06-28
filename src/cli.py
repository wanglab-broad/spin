from __future__ import annotations

import logging
from typing import Optional, Collection
import argparse

from anndata import AnnData

from spin import integrate, cluster


def main(
    adata_paths: Collection[str],
    write_path: str,
    batch_key: Optional[str],
    batch_labels: Optional[Collection[str]],
    n_nbrs: Collection[int],
    n_samples: Collection[Optional[int]],
    spatial_key: str,
    n_pcs: int,
    svd_solver: str,
    pca_key: str,
    region_key: str,
    umap_key: str,
    resolution: int,
    verbose: bool,
    random_state: int,
):
    """
    Spatially integrate and cluster SRT data using SPIN from the shell.
    
    Parameters
    ----------
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
    # Read data
    adatas = []
    for path in adata_paths:
        if verbose:
            logger.info(f'Reading {path}')
        adatas.append(sc.read_h5ad(path))

    # Integrate spatial features across samples
    adata = integrate(
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
    adata = cluster(
        adata,
        pca_key=pca_key,
        region_key=region_key,
        umap_key=umap_key,
        resolution=resolution,
        verbose=verbose,
    )

    # Write data
    if verbose:
        adata.write(write_path)
        logger.info(f'Written to {write_path}')


if __name__=='__main__':

    parser = argparse.ArgumentParser(
        description='SPatially INtegrate and cluster one or more spatially resolved \
                     transcriptomics datasets'
    )
    parser.add_argument('--adata_paths', type=str, nargs='+')
    parser.add_argument('--write_path', type=str)
    parser.add_argument('--batch_key', type=str, default=None)
    parser.add_argument('--batch_labels', type=str, nargs='+', default=None)
    parser.add_argument('--n_nbrs', type=int, nargs='+', default=[30])
    parser.add_argument('--n_samples', type=int, nargs='+', default=[None])
    parser.add_argument('--spatial_key', type=str, default='spatial')
    parser.add_argument('--n_pcs', type=int, default=50)
    parser.add_argument('--svd_solver', default='randomized')
    parser.add_argument('--pca_key', default='X_pca_spin')
    parser.add_argument('--region_key', default='region')
    parser.add_argument('--umap_key', default='X_umap_spin')
    parser.add_argument('--resolution', type=float, default=0.5)
    parser.add_argument('--verbose', type=bool, default=True)
    parser.add_argument('--random_state', type=int, default=0)
    args = parser.parse_args()

    main(args.adata_paths, args.write_path, args.batch_key, args.batch_labels,
         args.n_nbrs, args.n_samples, args.spatial_key, args.n_pcs, args.svd_solver,
         args.pca_key, args.region_key, args.umap_key, args.resolution, args.verbose,
         args.random_state)
