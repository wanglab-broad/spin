"""
Command-line interface for using SPIN from the shell.
"""

from __future__ import annotations

import argparse

from .spin import spin


def main(
    adata_paths,
    write_path,
    batch_key,
    batch_labels,
    n_nbrs,
    n_samples,
    spatial_key,
    n_pcs,
    svd_solver,
    pca_key,
    region_key,
    umap_key,
    resolution,
    verbose,
    random_state,
):
    spin(
        adata_paths=adata_paths,
        write_path=write_path,
        batch_key=batch_key,
        batch_labels=batch_labels,
        n_nbrs=n_nbrs,
        n_samples=n_samples,
        spatial_key=spatial_key,
        n_pcs=n_pcs,    
        svd_solver=svd_solver,
        pca_key=pca_key,
        region_key=region_key,
        umap_key=umap_key,
        resolution=resolution,
        verbose=verbose,
        random_state=random_state,
    )


# Parse arguments
parser = argparse.ArgumentParser(
    description='SPatially INtegrate and cluster one or more spatially resolved \
                    transcriptomics datasets'
)
parser.add_argument('--adata_paths', type=str, nargs='+', default=None)
parser.add_argument('--write_path', type=str, default=None)
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
