"""
Command-line interface for using SPIN from the shell.
"""

from __future__ import annotations

import argparse

from .spin import spin


def spin_cli():
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='SPatially INtegrate and cluster one or more spatially resolved \
                        transcriptomics datasets'
    )
    parser.add_argument('--adata_paths', type=str, nargs='+', default=None)
    parser.add_argument('--write_path', type=str, default=None)
    parser.add_argument('--batch_key', type=str, default=None)
    parser.add_argument('--batch_labels', type=str, nargs='+', default=None)
    parser.add_argument('--n_nbrs', type=int, nargs='+', default=30)
    parser.add_argument('--n_samples', type=int, nargs='+', default=None)
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

    # Run SPIN
    spin(
        adata_paths=args.adata_paths,
        write_path=args.write_path,
        batch_key=args.batch_key,
        batch_labels=args.batch_labels,
        n_nbrs=args.n_nbrs,
        n_samples=args.n_samples,
        spatial_key=args.spatial_key,
        n_pcs=args.n_pcs,    
        svd_solver=args.svd_solver,
        pca_key=args.pca_key,
        region_key=args.region_key,
        umap_key=args.umap_key,
        resolution=args.resolution,
        verbose=args.verbose,
        random_state=args.random_state,
    )

#main(args.adata_paths, args.write_path, args.batch_key, args.batch_labels,
        #args.n_nbrs, args.n_samples, args.spatial_key, args.n_pcs, args.svd_solver,
        #args.pca_key, args.region_key, args.umap_key, args.resolution, args.verbose,
        #args.random_state)
