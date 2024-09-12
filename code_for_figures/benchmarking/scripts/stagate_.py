import logging
import os
import argparse

import STAGATE_pyG_sub
import scanpy as sc
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import scipy as sp
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


def smooth(
    adata,
    k_physical,
    sample_rate,
    n_epochs,
    alpha,
    random_seed,
):

    # Get neighbors
    STAGATE_pyG_sub.Cal_Spatial_Net(
        adata,
        model='KNN',
        k_cutoff=k_physical,
        sample_rate=sample_rate,
    )

    # Train model
    adata = STAGATE_pyG_sub.train_STAGATE(
        adata,
        n_epochs=n_epochs,
        key_added='stagate',
        alpha=alpha,
        random_seed=random_seed,
    )

    return adata


def binary_search_leiden(
    adata,
    leiden_res,
    n_regions,
    jitter=0.1
):
    sc.tl.leiden(adata, resolution=leiden_res, neighbors_key='smoothed')
    n_clusters = len(adata.obs['leiden'].unique())
    lo, hi = 0, 5
    iters = 0
    while n_clusters != n_regions:
        if n_clusters < n_regions:
            lo = leiden_res
            leiden_res = (leiden_res+hi)/2
        elif n_clusters > n_regions:
            hi = leiden_res
            leiden_res = (leiden_res+lo)/2
        iters += 1
        if iters > 10:
            leiden_res += np.random.choice([-1,1]) * jitter
            hi += jitter
            lo -= max(0,jitter)
        sc.tl.leiden(adata, resolution=leiden_res, neighbors_key='smoothed')
        n_clusters = len(adata.obs['leiden'].unique())
        print(leiden_res, n_clusters, flush=True)
    print()


def cluster(
    adata,
    n_pcs,
    k_latent,
    n_regions,
    leiden_res,
):

    # Find latent neighbors
    nbrs = NearestNeighbors(n_neighbors=k_latent).fit(adata.obsm['stagate'])
    adata.obsp['smoothed_connectivities'] = nbrs.kneighbors_graph(mode='connectivity')
    adata.obsp['smoothed_distances'] = nbrs.kneighbors_graph(mode='distance')
    
    # Add metadata expected by Scanpy's Leiden and UMAP
    adata.uns['smoothed'] = {'connectivities_key': 'smoothed_connectivities',
                'distances_key': 'smoothed_distances',
                'params': {
                    'n_neighbors': k_latent,
                    'method': 'umap',
                    'random_state': 0,
                    'metric': 'euclidean',
                    'use_rep': 'stagate'
                }
    }

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=n_regions, n_init=10)
    kmeans.fit(adata.obsm['stagate'])
    adata.obs['kmeans'] = kmeans.labels_.astype(str)

    # Perform Leiden clustering
    binary_search_leiden(adata, leiden_res, n_regions)


def quantify(adata, n_regions):

    # Quantify spatial reconstruction
    Cp = adata.obsp['spatial_connectivities']
    Cl = adata.obsp['smoothed_connectivities']
    recon_score = Cp.multiply(Cl).count_nonzero() / (Cp+Cl).count_nonzero()

    # Quantify k-means accuracy
    kmeans_ari = adjusted_rand_score(adata.obs['region_true'], adata.obs['kmeans'])

    # Quantify Leiden accuracy
    leiden_ari = adjusted_rand_score(adata.obs['region_true'], adata.obs['leiden'])

    return recon_score, kmeans_ari, leiden_ari


def main(
    read_path,
    sample_rate,
    k_physical,
    k_latent,
    n_pcs,
    leiden_res,
    write_path,
    random_seed,
    decimals,
    n_epochs,
    alpha,
):

    # Initialize logger
    logger = logging.getLogger('STAGATE')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Set random seed
    np.random.seed(random_seed)

    # Load data
    logger.info(f'Loading data from {read_path}')
    adata = sc.read_h5ad(read_path)
    n_regions = len(adata.obs['region_true'].unique())

    # Smooth
    logger.info('Smoothing')
    logger.info(f'\tsample_rate={sample_rate}')
    logger.info(f'\trandom_seed={random_seed}')
    adata = smooth(adata, k_physical, sample_rate, n_epochs, alpha, random_seed)
    
    # Cluster
    logger.info('Clustering')
    cluster(adata, n_pcs, k_latent, n_regions, leiden_res)

    # Quantify results
    logger.info('Quantifying results')
    recon_score, kmeans_ari, leiden_ari = quantify(adata, n_regions)
    
    # Write results to disk
    os.makedirs(write_path, exist_ok=True)
    results = f'{np.round(recon_score, decimals=decimals)}_' + \
              f'{np.round(kmeans_ari, decimals=decimals)}_' + \
              f'{np.round(leiden_ari, decimals=decimals)}'
    filename_txt = f'{random_seed}_{results}.txt'
    write_path_txt = os.path.join(write_path, filename_txt)
    logger.info(f'Writing results to {write_path_txt}')
    f = open(write_path_txt, 'a')
    f.close()

    # Write example AnnData to disk
    if random_seed == 0:
        filename_adata = f'adata_{sample_rate}.h5ad'
        write_path_adata = os.path.join(write_path, filename_adata)
        logger.info(f'Saving AnnData to {write_path_adata}')
        adata.write(write_path_adata)

    logger.info('Done')


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--read_path', type=str)
    parser.add_argument('--sample_rate', type=float, default=None)
    parser.add_argument('--k_physical', type=int, default=50)
    parser.add_argument('--k_latent', type=int, default=15)
    parser.add_argument('--n_pcs', type=int)
    parser.add_argument('--leiden_res', type=float)
    parser.add_argument('--write_path', type=str)
    parser.add_argument('--random_seed', type=int)
    parser.add_argument('--decimals', type=int, default=5)
    parser.add_argument('--n_epochs', type=int, default=50)
    parser.add_argument('--alpha', type=float, default=0)
    args = parser.parse_args()

    main(
        args.read_path,
        args.sample_rate,
        args.k_physical,
        args.k_latent,
        args.n_pcs,
        args.leiden_res,
        args.write_path,
        args.random_seed,
        args.decimals,
        args.n_epochs,
        args.alpha
    )