import logging
import os
import argparse

from GraphST import GraphST_mod
import torch
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
    learning_rate,
    alpha,
    beta,
    gamma,
    random_seed,
):

    # Find latent neighbors
    nbrs = NearestNeighbors(n_neighbors=k_physical).fit(adata.obsm['spatial'])
    A = nbrs.kneighbors_graph()
    adata.obsp['spatial_connectivities'] = A.copy()

    # Subsample
    if sample_rate:
        adata.obsp['spatial_connectivities_subsampled'] = A.copy()
        for i in range(len(adata)):
            nbr_idxs = adata.obsp['spatial_connectivities'][i].nonzero()[1]
            idxs_to_drop = np.random.choice(
                nbr_idxs,
                size=int(len(nbr_idxs)*(1-sample_rate)),
                replace=False,
            )
            for j in idxs_to_drop:
                adata.obsp['spatial_connectivities_subsampled'][i,j] = 0
        adata.obsp['spatial_connectivities_subsampled'].eliminate_zeros()
    else:
        adata.obsp['spatial_connectivities_subsampled'] = np.eye(adata.shape[0])

    # Store in GraphST-readable format
    adata.obsm['adj'] = adata.obsp['spatial_connectivities_subsampled']
    adata.obsm['graph_neigh'] = adata.obsp['spatial_connectivities_subsampled']

    # Avoid preprocessing (method consistency argument; choosing simplicity)
    adata.var['highly_variable'] = np.ones(adata.shape[1], dtype=bool)

    # Train model
    device = torch.device('cpu')
    model = GraphST_mod.GraphST(
        adata,
        datatype='Stereo',
        device=device,
        epochs=n_epochs,
        learning_rate=learning_rate,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )
    adata = model.train()
    adata.obsm['X_smoothed'] = adata.obsm['emb'].copy()

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

    # Calculate latent representation
    pca = PCA(n_components=n_pcs)
    adata.obsm['X_smoothed_pca'] = pca.fit_transform(adata.obsm['X_smoothed'])
    
    # Find latent neighbors
    nbrs = NearestNeighbors(n_neighbors=k_latent).fit(adata.obsm['X_smoothed_pca'])
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
                    'use_rep': 'X_smoothed_pca'
                }
    }

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=n_regions, n_init=10)
    kmeans.fit(adata.obsm['X_smoothed_pca'])
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
    learning_rate,
    alpha,
    beta,
    gamma,
):

    # Initialize logger
    logger = logging.getLogger('GraphST')
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
    adata = smooth(
        adata,
        k_physical,
        sample_rate,
        n_epochs,
        learning_rate,
        alpha,
        beta,
        gamma,
        random_seed,
    )
    
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
    parser.add_argument('--n_epochs', type=int, default=100)
    parser.add_argument('--learning_rate', type=float, default=0.001)
    parser.add_argument('--alpha', type=float, default=10)
    parser.add_argument('--beta', type=float, default=1)
    parser.add_argument('--gamma', type=float, default=1)
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
        args.learning_rate,
        args.alpha,
        args.beta,
        args.gamma,
    )