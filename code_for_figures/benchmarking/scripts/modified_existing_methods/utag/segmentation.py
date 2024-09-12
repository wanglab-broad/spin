import typing as tp
import warnings
import os

import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import anndata
import parmap

from utag.types import Path, Array, AnnData
from utag.utils import sparse_matrix_dstack


def utag_mod(
    adata: AnnData,
    random_state: int,
    channels_to_use: tp.Sequence[str] = None,
    slide_key: tp.Optional[str] = "Slide",
    save_key: str = "UTAG Label",
    filter_by_variance: bool = False,
    max_dist: float = 0.1,
    normalization_mode: str = "l1_norm",
    keep_spatial_connectivity: bool = False,
    pca_kwargs: tp.Dict[str, tp.Any] = dict(n_comps=10),
    apply_umap: bool = False,
    umap_kwargs: tp.Dict[str, tp.Any] = dict(),
    apply_clustering: bool = True,
    clustering_method: tp.Sequence[str] = ["leiden", "kmeans"],
    resolution: float = 0.5,
    n_regions: int = 4,
    leiden_kwargs: tp.Dict[str, tp.Any] = None,
    parallel: bool = True,
    processes: int = None,
    sample_rate: float = 1.0
) -> AnnData:
    """
    Discover tissue architechture in single-cell imaging data
    by combining phenotypes and positional information of cells.

    Parameters
    ----------
    adata: AnnData
        AnnData object with spatial positioning of cells in obsm 'spatial' slot.
    channels_to_use: Optional[Sequence[str]]
        An optional sequence of strings used to subset variables to use.
        Default (None) is to use all variables.
    max_dist: float
        Maximum distance to cut edges within a graph.
        Should be adjusted depending on resolution of images.
        For imaging mass cytometry, where resolution is 1um, 20 often gives good results.
        Default is 20.
    slide_key: {str, None}
        Key of adata.obs containing information on the batch structure of the data.
        In general, for image data this will often be a variable indicating the image
        so image-specific effects are removed from data.
        Default is "Slide".
    save_key: str
        Key to be added to adata object holding the UTAG clusters.
        Depending on the values of `clustering_method` and `resolutions`,
        the final keys will be of the form: {save_key}_{method}_{resolution}".
        Default is "UTAG Label".
    filter_by_variance: bool
        Whether to filter vairiables by variance.
        Default is False, which keeps all variables.
    max_dist: float
        Recommended values are between 20 to 50 depending on magnification.
        Default is 20.
    normalization_mode: str
        Method to normalize adjacency matrix.
        Default is "l1_norm", any other value will not use normalization.
    keep_spatial_connectivity: bool
        Whether to keep sparse matrices of spatial connectivity and distance in the obsp attribute of the
        resulting anndata object. This could be useful in downstream applications.
        Default is not to (False).
    pca_kwargs: Dict[str, Any]
        Keyword arguments to be passed to scanpy.pp.pca for dimensionality reduction after message passing.
        Default is to pass n_comps=10, which uses 10 Principal Components.
    apply_umap: bool
        Whether to build a UMAP representation after message passing.
        Default is False.
    umap_kwargs: Dict[str, Any]
        Keyword arguments to be passed to scanpy.tl.umap for dimensionality reduction after message passing.
        Default is 10.0.
    apply_clustering: bool
        Whether to cluster the message passed matrix.
        Default is True.
    resolution: float
        What resolution should the methods in `clustering_method` be run at.
    leiden_kwargs: dict[str, Any]
        Keyword arguments to pass to scanpy.tl.leiden.
    parallel: bool
        Whether to run message passing part of algorithm in parallel.
        Will accelerate the process but consume more memory.
        Default is True.
    processes: int
        Number of processes to use in parallel.
        Default is to use all available (-1).

    Returns
    -------
    adata: AnnData
        AnnData object with UTAG domain predictions for each cell in adata.obs, column `save_key`.
    """
    ad = adata.copy()
    np.random.seed(random_state)

    if channels_to_use:
        ad = ad[:, channels_to_use]

    if filter_by_variance:
        ad = low_variance_filter(ad)

    if isinstance(clustering_method, list):
        clustering_method = [m.upper() for m in clustering_method]
    elif isinstance(clustering_method, str):
        clustering_method = [clustering_method.upper()]
    else:
        print(
            "Invalid Clustering Method. Clustering Method Should Either be a string or a list"
        )
        return
    assert all(m in ["LEIDEN", "PARC", "KMEANS"] for m in clustering_method)

    if "PARC" in clustering_method:
        from parc import PARC  # early fail if not available
    if "KMEANS" in clustering_method:
        from sklearn.cluster import KMeans

    print("Applying UTAG Algorithm...", flush=True)
    if slide_key:
        ads = [
            ad[ad.obs[slide_key] == slide].copy() for slide in ad.obs[slide_key].unique()
        ]
        ad_list = parmap.map(
            _parallel_message_pass,
            ads,
            radius=max_dist,

            coord_type="generic",
            set_diag=True,
            mode=normalization_mode,
            pm_pbar=True,
            pm_parallel=parallel,
            pm_processes=processes,
        )
        ad_result = anndata.concat(ad_list)
        if keep_spatial_connectivity:
            ad_result.obsp["spatial_connectivities"] = sparse_matrix_dstack(
                [x.obsp["spatial_connectivities"] for x in ad_list]
            )
            ad_result.obsp["spatial_distances"] = sparse_matrix_dstack(
                [x.obsp["spatial_distances"] for x in ad_list]
            )
    else:
        sq.gr.spatial_neighbors(ad, radius=max_dist, coord_type="generic", set_diag=True)
        ad_result = custom_message_passing(ad, mode=normalization_mode, sample_rate=sample_rate)

    if apply_clustering:
        if "n_comps" in pca_kwargs:
            if pca_kwargs["n_comps"] > ad_result.shape[1]:
                pca_kwargs["n_comps"] = ad_result.shape[1] - 1
                print(
                    f"Overwriding provided number of PCA dimensions to match number of features: {pca_kwargs['n_comps']}"
                )
        pca_kwargs["n_comps"] = int(pca_kwargs["n_comps"])
        sc.tl.pca(ad_result, **pca_kwargs)
        sc.pp.neighbors(ad_result)

        if apply_umap:
            print("Running UMAP on Input Dataset...", flush=True)
            sc.tl.umap(ad_result, **umap_kwargs)

        res_key1 = "leiden"
        res_key3 = "kmeans"
        
        if "LEIDEN" in clustering_method:
            
            def binary_search_leiden(adata, leiden_res, n_regions, kwargs):
                sc.tl.leiden(adata, resolution=leiden_res)
                n_clusters = len(adata.obs['leiden'].unique())
                lo, hi = 0, 5
                while n_clusters != n_regions:
                    if n_clusters < n_regions:
                        lo = leiden_res
                        leiden_res = (leiden_res+hi)/2
                    elif n_clusters > n_regions:
                        hi = leiden_res
                        leiden_res = (leiden_res+lo)/2
                    sc.tl.leiden(adata, resolution=leiden_res)
                    n_clusters = len(adata.obs['leiden'].unique())
                    print(leiden_res, n_clusters, flush=True)

            print(f"Applying Leiden Clustering...", flush=True)
            kwargs = dict()
            kwargs.update(leiden_kwargs or {})

            binary_search_leiden(ad_result, resolution, n_regions, kwargs)

            add_probabilities_to_centroid(ad_result, res_key1)

        if "KMEANS" in clustering_method:
            print(f"Applying K-means Clustering...", flush=True)
            kmeans = KMeans(n_clusters=n_regions, random_state=random_state).fit(ad_result.obsm["X_pca"])
            ad_result.obs[res_key3] = pd.Categorical(kmeans.labels_.astype(str))
            add_probabilities_to_centroid(ad_result, res_key3)

    return ad_result


def _parallel_message_pass(
    ad: AnnData,
    radius: float,
    coord_type: str,
    set_diag: bool,
    mode: str,
    sample_rate: float,
):
    sq.gr.spatial_neighbors(ad, radius=radius, coord_type=coord_type, set_diag=set_diag)
    ad = custom_message_passing(ad, mode=mode, sample_rate=sample_rate)
    return ad


def custom_message_passing(adata: AnnData, mode: str = "l1_norm", sample_rate=1) -> AnnData:
    from scipy.linalg import sqrtm
    from scipy.sparse import csr_matrix

    if sample_rate == 0:
        adata.obsp['spatial_connectivities'] = np.eye(adata.shape[0])
        #print(adata.obsp['spatial_connectivities'], flush=0)
    elif sample_rate < 1:
        A = adata.obsp['spatial_connectivities']
        n_cells = A.shape[0]
        A_subsampled = np.zeros((n_cells,n_cells))

        for i in range(n_cells):
            nbr_idxs = A[i].nonzero()[1]
            n_samples = int(len(nbr_idxs) * sample_rate)
            nbr_idxs_subsampled = np.random.choice(nbr_idxs, size=n_samples, replace=False)
            for j in nbr_idxs_subsampled:
                A_subsampled[i,j] = 1
        print(len(nbr_idxs), flush=True)
        print(n_samples, flush=True)
        adata.obsp['spatial_connectivities'] = A_subsampled
    
    adata.obsp['spatial_connectivities'] = csr_matrix(adata.obsp['spatial_connectivities'])

    if mode == "l1_norm":
        A = adata.obsp["spatial_connectivities"]
        A_mod = np.asarray(A + np.eye(A.shape[0]))
        from sklearn.preprocessing import normalize
        affinity = normalize(A_mod, axis=1, norm="l1")
    else:
        # Plain A_mod multiplication
        A = adata.obsp["spatial_connectivities"]
        affinity = A

    adata.obsm['X_smoothed'] = affinity @ adata.X
    return adata


def low_variance_filter(adata: AnnData) -> AnnData:
    return adata[:, adata.var["std"] > adata.var["std"].median()]


def add_probabilities_to_centroid(
    adata: AnnData, col: str, name_to_output: str = None
) -> AnnData:
    from utag.utils import z_score
    from scipy.special import softmax

    if name_to_output is None:
        name_to_output = col + "_probabilities"

    mean = z_score(adata.to_df()).groupby(adata.obs[col]).mean()
    probs = softmax(adata.to_df() @ mean.T, axis=1)
    adata.obsm[name_to_output] = probs
    return adata


def evaluate_performance(
    adata: AnnData,
    batch_key: str = "Slide",
    truth_key: str = "DOM_argmax",
    pred_key: str = "cluster",
    method: str = "rand",
) -> Array:
    assert method in ["rand", "homogeneity"]
    from sklearn.metrics import rand_score, homogeneity_score

    score_list = []
    for key in adata.obs[batch_key].unique():
        batch = adata[adata.obs[batch_key] == key]
        if method == "rand":
            score = rand_score(batch.obs[truth_key], batch.obs[pred_key])
        elif method == "homogeneity":
            score = homogeneity_score(batch.obs[truth_key], batch.obs[pred_key])
        score_list.append(score)
    return score_list
