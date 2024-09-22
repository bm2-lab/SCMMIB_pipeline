# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: python (maxfuse)
#     language: python
#     name: maxfuse
# ---

# %%
import os
import sys
import json
import match
import metrics
import numpy as np
from numpy import std

import pandas as pd
import scanpy as sc
import anndata as ad
import maxfuse as mf
from scipy.io import mmread
from scipy.sparse import csr_matrix

from anndata import read_h5ad
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay


# %%
def MaxFuse_module(input_path, 
                   output_path, 
                   config):
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    adt = sc.read_h5ad(os.path.join(input_path, config['adt_h5ad_filename']))
    rna.var_names_make_unique()
    adt.var_names_make_unique()

    # Load Metadata
    metadata = pd.read_csv(os.path.join(input_path, config['metadata']), header=0)
    metadata.index=metadata[config['barcode_key']].values

    shared_features = rna.var_names.intersection(adt.var_names)

    # subset shared features
    rna_shared = rna[:, shared_features].copy()
    adt_shared = adt[:, shared_features].copy()

    # Get target_sum for normlization of shared matrix
    rna_counts = rna_shared.X.toarray().sum(axis=1)
    adt_counts = np.array(adt_shared.X.sum(axis=1))
    target_sum = (np.median(rna_counts.copy()) + np.median(adt_counts.copy())) / 2
    target_sum

    # process rna_shared
    sc.pp.normalize_total(rna_shared, target_sum=target_sum)
    sc.pp.log1p(rna_shared)
    sc.pp.scale(rna_shared)
    rna_shared = rna_shared.X.copy()

    # process adt_shared
    sc.pp.normalize_total(adt_shared, target_sum=target_sum)
    sc.pp.log1p(adt_shared)
    sc.pp.scale(adt_shared)
    adt_shared = adt_shared.X.copy()

    # process rna
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna)
    # only retain highly variable genes
    rna = rna[:, rna.var.highly_variable].copy()
    sc.pp.scale(rna)
    rna_active = rna.X

    # process adt
    sc.pp.normalize_total(adt)
    sc.pp.log1p(adt)
    sc.pp.scale(adt)
    adt_active = adt.X

    # preparations
    # call constructor for Fusor object
    fusor = mf.model.Fusor(
        shared_arr1=rna_shared,
        shared_arr2=adt_shared,
        active_arr1=rna_active,
        active_arr2=adt_active,
        labels1=None,
        labels2=None
    )

    # split_into_batches: 
    # reduce computational complexity
    fusor.split_into_batches(
        max_outward_size=5000,
        matching_ratio=3,
        metacell_size=2,
        verbose=True
    )

    # construct appropriate nearest-neighbor graphs 
    # for each modality with all features available
    fusor.construct_graphs(
        n_neighbors1=15,
        n_neighbors2=15,
        svd_components1=min(30, len(shared_features)-1),
        svd_components2=min(30, len(shared_features)-1),
        resolution1=2,
        resolution2=2,
        # if two resolutions differ less than resolution_tol
        # then we do not distinguish between then
        resolution_tol=0.1,
        verbose=True
    )

    # finding initial pivots using shared arrays
    fusor.find_initial_pivots(
        wt1=0.7, wt2=0.7,
        svd_components1=min(25, len(shared_features)-1), svd_components2=min(20, len(shared_features)-1)
    )

    # finding refined pivots
    fusor.refine_pivots(
        wt1=0.7, wt2=0.7,
        svd_components1=min(30, len(shared_features)-1), svd_components2=min(30, len(shared_features)-1),
        cca_components=min(20, len(shared_features)-1),
        n_iters=3,
        randomized_svd=False, 
        svd_runs=1,
        verbose=True
    )

    # filters away unreliable pivots
    fusor.filter_bad_matches(target='pivot', filter_prop=0.3)

    # get joint embedding
    rna_cca, adt_cca = fusor.get_embedding(
        active_arr1=fusor.active_arr1,
        active_arr2=fusor.active_arr2
    )

    # save latent
    rna_latent = pd.DataFrame(data=rna_cca, index=rna.obs.index)
    rna_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-MaxFuse-RNA-latent_1.csv"))

    adt_latent = pd.DataFrame(data=adt_cca, index=adt.obs.index)    
    adt_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-MaxFuse-ADT-latent_2.csv"))

    multi_latent = pd.merge(rna_latent, adt_latent, left_index=True, right_index=True, how='outer')
    multi_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-MaxFuse-multi-latent.csv"))


    # mf.metrics.get_foscttm(
    #     dist=mf.utils.cdist_correlation(rna_cca[:,:dim_use], adt_cca[:,:dim_use]),
    #     true_matching='identity'
    # )

#     # UMAP
#     cca_adata = ad.AnnData(
#         np.concatenate((rna_cca[:,:15], adt_cca[:,:15]), axis=0), 
#         dtype=np.float32
#     )
#     cca_adata.obs['data_type'] = ['rna'] * rna_cca.shape[0] + ['adt'] * adt_cca.shape[0]
#     cell_type = rna.obs['cell_type'].to_numpy()
#     cca_adata.obs['cell_type'] = list(cell_type) * 2
#     sc.pp.neighbors(cca_adata, n_neighbors=15)
#     sc.tl.umap(cca_adata)
#     # sc.pl.umap(cca_adata, color='data_type')

#     rna_adata = cca_adata[cca_adata.obs['data_type'] == 'rna']
#     rna_umap = pd.DataFrame(data=rna_adata.obsm['X_umap'], columns=["UMAP1","UMAP2"] , index=rna.obs.index)
#     rna_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-MaxFuse-RNA-umap.csv"))

#     adt_adata = cca_adata[cca_adata.obs['data_type'] == 'adt']
#     adt_umap = pd.DataFrame(data=adt_adata.obsm['X_umap'], columns=["UMAP1","UMAP2"] , index=adt.obs.index)
#     adt_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-MaxFuse-ADT-umap.csv"))

# %%
MaxFuse_module(input_path = sys.argv[1],
               output_path = sys.argv[2],
               config = json.load(open(sys.argv[3]))
              )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node"
# output_path = "/home/wsg/BM/results/task/scRNA+ADT/accuracy/SPATIAL/lymph_node/rep_1/run_MaxFuse"
# config = json.load(open("/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node/lymph_node.json"))

# %%
# MaxFuse_module(input_path,
#                output_path,
#                config) 

# %%

# %%

# %%
