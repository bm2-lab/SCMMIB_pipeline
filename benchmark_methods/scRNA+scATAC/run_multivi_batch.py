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
#     display_name: python (scvi-env)
#     language: python
#     name: scvi-env
# ---

# %%
import os
import sys 
import json
import time
import pickle

import scvi
import torch
import random
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

# %%
from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def multivi_module(input_path, output_path, dataset, max_cpu=30):
    
    # Make Dir
    if not os.path.exists(output_path):
        os.system("mkdir -p {}".format(output_path))
    torch.set_num_threads(max_cpu)
    
    # Load Data
    rna = sc.read_h5ad("{}/{}-{}-{}-RNA-counts.h5ad".format(input_path, 
                                                            dataset["data_name"], 
                                                            dataset["data_type"], 
                                                            dataset["task_name"]))

    atac = sc.read_h5ad("{}/{}-{}-{}-ATAC-peaks.h5ad".format(input_path, 
                                                            dataset["data_name"], 
                                                            dataset["data_type"], 
                                                            dataset["task_name"]))


    # Filter Data
    ## RNA
    sc.pp.filter_genes(rna, min_cells=int(rna.shape[0] * 0.001))
    # sc.pp.filter_cells(rna, min_genes=3)
    
    ## ATAC
    sc.pp.filter_genes(atac, min_cells=int(rna.shape[0] * 0.001))
    # sc.pp.filter_cells(atac, min_genes=3)

    # Merge Data
    multi = anndata.concat([rna, atac], axis=1, join="outer")
    multi.var_names_make_unique()
    multi.obs = rna.obs

    # Add Metadata
    multi.var['modality'] = np.repeat(["Gene Expression", "Peaks"], [rna.shape[1], atac.shape[1]], axis=0)
    multi.obs['modality'] = np.repeat(["BMMC"], rna.shape[0], axis=0)

#     multi = sc.read_h5ad("{}/{}-{}-{}-multi-filtered.h5ad".format(input_path, dataset["data_name"],
#                                                     "multivi", dataset["task_type"]))

    # Train Data
#     batch_key='batch', 
    scvi.model.MULTIVI.setup_anndata(multi, batch_key='modality', categorical_covariate_keys=['batch'])
    mvi = scvi.model.MULTIVI(
        multi,
        n_genes=(multi.var['modality']=='Gene Expression').sum(),
        n_regions=(multi.var['modality']=='Peaks').sum(),
    )
    # mvi.view_anndata_setup()
    mvi.train(30, save_best=True, use_gpu=True)

    # Save Results
    ## prepare latent
    latent = mvi.get_latent_representation() # 20 dimension
    latent = pd.DataFrame(latent, index=multi.obs_names)
    ## save latent
    latent.to_csv("{}/{}-{}-{}-{}-{}-multi-latent.csv".format(output_path,
                                                        dataset["data_name"], 
                                                        dataset["data_type"], 
                                                        dataset["task_name"],
                                                        dataset["task_type"],
                                                        "multivi_batch"))

    # Visulize Data
    multi.obsm["X_multi_vi"] = latent
    sc.pp.neighbors(multi, use_rep="X_multi_vi")
    sc.tl.umap(multi, min_dist=0.2)
    sc.tl.louvain(multi)

    ## prepare UMAP
    umap = pd.DataFrame(data=multi.obsm["X_umap"],
                        columns=["UMAP1", "UMAP2"],
                        index=multi.obs_names)
    umap.insert(2, "cluster", multi.obs['louvain'].values)
    ## save UMAP
    umap.to_csv("{}/{}-{}-{}-{}-{}-multi-umap.csv".format(output_path,
                                                    dataset["data_name"], 
                                                    dataset["data_type"], 
                                                    dataset["task_name"],
                                                    dataset["task_type"],
                                                    "multivi_batch"))

    ## prepare imputation
    ## RNA
    # rna_imputation = vae.get_normalized_expression()
    ## ATAC
    # atac_imputation = vae.get_accessibility_estimates()

    ## save imputation
    ## RNA
    # rna_imputation.to_csv("{}/{}-{}-{}-RNA-imputation.csv".format(output_path, 
    #                                                               dataset["data_name"],
    #                                                               "multivi", 
    #                                                               dataset["task_type"]))
    ## ATAC
    # atac_imputation.to_csv("{}/{}-{}-{}-ATAC-imputation.csv".format(output_path, 
    #                                                                 dataset["data_name"], 
    #                                                                 "multivi", 
    #                                                                 dataset["task_type"]))

# %%
multivi_module(input_path=sys.argv[1], 
               output_path=sys.argv[2], 
               dataset=json.load(open(sys.argv[3])),
               max_cpu=30)

# %%
# input_path = "/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/c1k"
# output_path = "/home/wsg/BM/pipeline/results/BMMC/RNA+ATAC/multiome/scalability/pair/c1k/run_multivi"
# dataset=json.load(open("/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/c1k.json"))


# %%
# multivi_module(input_path,
#                output_path,
#                dataset,
#                max_cpu=30)

# %%
