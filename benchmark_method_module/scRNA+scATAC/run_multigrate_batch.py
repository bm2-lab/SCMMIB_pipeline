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
#     display_name: python (multigrate)
#     language: python
#     name: multigrate
# ---

# %%
import warnings
warnings.filterwarnings('ignore')

import os
import sys 
import json
import pickle

import torch
import scanpy as sc
import numpy as np
import pandas as pd
import muon
import multigrate as mtg
# https://multigrate.readthedocs.io/en/latest/notebooks/paired_integration_multiome.html

# %%
from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def multigrate_module(input_path,
                      output_path,
                      dataset,
                      max_cpu=30):
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
    
    # RNA preprocess
    rna.layers['counts'] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor='seurat', batch_key='batch')
    rna
    rna_hvg = rna[:, rna.var.highly_variable].copy()
    rna_hvg
    
    # ATAC preprocess
    atac.layers['counts'] = atac.X.copy()
    muon.atac.pp.tfidf(atac, scale_factor=1e4)
    atac
    atac.layers['tf-idf'] = atac.X.copy()
    atac.X = atac.layers['counts'].copy()
    sc.pp.normalize_total(atac, target_sum=1e4)
    sc.pp.log1p(atac)
    sc.pp.highly_variable_genes(atac, n_top_genes=30000, flavor='seurat', batch_key='batch')
    atac.layers['log-norm'] = atac.X.copy()
    atac_hvf = atac[:, atac.var.highly_variable].copy()
    atac_hvf
    
    # Merge Data
    adata = mtg.data.organize_multiome_anndatas(
        adatas = [[rna_hvg], [atac_hvf]],           # a list of anndata objects per modality, RNA-seq always goes first
        layers = [['counts'], ['log-norm']],    # if need to use data from .layers, if None use .X
    )
#     adata
    
    # Setup anndata
    mtg.model.MultiVAE.setup_anndata(
        adata,
        rna_indices_end=4000, # how many features in the rna-seq modality
        categorical_covariate_keys=["batch"]
    )
    
    # Initialize Model
    model = mtg.model.MultiVAE(
        adata,
        losses=['nb', 'mse'],
    )
    
    # Train model
    model.train()
    
    model.plot_losses()
    model.get_latent_representation()
    
    ## latent
    latent = pd.DataFrame(data=adata.obsm['latent'],
                          index=adata.obs_names)
    type(latent)
    
    latent.to_csv("{}/{}-{}-{}-{}-{}-multi-latent.csv".format(output_path,
                                                       dataset["data_name"],
                                                       dataset["data_type"],
                                                       dataset["task_name"],
                                                       dataset["task_type"],
                                                       "multigrate_batch"))
    
    # Visualize Result
    sc.pp.neighbors(adata, use_rep='latent')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['cell_type', 'batch'], frameon=False, ncols=1)
    
    ## UMAP
    umap = pd.DataFrame(data=adata.obsm["X_umap"],
                      columns=["UMAP1", "UMAP2"],
                      index=adata.obs_names)
    umap.insert(2, "cluster", adata.obs['cell_type'].to_list())
    umap.to_csv("{}/{}-{}-{}-{}-{}-multi-umap.csv".format(output_path,
                                                    dataset["data_name"],
                                                    dataset["data_type"],
                                                    dataset["task_name"],
                                                    dataset["task_type"],
                                                    "multigrate_batch"))



# %%
multigrate_module(input_path=sys.argv[1], 
                  output_path=sys.argv[2],
                  dataset=json.load(open(sys.argv[3])),
                  max_cpu=30)

# %%
# input_path = "/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/c1k"
# output_path = "/home/wsg/BM/pipeline/results/BMMC/RNA+ATAC/multiome/scalability/pair/c1k/run_multigrate"
# dataset=json.load(open("/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/c1k.json"))

# %%
# multigrate_module(input_path,
#                   output_path,
#                   dataset,
#                   max_cpu=30)
