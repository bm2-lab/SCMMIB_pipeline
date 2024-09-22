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
                      config):
    
    torch.set_num_threads(5)

    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    # RNA preprocess
    rna.layers['counts'] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000)
    rna_hvg = rna[:, rna.var.highly_variable].copy()

    # ATAC preprocess
    atac.layers['counts'] = atac.X.copy()
    muon.atac.pp.tfidf(atac, scale_factor=1e4)
    atac.layers['tf-idf'] = atac.X.copy()
    atac.X = atac.layers['counts'].copy()
    sc.pp.normalize_total(atac, target_sum=1e4)
    sc.pp.log1p(atac)
    sc.pp.highly_variable_genes(atac, n_top_genes=30000)
    atac.layers['log-norm'] = atac.X.copy()
    atac_hvf = atac[:, atac.var.highly_variable].copy()

    # Merge Data
    adata = mtg.data.organize_multiome_anndatas(
        adatas = [[rna_hvg], [atac_hvf]],           # a list of anndata objects per modality, RNA-seq always goes first
        layers = [['counts'], ['log-norm']],        # if need to use data from .layers, if None use .X
    )

    # Setup anndata
    mtg.model.MultiVAE.setup_anndata(
        adata,
        rna_indices_end=3000,          # how many features in the rna-seq modality
        categorical_covariate_keys=None
    )

    # Initialize Model
    model = mtg.model.MultiVAE(
        adata,
        losses=['nb', 'mse'],
    )

    # Train model
    model.train() # default lr = 0.0005
    # model.train(lr = 0.0002) # set when n_cells >= 100k

    # model.plot_losses()
    model.get_latent_representation()

    ## latent
    latent = pd.DataFrame(adata.obsm['latent'], index=adata.obs_names)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-multigrate-multi-latent.csv"))

    # Visualize Result
    sc.pp.neighbors(adata, use_rep='latent')
    sc.tl.umap(adata)
    # sc.pl.umap(adata, color=['cell_type', 'batch'], frameon=False, ncols=1)

    ## UMAP
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"],
                        index=adata.obs_names)
    umap.insert(2, "cluster", adata.obs['cell_type'].to_list())
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-multigrate-multi-umap.csv"))
    
    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-multigrate-gpu_memory.csv'), header=["gpu_memory"])


# %%
multigrate_module(input_path = sys.argv[1], 
                  output_path = sys.argv[2],
                  config = json.load(open(sys.argv[3]))
                 )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/HSPC/p10/rep_1/run_multigrate"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
# multigrate_module(input_path,
#                   output_path,
#                   config)
