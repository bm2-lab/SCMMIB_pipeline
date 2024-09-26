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
#     display_name: python (benchmark)
#     language: python
#     name: benchmark
# ---

# %%
import os
import sys 
import json

import uniport as up
import numpy as np
import pandas as pd
import torch
import scanpy as sc
import os
import episcanpy as epi
from sklearn.preprocessing import MinMaxScaler
import scipy.sparse

# https://uniport.readthedocs.io/en/latest/examples/vertical_snare_cell_line.html

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def uniPort_module(input_path, 
                   output_path, 
                   config):
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    # Set up
    atac.obs['domain_id'] = 0
    atac.obs['domain_id'] = atac.obs['domain_id'].astype('category')
    atac.obs['source'] = 'ATAC'

    rna.obs['domain_id'] = 1
    rna.obs['domain_id'] = rna.obs['domain_id'].astype('category')
    rna.obs['source'] = 'RNA'

    # preprocess 
    ## RNA
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000, inplace=False, subset=True)
    # up.batch_scale(rna) # chunk_size不可以被cell数量整除 chunk_size=20000 改为 chunk_size=20001
    up.batch_scale(rna, chunk_size=20001)

    ## ATAC
    atac.X[atac.X>1] = 1
    epi.pp.select_var_feature(atac, nb_features=2000, show=False, copy=False)
    sc.pp.normalize_total(atac)
    # up.batch_scale(atac) # chunk_size不可以被cell数量整除 chunk_size=20000 改为 chunk_size=20001
    up.batch_scale(atac, chunk_size=20001)

    # Run uniport
    # comb_tab = up.Run(adatas=[atac, rna], mode='v', lr=0.0001, iteration=10000, num_workers=False) # cell数量50k时报错，调低学习率后解决 (lr=0.001 改为 lr=0.0001)
    comb_tab = up.Run(adatas=[atac, rna], mode='v', lr=0.0001, iteration=10000, num_workers=False)

    # cluster
    sc.pp.neighbors(comb_tab, use_rep='latent')
    sc.tl.umap(comb_tab, min_dist=0.1)
    sc.tl.louvain(comb_tab)

    # Save results
    ## latent
    latent = pd.DataFrame(comb_tab.obsm["latent"], index=comb_tab.obs_names)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-uniPort-multi-latent.csv"))

    ## UMAP
    umap = pd.DataFrame(comb_tab.obsm["X_umap"], columns=["UMAP1","UMAP2"], index=comb_tab.obs_names)
    umap.insert(2, "cluster", comb_tab.obs['louvain'].values)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-uniPort-multi-umap.csv"))

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    if len(gpu_memory):
        gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-uniPort-gpu_memory.csv'), header=["gpu_memory"])


# %%
uniPort_module(input_path = sys.argv[1], 
               output_path = sys.argv[2], 
               config = json.load(open(sys.argv[3]))
              )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/10x_PBMC/RNA+ATAC/RawData"
# output_path = "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/10x_PBMC/RawData/rep_1/run_uniPort"
# config = json.load(open("/home/wsg/BM/data/10x_PBMC/RNA+ATAC/RawData/RawData.json"))

# %%
# uniPort_module(input_path,
#                output_path,
#                config
#               )

# %%
