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
#     display_name: python (SCALEX)
#     language: python
#     name: scalex
# ---

# %%
import os
import sys
import glob
import json

import torch
import scalex
from scalex import SCALEX
from scalex.plot import embedding
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

# https://github.com/jsxlei/SCALEX

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def SCALEX_module(input_path,
                  output_path,
                  config,
                  gpu_index
                 ):
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    data_list = [os.path.join(input_path, config['rna_h5ad_filename']), 
                 os.path.join(input_path, config['gam_h5ad_filename'])]

    # run_SCALEX = "/home/wsg/software/miniconda3/envs/benchmark/bin/SCALEX.py --data_list " + scRNA_Counts + ' ' + Gene_Activity_Matrix + ' --batch_categories RNA ATAC --min_features 0 --min_cells 0 -o ' + output_path + ' -g '+ str(gpu_index)
    # run_SCALEX = "/home/wsg/software/miniconda3/envs/benchmark/bin/SCALEX.py --data_list " + scRNA_Counts + ' ' + Gene_Activity_Matrix + ' --chunk_size 20001 --batch_categories RNA ATAC --min_features 0 --min_cells 0 -o ' + output_path + ' -g '+ str(gpu_index)

    adata = SCALEX(data_list, batch_categories=['RNA', 'ATAC'], chunk_size=20001, min_features=0, min_cells=0, gpu=gpu_index, ignore_umap=True)

    # sc.pl.umap(adata,color=['batch', 'leiden'],legend_fontsize=10, ncols=2)

    # Split adata
    rna = adata[adata.obs['batch'] == 'RNA', ]
    atac = adata[adata.obs['batch'] == 'ATAC', ]

    # Save Results
    ## latent
    rna_latent = pd.DataFrame(data=rna.obsm['latent'],
                              index=rna.obs_names)
    rna_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-SCALEX-RNA-latent.csv"))

    atac_latent = pd.DataFrame(data=atac.obsm['latent'],
                              index=atac.obs_names)
    atac_latent.index = [barcode.replace('.', '-') for barcode in atac_latent.index]
    atac_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-SCALEX-ATAC-latent.csv"))

#     ## UMAP
#     rna_umap = pd.DataFrame(data=rna.obsm['X_umap'],
#                             columns=["UMAP1", "UMAP2"],
#                             index=rna.obs_names)
#     rna_umap.insert(2, "cluster", rna.obs['leiden'].values)
#     rna_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-SCALEX-RNA-umap.csv"))

#     atac_umap = pd.DataFrame(data=atac.obsm['X_umap'],
#                             columns=["UMAP1", "UMAP2"],
#                             index=atac.obs_names)
#     atac_umap.insert(2, "cluster", rna.obs['leiden'].values)
#     atac_umap.index = [barcode.replace('.', '-') for barcode in atac_umap.index]
#     atac_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-SCALEX-ATAC-umap.csv"))
    
    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-SCALEX-gpu_memory.csv'), header=["gpu_memory"])


# %%
SCALEX_module(input_path = sys.argv[1], 
              output_path = sys.argv[2],
              config = json.load(open(sys.argv[3])),
              gpu_index = gpu_index
             )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
# output_path = "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/HSPC/p10/rep_1/run_SCALEX"
# config = json.load(open("/home/wsg/BM/data/HSPC/RNA+ATAC/p10/p10.json"))

# %%
# SCALEX_module(input_path,
#               output_path,
#               config,
#               gpu_index)
