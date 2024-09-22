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
import anndata as ad

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
                  gpu_index):
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

    data_list = [rna, adt]
    multi = ad.concat(data_list, axis = 1, join = 'outer')
    del multi.raw
    # multi.var_names.value_counts()
    data_path = os.path.join(output_path, "CITE-seq_merged_matrix.h5ad")
    multi.write_h5ad(data_path)

    adata = SCALEX(data_path, batch_categories=['RNA'], chunk_size=20001, min_features=0, min_cells=0, gpu=gpu_index)

    os.remove(data_path)

    # Save Results
    ## latent
    latent = pd.DataFrame(data=adata.obsm['latent'],
                              index=adata.obs_names)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-SCALEX-multi-latent.csv"))

    ## UMAP
    umap = pd.DataFrame(data=adata.obsm['X_umap'],
                            columns=["UMAP1", "UMAP2"],
                            index=adata.obs_names)
    umap.insert(2, "cluster", adata.obs['leiden'].values)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-SCALEX-multi-umap.csv"))

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
# input_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1"
# output_path = "/home/wsg/BM/results/task/scRNA+ADT/accuracy/SPATIAL/lymph_node_A1/rep_1/run_SCALEX"
# config = json.load(open("/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/lymph_node.json"))

# %%
# SCALEX_module(input_path, 
#               output_path, 
#               config,
#               gpu_index = gpu_index)

# %%
