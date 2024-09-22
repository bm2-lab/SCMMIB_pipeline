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
import pickle

import torch
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io

from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt

# https://github.com/epurdom/cobolt/blob/master/docs/tutorial.ipynb

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def cobolt_module(input_path,
                  output_path,
                  config):
    torch.set_num_threads(5)
    
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    barcodes = rna.obs_names.values.astype('str')
    features = rna.var_names.values.astype('str')
    rna_data = SingleData("GeneExpr", "SCMMIB", features, rna.X.copy(), barcodes)

    barcodes = atac.obs_names.values.astype('str')
    features = atac.var_names.values.astype('str')
    atac_data = SingleData("ChromAccess", "SCMMIB", features, atac.X.copy(), barcodes)

#     # Quality filtering on features.
#     rna_data.filter_features(upper_quantile=0.99, lower_quantile=0.7)
#     atac_data.filter_features(upper_quantile=0.99, lower_quantile=0.7)

    # Merge Data
    multi_dt = MultiomicDataset.from_singledata(rna_data, atac_data)

    # Train Data
    # Try a smaller learning rate.
    model = Cobolt(dataset=multi_dt, lr=0.001, n_latent=100)
    model.train(num_epochs=20)

    ## Calculate the corrected latent variables
    model.calc_all_latent()

    ## Get the latent variables (10 dimensions)
    latent_embed = model.get_all_latent()

    # Cluster Data
    model.clustering(algo="louvain", resolution=0.5)
    latent_clusters = model.get_clusters("louvain", resolution=0.5)

    # Visualize Data
    model.scatter_plot(reduc="UMAP", algo="louvain", resolution=0.5, s=0.2)
    umap_reduc = model.reduction["UMAP2"]["embedding"]

    # Save Results
    ## save UMAP
    umap = pd.DataFrame(np.column_stack((umap_reduc, latent_clusters)), columns=["UMAP1", "UMAP2", 'cluster'], index = latent_embed[1])
    umap.index = [index.split('~')[1] for index in umap.index]
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-cobolt-multi-umap.csv"))

    ## save latent
    latent = pd.DataFrame(latent_embed[0], index=latent_embed[1])
    latent.index = [index.split('~')[1] for index in latent.index]
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-cobolt-multi-latent.csv"))
    
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
        gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-cobolt-gpu_memory.csv'), header=["gpu_memory"])


# %%
cobolt_module(input_path = sys.argv[1], 
              output_path = sys.argv[2],
              config = json.load(open(sys.argv[3]))
             )


# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
# output_path = "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/HSPC/p10/rep_1/run_test"
# config = json.load(open("/home/wsg/BM/data/HSPC/RNA+ATAC/p10/p10.json"))


# %%
# cobolt_module(input_path,
#               output_path,
#               config)

# %%
