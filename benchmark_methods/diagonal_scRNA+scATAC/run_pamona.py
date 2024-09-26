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
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import torch

from pamona import Pamona

# %%
from nvitop import Device
devices = Device.all() 
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def pamona_module(input_path,
                  output_path,
                  config):
    
    torch.set_num_threads(5)

    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    gam = sc.read_h5ad(os.path.join(input_path, config['gam_h5ad_filename']))

    # Load Metadata
    metadata = pd.read_csv(os.path.join(input_path, config['metadata']), header=0)
    metadata.index=metadata[config['barcode_key']].values

    # RNA preprocess
    rna.layers['counts'] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    # sc.pp.highly_variable_genes(rna, n_top_genes=3000, batch_key='batch')
    sc.pp.highly_variable_genes(rna, n_top_genes=3000)
    rna_hvg = rna[:, rna.var.highly_variable].copy()

    ## ATAC
    sc.pp.filter_genes(atac, min_cells=int(rna.shape[0] * 0.001))
    sc.pp.filter_genes(atac, max_cells=int(rna.shape[0] * 0.1))

    Pa = Pamona.Pamona(Lambda=10, n_neighbors=40) 
    integrated_data, T = Pa.run_Pamona([rna_hvg.X.A, atac.X.A])    

    # Save Results
    df = pd.DataFrame(integrated_data[0])
    df.index = rna.obs_names
    df.to_csv(os.path.join(output_path, config["output_prefix"] + "-pamona-RNA-latent.csv"))


    df = pd.DataFrame(integrated_data[1])
    df.index = atac.obs_names
    df.to_csv(os.path.join(output_path, config["output_prefix"] + "-pamona-ATAC-latent.csv"))
    
    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-pamona-gpu_memory.csv'), header=["gpu_memory"])

# %%
pamona_module(input_path = sys.argv[1], 
              output_path = sys.argv[2],
              config = json.load(open(sys.argv[3]))
             )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/p10"
# output_path = "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/BMMC/p10_GAM/rep_1/run_pamona"
# config = json.load(open("/home/wsg/BM/data/BMMC/RNA+ATAC/p10/p10.json"))

# %%
# pamona_module(input_path,
#               output_path,
#               config)

# %%
