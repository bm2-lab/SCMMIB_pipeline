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
#     display_name: python (multivi)
#     language: python
#     name: multivi
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
from sklearn.metrics import precision_recall_curve, auc

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def multivi_module(input_path, 
                   output_path, 
                   config):
    torch.set_num_threads(5)
    
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))


    # Filter Data
    # filter features to remove those that appear in fewer than 1% of the cells
    ## RNA
    sc.pp.filter_genes(rna, min_cells=int(rna.shape[0] * 0.001))
    # sc.pp.filter_cells(rna, min_genes=3)

    ## ATAC
    sc.pp.filter_genes(atac, min_cells=int(rna.shape[0] * 0.001))
    # sc.pp.filter_cells(atac, min_genes=3)

    # Merge Data
    multi = anndata.concat([rna, atac], axis=1, join="outer")
    multi.var_names_make_unique()

    # Add Metadata
    multi.var['modality'] = np.repeat(["Gene Expression", "Peaks"], [rna.shape[1], atac.shape[1]], axis=0)
    multi.obs['modality'] = np.repeat(["paired"], multi.shape[0], axis=0)

    # Train Data
    scvi.model.MULTIVI.setup_anndata(multi, batch_key='modality')
    mvi = scvi.model.MULTIVI(
        multi,
        n_genes=(multi.var['modality']=='Gene Expression').sum(),
        n_regions=(multi.var['modality']=='Peaks').sum(),
    )
    # mvi.view_anndata_setup()
    mvi.train(30, save_best=True, use_gpu=gpu_index)

    # Save Results
    ## save latent
    latent = mvi.get_latent_representation() # 20 dimension
    latent = pd.DataFrame(latent, index=multi.obs_names)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-multivi-multi-latent.csv"))

    # Visulize Data
    multi.obsm["X_multi_vi"] = latent
    sc.pp.neighbors(multi, use_rep="X_multi_vi")
    sc.tl.umap(multi, min_dist=0.2)
    sc.tl.louvain(multi)

    ## save UMAP
    umap = pd.DataFrame(multi.obsm["X_umap"], columns=["UMAP1", "UMAP2"],index=multi.obs_names)
    umap.insert(2, "cluster", multi.obs['louvain'].values)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-multivi-multi-umap.csv"))

    from sklearn.metrics import precision_recall_curve, auc

    # save imputation
    ## RNA
    imputation_rna = mvi.get_normalized_expression() # imputation
    imputation_rna.to_csv(os.path.join(output_path, config["output_prefix"] + "-multivi-imputation-rna.csv.gz"), compression='gzip')

    ## ATAC
    imputation_atac = mvi.get_accessibility_estimates() # imputation
    imputation_atac.to_csv(os.path.join(output_path, config["output_prefix"] + "-multivi-imputation-atac.csv.gz"), compression='gzip')

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-multivi-gpu_memory.csv'), header=["gpu_memory"])


# %%
multivi_module(input_path = sys.argv[1], 
               output_path = sys.argv[2], 
               config = json.load(open(sys.argv[3]))
              )

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/data/test/run_multivi"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
# multivi_module(input_path, 
#                output_path, 
#                config)

# %%
