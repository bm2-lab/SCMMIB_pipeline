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

# from scvi.model.MULTIVI import setup_anndata, synthetic_iid, transfer_anndata_setup
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.io as sp_io
from anndata import AnnData
import numpy as np
import scipy.sparse as sparse
import torch

import anndata
import matplotlib.pyplot as plt
import mudata as md
import muon
import scanpy as sc
import scvi

# %%
from nvitop import Device
devices = Device.all() 
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def totalVI_module(input_path, 
                   output_path, 
                   config):
    torch.set_num_threads(5)

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

    adata = rna
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["protein_expression"]=adt.X.todense()
    adata.uns['protein_names']=adt.var_names

    # Set up
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata,
                                n_top_genes=3000,
                                flavor="seurat_v3",
    #                                 batch_key="batch",
                                subset=True,
                                layer="counts")

    scvi.model.TOTALVI.setup_anndata(
        adata,
        protein_expression_obsm_key="protein_expression",
        layer="counts"
    #         batch_key="batch"
    )

    # Prepare and run model
    vae = scvi.model.TOTALVI(adata, latent_distribution="normal")
    vae.train(use_gpu=True)

    # Save Results
    ## model
    vae.save('{}/totalVI_model.pkl'.format(output_path), overwrite=True)

    ## prepare
    adata.obsm["X_totalVI"] = vae.get_latent_representation()
    rna, protein = vae.get_normalized_expression(
        n_samples=25,
        return_mean=True,
    )
    adata.layers["denoised_rna"], adata.obsm["denoised_protein"] = rna, protein

    ## latent
    latent = vae.get_latent_representation()
    prior_adata = anndata.AnnData(X=adt.X)
    prior_adata.obsm["X_totalVI"] = latent
    ## save latent
    latent = pd.DataFrame(latent, index=adata.obs_names)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-totalVI-multi-latent.csv"))

    ## umap
    sc.pp.neighbors(prior_adata, use_rep="X_totalVI", n_neighbors=30)
    sc.tl.umap(prior_adata, min_dist=0.3)
    # save umap
    umap = pd.DataFrame(data=prior_adata.obsm["X_umap"], columns=["UMAP1","UMAP2"] , index=adata.obs_names)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-totalVI-multi-umap.csv"))

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-totalVI-gpu_memory.csv'), header=["gpu_memory"])



# %%
totalVI_module(input_path = sys.argv[1],
               output_path = sys.argv[2],
               config = json.load(open(sys.argv[3]))
              )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1"
# output_path = "/home/wsg/BM/results/task/scRNA+ADT/accuracy/SPATIAL/lymph_node_A1/rep_1/run_totalVI"
# config = json.load(open("/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/lymph_node.json"))

# %%
# totalVI_module(input_path,
#                output_path,
#                config) 
