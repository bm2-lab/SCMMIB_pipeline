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
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def totalVI_module(input_path, output_path, dataset, max_cpu=30):
#     torch.set_num_threads(30)
    # Make Dir
    if not os.path.exists(output_path):
        os.system("mkdir -p {}".format(output_path))
    
    # Load Data
    rna = sc.read_h5ad("{}/{}-{}-{}-RNA-counts.h5ad".format(input_path, 
                                                            dataset["data_name"],
                                                            dataset["data_type"],
                                                            dataset["task_name"]))

    adt = sc.read_h5ad("{}/{}-{}-{}-ADT-counts.h5ad".format(input_path,
                                                            dataset["data_name"],
                                                            dataset["data_type"],
                                                            dataset["task_name"]))

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
                                # batch_key="batch",
                                subset=True,
                                layer="counts")

    scvi.model.TOTALVI.setup_anndata(
        adata,
        protein_expression_obsm_key="protein_expression",
        layer="counts",
        batch_key="batch"
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
    ## umap
    sc.pp.neighbors(prior_adata, use_rep="X_totalVI", n_neighbors=30)
    sc.tl.umap(prior_adata, min_dist=0.3)
    ## save latent
    latent = pd.DataFrame(latent, index=adata.obs_names)
    latent.to_csv("{}/{}-{}-{}-{}-{}-multi-latent.csv".format(output_path,
                                                           dataset["data_name"],
                                                           dataset["data_type"],
                                                           dataset["task_name"],
                                                           dataset["task_type"],
                                                           'totalVI'))
    # save umap
    umap = pd.DataFrame(data=prior_adata.obsm["X_umap"], columns=["UMAP1","UMAP2"] , index=adata.obs_names)
    umap.to_csv("{}/{}-{}-{}-{}-{}-multi-umap.csv".format(output_path,
                                                       dataset["data_name"],
                                                       dataset["data_type"],
                                                       dataset["task_name"],
                                                       dataset["task_type"],
                                                       'totalVI'))


# %%
totalVI_module(input_path=sys.argv[1],
               output_path=sys.argv[2],
               dataset=json.load(open(sys.argv[3])),
               max_cpu=30)

# %%
# input_path = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/p10/downsample/R10_A10"
# output_path = "/home/wsg/BM/pipeline/results/BMMC/RNA+ADT/CITE-seq/downsample/p10/R10_A10/run_totalVI_batch"
# dataset=json.load(open("/NFS_home/NFS_home_2/wsg/BM/pipeline/results/BMMC/RNA+ADT/CITE-seq/downsample/p10/cromwell-executions/main/9508ca2b-21e9-4636-9807-de419d0ab88e/call-run_totalVI_batch/execution/write_json_fb872a6d0d3de8886957fa3a4485d649.tmp"))

# %%
# totalVI_module(input_path,
#                output_path,
#                dataset,
#                max_cpu=30) 

# %%
