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
import episcanpy as epi
from sklearn.preprocessing import MinMaxScaler
import scipy.sparse

# https://uniport.readthedocs.io/en/latest/examples/PBMC/pbmc_integration.html


# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def uniPort_module(input_path,
                   output_path,
                   config,
                   gpu_index
                  ):
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    gam = sc.read_h5ad(os.path.join(input_path, config['gam_h5ad_filename']))

    # setup object
    gam.obs['domain_id'] = 0
    gam.obs['domain_id'] = gam.obs['domain_id'].astype('category')
    gam.obs['source'] = 'ATAC'
    rna.obs['domain_id'] = 1
    rna.obs['domain_id'] = rna.obs['domain_id'].astype('category')
    rna.obs['source'] = 'RNA'

    # filter
    up.filter_data(gam, min_features=0, min_cells=200)
    up.filter_data(rna, min_features=0, min_cells=200)

    # merge
    comp = gam.concatenate(rna, join='inner', batch_key='domain_id')

    # preprocess
    sc.pp.normalize_total(comp)
    sc.pp.log1p(comp)
    sc.pp.highly_variable_genes(comp, n_top_genes=2000, inplace=False, subset=True)
    # up.batch_scale(comp)
    up.batch_scale(comp, chunk_size=20001)

    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, inplace=False, subset=True)
    # up.batch_scale(rna)
    up.batch_scale(rna, chunk_size=20001)

    sc.pp.normalize_total(gam)
    sc.pp.log1p(gam)
    sc.pp.highly_variable_genes(gam, n_top_genes=2000, inplace=False, subset=True)
    # up.batch_scale(gam)
    up.batch_scale(gam, chunk_size=20001)

    # run uniPort
    uniprot_model = up.Run(adatas=[gam,rna], adata_cm=comp, lambda_s=1.0, gpu=gpu_index)

    # cluster
    sc.pp.neighbors(uniprot_model, use_rep='latent')
    sc.tl.umap(uniprot_model, min_dist=0.1)
    sc.tl.louvain(uniprot_model)

    # save results
    def replace_last(string, old, new):
        return new.join(string.rsplit(old, 1))
    ## latent
    df=pd.DataFrame(uniprot_model.obsm['latent'],index=uniprot_model.obs.index.values)

    atac_latent = df.iloc[int(df.shape[0]/2):,]
    atac_latent.index = [replace_last(barcode, '-1', '') for barcode in atac_latent.index]
    atac_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-uniPort-ATAC-latent.csv"))

    rna_latent = df.iloc[:int(df.shape[0]/2),]
    rna_latent.index = [replace_last(barcode, '-0', '') for barcode in rna_latent.index]
    rna_latent.index = [barcode.replace('.', '-') for barcode in rna_latent.index]
    rna_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-uniPort-RNA-latent.csv"))

    ## umap
    df=pd.DataFrame(uniprot_model.obsm['X_umap'],
                    index=uniprot_model.obs.index.values)
    df.insert(2,"cluster",uniprot_model.obs['louvain'].values)

    atac_umap = df.iloc[int(df.shape[0]/2):,]
    atac_umap.index = [replace_last(barcode, '-1', '') for barcode in atac_umap.index]
    atac_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-uniPort-ATAC-umap.csv"))

    rna_umap = df.iloc[:int(df.shape[0]/2),]
    rna_umap.index = [replace_last(barcode, '-0', '') for barcode in rna_umap.index]
    rna_umap.index = [barcode.replace('.', '-') for barcode in rna_umap.index]
    rna_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-uniPort-RNA-umap.csv"))
    
    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-uniPort-gpu_memory.csv'), header=["gpu_memory"])


# %%
uniPort_module(input_path = sys.argv[1], 
               output_path = sys.argv[2], 
               config = json.load(open(sys.argv[3])),
               gpu_index = gpu_index
              )

# %%

# %%
# input_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/p10"
# output_path = "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/BMMC/p10_GAM/rep_1/run_uniPort"
# config = json.load(open("/home/wsg/BM/data/BMMC/RNA+ATAC/p10/p10.json"))

# %%
# uniPort_module(input_path,
#                output_path,
#                config,
#                gpu_index)

# %%
