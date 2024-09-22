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

import numpy as np
import scanpy as sc
import pandas as pd
from time import time
from math import ceil
from scanpy import read
from copy import deepcopy
from matplotlib import pyplot
from anndata import AnnData, read_h5ad
from scipy.stats import spearmanr, gamma, poisson

import torch
from torch import tensor
from torch.utils.data import DataLoader, TensorDataset

from sciPENN.sciPENN_API import sciPENN_API

# %%
from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def sciPENN_batch_module(input_path, 
                         output_path, 
                         config):
    torch.set_num_threads(5)
    
    # Make Dir
    if not os.path.exists(output_path):
        os.system("mkdir -p {}".format(output_path))
    
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    adt = sc.read_h5ad(os.path.join(input_path, config['adt_h5ad_filename']))
    rna.X = rna.X.astype(np.float32)
    adt.X = adt.X.astype(np.float32)
    rna.var_names_make_unique()
    adt.var_names_make_unique()

    # Load Metadata
    metadata = pd.read_csv(os.path.join(input_path, config['metadata']), header=0)
    metadata.index=metadata[config['barcode_key']].values

    # RNA preprocess
    rna.layers['counts'] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    #sc.pp.highly_variable_genes(rna, n_top_genes=1000, batch_key='batch')
    sc.pp.highly_variable_genes(rna, n_top_genes=3000)
    rna_hvg = rna[:, rna.var.highly_variable].copy()

    rna_tab = sc.AnnData(rna_hvg.X)
    obs = pd.DataFrame(rna_hvg.obs.index)
    obs['batch'] = rna_hvg.obs.batch.values
    obs['donor'] = rna_hvg.obs.batch.values
    rna_tab.obs = obs
    rna_tab.var = pd.DataFrame(rna_hvg.var.index)
    rna_tab

    # print(adt.shape)
    # adt = adt[:, adt.var.features.isin(['CD86','CD274','CD270','CD155'])]
    # print(adt.shape)

    adt_tab = sc.AnnData(adt.X)
    obs = pd.DataFrame(adt.obs.index)
    obs['batch'] = adt.obs.batch.values
    obs['donor'] = adt.obs.batch.values
    adt_tab.obs = obs
    adt_tab.var = pd.DataFrame(adt.var.index)
    adt_tab

    rna_tab.obs['batch'] = 0
    adt_tab.obs['batch'] = 0

    sciPENN = sciPENN_API(gene_trainsets = [rna_tab], # 多个数据集为[gene_tab1, gene_tab2]
                          protein_trainsets = [adt_tab], min_genes=0, train_batchkeys= ['donor'])

    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, 
                  decay_step = 0.1, lr = 10**(-3), weights_dir = "output", load = True)

    embedding = sciPENN.embed()

    embedding.obsm["latent"]=embedding.X
    sc.pp.neighbors(embedding, use_rep="latent", n_neighbors=30)
    sc.tl.umap(embedding, min_dist=0.3)
    sc.tl.louvain(embedding)

    latent = pd.DataFrame(embedding.X,index=embedding.obs.iloc[:,0])
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-sciPENN_batch-multi-latent.csv"))
    
    umap = pd.DataFrame(data=embedding.obsm["X_umap"], columns=["UMAP1","UMAP2"] , index=embedding.obs.iloc[:,0] )
    umap.insert(2,"cluster",embedding.obs['louvain'].values)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-sciPENN_batch-multi-umap.csv"))
    
    os.system("rm -r output")


# %%
sciPENN_batch_module(input_path=sys.argv[1],
                     output_path=sys.argv[2],
                     config=json.load(open(sys.argv[3]))
                    )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/BMMC/RNA+ADT/p10"
# output_path = "/home/wsg/BM/results/task/scRNA+ADT/accuracy/BMMC/p10/rep_1/run_sciPENN_batch"
# config=json.load(open("/home/wsg/BM/data/BMMC/RNA+ADT/p10/p10.json"))

# %%
# sciPENN_batch_module(input_path,
#                      output_path,
#                      config
#                     )

# %%
