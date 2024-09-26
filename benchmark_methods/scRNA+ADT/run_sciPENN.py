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
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def sciPENN_module(input_path, 
                   output_path, 
                   config):
    torch.set_num_threads(5)

    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
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
    sc.pp.highly_variable_genes(rna, n_top_genes=3000)
    rna = rna[:, rna.var.highly_variable].copy()

    rna.obs['batch'] = 0
    adt.obs['batch'] = 0

    sciPENN = sciPENN_API(gene_trainsets = [rna],  # 多个数据集为[gene_tab1, gene_tab2]
                          protein_trainsets = [adt], min_genes=0)

    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, 
                  decay_step = 0.1, lr = 10**(-3), weights_dir = "output", load = True)

    embedding = sciPENN.embed()

    embedding.obsm["latent"]=embedding.X
    sc.pp.neighbors(embedding, use_rep="latent", n_neighbors=30)
    sc.tl.umap(embedding, min_dist=0.3)
    sc.tl.louvain(embedding)

    os.system("rm -r output")

    latent = pd.DataFrame(embedding.X,index=embedding.obs.iloc[:,0])
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-sciPENN-multi-latent.csv"))

    umap = pd.DataFrame(data=embedding.obsm["X_umap"], columns=["UMAP1","UMAP2"] , index=embedding.obs.iloc[:,0] )
    umap.insert(2,"cluster",embedding.obs['louvain'].values)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-sciPENN-multi-umap.csv"))

    # gpu_memory
    pid = os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-sciPENN-gpu_memory.csv'), header=["gpu_memory"])


# %%
sciPENN_module(input_path = sys.argv[1],
               output_path = sys.argv[2],
               config = json.load(open(sys.argv[3]))
              )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
# output_path = "/home/wsg/BM/results/task/scRNA+ADT/accuracy/HSPC/p10/rep_1/run_sciPENN"
# config = json.load(open("/home/wsg/BM/data/HSPC/RNA+ADT/p10/p10.json"))

# %%
# sciPENN_module(input_path,
#                output_path,
#                config) 

# %%
