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
#     display_name: python (SpatialGlue)
#     language: python
#     name: spatialglue
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

import SpatialGlue
from SpatialGlue import preprocess

# https://spatialglue-tutorials.readthedocs.io/en/latest/Tutorial%202_data%20integration%20for%20mouse%20thymus%20Stereo-CITE-seq.html#

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def SpatialGlue_module(input_path,
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

    # Specify data type
    data_type = 'CITE-seq'

    # Fix random seed
    from SpatialGlue.preprocess import fix_seed
    random_seed = 2022
    fix_seed(random_seed)

    # Pre-processing data
    from SpatialGlue.preprocess import clr_normalize_each_cell, pca
    # RNA
    sc.pp.filter_genes(rna, min_cells=10)
    # sc.pp.filter_cells(rna, min_genes=80)

    sc.pp.filter_genes(adt, min_cells=50)
    adt = adt[rna.obs_names].copy()

    sc.pp.highly_variable_genes(rna, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)

    rna_high =  rna[:, rna.var['highly_variable']]
    rna.obsm['feat'] = pca(rna_high, n_comps=adt.n_vars-1)

    # Protein
    adt = clr_normalize_each_cell(adt)
    adt.obsm['feat'] = pca(adt, n_comps=adt.n_vars-1)

    # Constructing neighbor graph
    from SpatialGlue.preprocess import construct_neighbor_graph
    data = construct_neighbor_graph(rna, adt, datatype=data_type)

    # Training the model
    # define model
    from SpatialGlue.SpatialGlue_pyG import Train_SpatialGlue
    model = Train_SpatialGlue(data, datatype=data_type, device=gpu_index)

    # train model
    output = model.train()

    adata = rna.copy()
    adata.obsm['emb_latent_omics1'] = output['emb_latent_omics1']
    adata.obsm['emb_latent_omics2'] = output['emb_latent_omics2']
    adata.obsm['SpatialGlue'] = output['SpatialGlue']
    adata.obsm['alpha'] = output['alpha']
    adata.obsm['alpha_omics1'] = output['alpha_omics1']
    adata.obsm['alpha_omics2'] = output['alpha_omics2']

    # Save Results
    ## save latent
    latent = pd.DataFrame(data=adata.obsm['SpatialGlue'],  index=adata.obs.index)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-multi-latent.csv"))


# %%
SpatialGlue_module(input_path = sys.argv[1], 
                   output_path = sys.argv[2],
                   config = json.load(open(sys.argv[3]))
                  )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1"
# output_path = "/home/wsg/BM/results/task/spatial_scRNA+ADT/accuracy/lymph_node_A1/rep_1/run_SpatialGlue"
# config = json.load(open("/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/lymph_node.json"))

# %%
# SpatialGlue_module(input_path,
#                    output_path,
#                    config)
