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
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

from scipy import io

# %% [markdown]
# # Data Manipulation: convert h5ad to mtx

# %% [markdown]
# ## multiome

# %%
HSPC_METADATA = '/home/wsg/BM/data/HSPC/RawData/Multimodal_Single-Cell_Integration/metadata.csv'
metadata = pd.read_csv(HSPC_METADATA)

metadata_multi = metadata[metadata.technology=="multiome"]
metadata_cite = metadata[metadata.technology=="citeseq"]

# %%
metadata_multi.shape, metadata_cite.shape

# %%
# we only use train_data because test_data only have ATAC modality

# HSPC train_data ATAC-seq peak counts transformed
HSPC_multi_ATAC_path ='/home/wsg/BM/data/HSPC/RawData/Raw_Counts/train_multi_inputs_raw.h5'

# HSPC train_data RNA gene expression levels as library-size normalized and log1p transformed counts for the same cells
HSPC_multi_RNA_path ='/home/wsg/BM/data/HSPC/RawData/Raw_Counts/train_multi_targets_raw.h5'

# %%
# HSPC_multi_RNA_counts = pd.read_hdf(HSPC_multi_RNA_path)
HSPC_multi_RNA_counts.shape

# %%
# HSPC_multi_ATAC_counts = pd.read_hdf(HSPC_multi_ATAC_path)
HSPC_multi_ATAC_counts.shape

# %%
sum(HSPC_multi_RNA_counts.index == HSPC_multi_ATAC_counts.index)

# %%
metadata_multi_train = metadata_multi[metadata_multi['cell_id'].isin(HSPC_multi_RNA_counts.index)]
metadata_multi_train.shape

metadata_multi_train['barcode'] = metadata_multi_train['cell_id']
metadata_multi_train.set_index('cell_id', inplace=True)

# %%
print(metadata_multi_train.shape)
print(metadata_multi_train['cell_type'].value_counts())

# %%
np.random.seed(1234)
metadata_multi_train_p10 = metadata_multi_train.sample(frac=0.1)
metadata_multi_train_p10['cell_type'].value_counts()

# %%
p10_condition = HSPC_multi_RNA_counts.index.isin(metadata_multi_train_p10['barcode'])
HSPC_multi_RNA_counts_p10 = HSPC_multi_RNA_counts[p10_condition]

# %%
adata_rna = sc.AnnData(X=HSPC_multi_RNA_counts_p10.values, 
                           obs=metadata_multi_train_p10,
                           var=pd.DataFrame(index=HSPC_multi_RNA_counts_p10.columns))
adata_rna

# %%
output_path = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
# save hd5
adata_rna.write_h5ad("{}/HSPC-multiome-p10-RNA-counts.h5ad".format(output_path))

# %%

# %%
HSPC_multi_ATAC_counts_p10 = HSPC_multi_ATAC_counts[p10_condition]

# %%
adata_atac = sc.AnnData(X=HSPC_multi_ATAC_counts_p10.values, 
                           obs=metadata_multi_train_p10,
                           var=pd.DataFrame(index=HSPC_multi_ATAC_counts_p10.columns))
adata_atac

# %%
output_path = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
# save hd5
adata_atac.write_h5ad("{}/HSPC-multiome-p10-ATAC-peaks.h5ad".format(output_path))

# %%

# %%

# %%

# %% [markdown]
# ## CITE-seq

# %%
HSPC_METADATA = '/home/wsg/BM/data/HSPC/RawData/Multimodal_Single-Cell_Integration/metadata.csv'
metadata = pd.read_csv(HSPC_METADATA)

metadata_multi = metadata[metadata.technology=="multiome"]
metadata_cite = metadata[metadata.technology=="citeseq"]

# %%
metadata_multi.shape, metadata_cite.shape

# %%
# we only use train_data because test_data only have RNA modality

# HSPC train_data RNA gene expression levels as library-size normalized and log1p transformed counts (gene expression levels)
HSPC_cite_RNA_path ='/home/wsg/BM/data/HSPC/RawData/Raw_Counts/train_cite_inputs_raw.h5'

# HSPC train_data surface protein levels for the same cells that have been dsb normalized
HSPC_cite_ADT_path ='/home/wsg/BM/data/HSPC/RawData/Raw_Counts/train_cite_targets_raw.h5'

# %%
HSPC_cite_RNA_counts = pd.read_hdf(HSPC_cite_RNA_path)
HSPC_cite_RNA_counts.shape

# %%
HSPC_cite_ADT_counts = pd.read_hdf(HSPC_cite_ADT_path)
HSPC_cite_ADT_counts.shape

# %%
sum(HSPC_cite_RNA_counts.index == HSPC_cite_ADT_counts.index)

# %%
metadata_cite_train = metadata_cite[metadata_cite['cell_id'].isin(HSPC_cite_RNA_counts.index)]
metadata_cite_train.shape

metadata_cite_train['barcode'] = metadata_cite_train['cell_id']
metadata_cite_train.set_index('cell_id', inplace=True)

# %%
print(metadata_cite_train.shape)
print(metadata_cite_train['cell_type'].value_counts())

# %%
np.random.seed(1234)
metadata_cite_train_p10 = metadata_cite_train.sample(frac=0.1)
metadata_cite_train_p10['cell_type'].value_counts()

# %%
p10_condition = HSPC_cite_RNA_counts.index.isin(metadata_cite_train_p10['barcode'])
HSPC_cite_RNA_counts_p10 = HSPC_cite_RNA_counts[p10_condition]

# %%
adata_rna = sc.AnnData(X=HSPC_cite_RNA_counts_p10.values, 
                       obs=metadata_cite_train_p10,
                       var=pd.DataFrame(index=HSPC_cite_RNA_counts_p10.columns))
adata_rna

# %%
output_path = "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
# save hd5
adata_rna.write_h5ad("{}/HSPC-cite-p10-RNA-counts.h5ad".format(output_path))

# %%
adata_rna.obs.to_csv(output_path + '/metadata.csv')

# %%

# %%
HSPC_cite_ADT_counts_p10 = HSPC_cite_ADT_counts[p10_condition]

# %%
adata_adt = sc.AnnData(X=HSPC_cite_ADT_counts_p10.values, 
                        obs=metadata_cite_train_p10,
                        var=pd.DataFrame(index=HSPC_cite_ADT_counts_p10.columns))
adata_adt

# %%
output_path = "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
# save hd5
adata_adt.write_h5ad("{}/HSPC-multiome-p10-ADT-counts.h5ad".format(output_path))
