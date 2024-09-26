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

# %%
COVID_H5AD_PATH = "/home/wsg/BM/data/COVID19/E-MTAB-10026/covid_portal_210320_with_raw.h5ad"
covid = sc.read_h5ad(COVID_H5AD_PATH)

# %%
print(type(covid.X))
covid

# %%
covid_rna = covid[:, covid.var["feature_types"] == "Gene Expression"]
covid_adt = covid[:, covid.var["feature_types"] == "Antibody Capture"]

# %%
print(covid_rna)
print(covid_adt)

# %%
print(covid_rna.X[:10,:10].todense())
print(type(covid_rna.X))

# %%
print(covid_rna.layers['raw'][:10,:10].todense())
print(type(covid_rna.layers['raw'].copy()))

# %%
covid_rna.layers['processed'] = covid_rna.X.copy()
covid_rna.X = covid_rna.layers['raw'].copy()
covid_rna.X.todense()

# %%
del covid_rna.layers['processed']
del covid_rna.layers["raw"]

# %%
covid_rna.X

# %%
# covid_rna.X = covid_rna.X.astype(np.int32)
print(covid_rna.X[:10,:10].todense())
print(type(covid_rna.X.copy()))

# %%
covid_rna.var_names

# %%
output_path = "/home/wsg/BM/data/COVID19/RNA+ADT/RawData"
covid_rna.write_h5ad("{}/COVID19-CITE_seq-raw-RNA-counts.h5ad".format(output_path))

# %%
output_path = "/home/wsg/BM/data/COVID19/RNA+ADT/RawData"
covid_rna = sc.read_h5ad("{}/COVID19-CITE_seq-raw-RNA-counts.h5ad".format(output_path))

# %%
covid_rna.obs['barcode'] = covid_rna.obs_names

# %%
covid_rna.obs.to_csv("{}/metadata.csv".format(output_path))

# %%

# %%

# %%

# %%
adt_names = [name.split('AB_')[1] for name in covid_adt.var_names if len(name.split('AB_')) > 1]
covid_adt.var_names = adt_names
covid_adt.var_names

# %%
covid_adt.layers['processed'] = covid_adt.X.copy()
covid_adt.X = covid_adt.layers['raw'].copy()
covid_adt.X.todense()

# %%
del covid_adt.layers['processed']
del covid_adt.layers["raw"]

# %%
# covid_adt.X = covid_adt.X.astype(np.int32)
print(covid_adt.X[:10,:10].todense())
print(type(covid_adt.X.copy()))

# %%
output_path = "/home/wsg/BM/data/COVID19/RNA+ADT/RawData"
# save hd5
covid_adt.write_h5ad("{}/COVID19-CITE_seq-raw-ADT-counts.h5ad".format(output_path))

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
# Save RNA
## makr dir
# !mkdir /home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-RNA-counts.mtx
## save X to mtx
io.mmwrite('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-RNA-counts.mtx/matrix', covid_rna.X.T)
## save barcodes
with open('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-RNA-counts.mtx/barcodes.tsv', 'w') as f:
    for item in covid_rna.obs_names:
        f.write(item + '\n')      
## save features
with open('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-RNA-counts.mtx/features.tsv', 'w') as f:
    for item in covid_rna.var_names:
        f.write(item + '\n')
## gzip file
# !gzip /home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-RNA-counts.mtx/*
## save metadata
covid_rna.obs.to_csv('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/metadata.csv')

# %%
# Save ADT
## makr dir
# !mkdir /home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-ADT-counts.mtx
## save X to mtx
io.mmwrite('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-ADT-counts.mtx/matrix', covid_adt.X.T)
## save barcodes
with open('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-ADT-counts.mtx/barcodes.tsv', 'w') as f:
    for item in covid_adt.obs_names:
        f.write(item + '\n')      
## save features
with open('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-ADT-counts.mtx/features.tsv', 'w') as f:
    for item in covid_adt.var_names:
        f.write(item + '\n')
## gzip file
# !gzip /home/wsg/BM/data/COVID19/RNA+ADT/RawData/COVID19-CITE_seq-raw-ADT-counts.mtx/*
## save metadata
covid_adt.obs.to_csv('/home/wsg/BM/data/COVID19/RNA+ADT/RawData/metadata.csv')

# %%

# %%
