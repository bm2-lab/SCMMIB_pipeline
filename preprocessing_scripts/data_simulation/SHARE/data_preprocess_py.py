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


# %%
def save_pickle(data, file_name):
    f = open(file_name, "wb")
    pickle.dump(data, f)
    f.close()

def load_pickle(file_name):
    f = open(file_name, "rb+")
    data = pickle.load(f)
    f.close()
    return data



# %%
import collections
import scipy.sparse as sp_sparse
import tables
 
CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
 
def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
         
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()
            
        return CountMatrix(feature_ref, barcodes, matrix)


# %%
def concatenate_csc_matrices_by_columns(matrix1, matrix2):
    new_data = np.concatenate((matrix1.data, matrix2.data))
    new_indices = np.concatenate((matrix1.indices, matrix2.indices))
    new_ind_ptr = matrix2.indptr + len(matrix1.data)
    new_ind_ptr = new_ind_ptr[1:]
    new_ind_ptr = np.concatenate((matrix1.indptr, new_ind_ptr))

    return scipy.sparse.csc_matrix((new_data, new_indices, new_ind_ptr))


# %% [markdown]
# # Data Manipulation: convert h5ad to mtx

# %% [markdown]
# ## RNA + ATAC

# %%
RNA_PATH = "/home/wsg/BM/data/SHARE/GB/shareseq_mouse_skin_rna.h5ad"
ATAC_PATH = "/home/wsg/BM/data/SHARE/GB/shareseq_mouse_skin_atac.h5ad"

# %%
import scanpy as sc
from scipy import io

rna = sc.read_h5ad(RNA_PATH)
atac = sc.read_h5ad(ATAC_PATH)

# %%
print(rna)
print(atac)

# %%
rna.obs
# adata.obs['bc'] = 
rna.obs['rna.bc'].str.replace(',', '.')

# %%
atac.obs

# %%
rna.obs['rna.bc'].values.isin(atac.obs.index)

# %%
atac.obs.index.isin(rna.obs['rna.bc'].str.replace(',', '.'))

# %%

# %%

# %%

# %%
test = sc.read_h5ad("/home/wsg/BM/data/SHARE/GB/mouse_skin_shareseq_atac_10k.h5ad")

# %%
test.var

# %%

# %%

# %%

# %%

# %%

# %%

# %%
BMMC_rna.layers['processed'] = BMMC_rna.X.copy()
BMMC_rna.X = BMMC_rna.layers['counts'].copy()
BMMC_rna.X.todense()

# %%
BMMC_atac.layers['binarized'] = BMMC_atac.X.copy()
BMMC_atac.X = BMMC_atac.layers['counts'].copy()
BMMC_atac.X.todense()

# %%
# print(BMMC_atac.X.T)
print(BMMC_atac.layers['binarized'].T)

# %%
output_path = "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/RawData"

# save hd5
BMMC_rna.write_h5ad("{}/BMMC-multiome-raw-RNA-counts.h5ad".format(output_path))
# save hd5
BMMC_atac.write_h5ad("{}/BMMC-multiome-raw-ATAC-peaks.h5ad".format(output_path))

# %%
# Save RNA
## makr dir
# !mkdir RNA+ATAC/multiome/RawData/BMMC-multiome-raw-RNA-counts.mtx
## save X to mtx
io.mmwrite('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-RNA-counts.mtx/matrix', BMMC_rna.X.T)
## save barcodes
with open('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-RNA-counts.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_rna.obs_names:
        f.write(item + '\n')      
## save features
with open('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-RNA-counts.mtx/features.tsv', 'w') as f:
    for item in BMMC_rna.var_names:
        f.write(item + '\n')
## gzip file
# !gzip RNA+ATAC/multiome/RawData/BMMC-multiome-raw-RNA-counts.mtx/*
## save metadata
BMMC_rna.obs.to_csv('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-RNA-counts.mtx/metadata.csv')

# %%
# Save ATAC
## makr dir
# !mkdir RNA+ATAC/multiome/RawData/BMMC-multiome-raw-ATAC-peaks.mtx
## save X to mtx
io.mmwrite('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-ATAC-peaks.mtx/matrix', BMMC_atac.X.T)
## save barcodes
with open('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-ATAC-peaks.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_atac.obs_names:
        f.write(item + '\n')      
## save features
with open('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-ATAC-peaks.mtx/features.tsv', 'w') as f:
    for item in BMMC_atac.var_names:
        f.write(item + '\n')
## gzip file
# !gzip RNA+ATAC/multiome/RawData/BMMC-multiome-raw-ATAC-peaks.mtx/*
## save metadata
BMMC_atac.obs.to_csv('RNA+ATAC/multiome/RawData/BMMC-multiome-raw-ATAC-peaks.mtx/metadata.csv')

# %%
# Save binarized ATAC
## makr dir
# !mkdir RNA+ATAC/multiome/RawData/BMMC-multiome-binarized-ATAC-peaks.mtx
## save X to mtx
io.mmwrite('RNA+ATAC/multiome/RawData/BMMC-multiome-binarized-ATAC-peaks.mtx/matrix', BMMC_atac.layers['binarized'].T)
## save barcodes
with open('RNA+ATAC/multiome/RawData/BMMC-multiome-binarized-ATAC-peaks.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_atac.obs_names:
        f.write(item + '\n')      
## save features
with open('RNA+ATAC/multiome/RawData/BMMC-multiome-binarized-ATAC-peaks.mtx/features.tsv', 'w') as f:
    for item in BMMC_atac.var_names:
        f.write(item + '\n')
## gzip file
# !gzip RNA+ATAC/multiome/RawData/BMMC-multiome-binarized-ATAC-peaks.mtx/*
## save metadata
BMMC_atac.obs.to_csv('RNA+ATAC/multiome/RawData/BMMC-multiome-binarized-ATAC-peaks.mtx/metadata.csv')

# %%

# %%
# print(BMMC_multi.obs.Site.value_counts())
# print(BMMC_multi.obs.Samplename.value_counts())
# print(BMMC_multi.obs.batch.value_counts())
# print(BMMC_multi.obs.Modality.value_counts())
# print(BMMC_multi.obs.cell_type.value_counts())

# print(BMMC_multi.obs.DonorAge.value_counts())
# print(BMMC_multi.obs.DonorBloodType.value_counts())
# print(BMMC_multi.obs.DonorBMI.value_counts())
# print(BMMC_multi.obs.DonorGender.value_counts())
# print(BMMC_multi.obs.DonorID.value_counts())
# print(BMMC_multi.obs.DonorNumber.value_counts())
# print(BMMC_multi.obs.DonorRace.value_counts())
# print(BMMC_multi.obs.DonorSmoker.value_counts())


# %%
import scanpy as sc
from scipy import io
BMMC_multi = sc.read_h5ad(BMMC_multi_h5ad)

# print(BMMC_multi.obs.DonorSmoker.value_counts())
BMMC_rna = BMMC_multi[:, BMMC_multi.var["feature_types"] == "GEX"]
BMMC_atac = BMMC_multi[:, BMMC_multi.var["feature_types"] == "ATAC"]
print(BMMC_rna.var.feature_types)
print(BMMC_atac.var.feature_types)

# %%
print(BMMC_rna.X.todense())
BMMC_rna.layers['processed'] = BMMC_rna.X.copy()
BMMC_rna.X = BMMC_rna.layers['counts'].copy()
print(BMMC_rna.X.todense())

# %%
print(BMMC_atac.X.todense())
BMMC_atac.layers['processed'] = BMMC_atac.X.copy()
BMMC_atac.X = BMMC_atac.layers['counts'].copy()
print(BMMC_atac.X.todense())

# %%
# !ls /Data/wangsg/BM/pipeline/results/BMMC/data_preprocess

# %%
output_path = "/Data/wangsg/BM/pipeline/results/BMMC/data_preprocess"

# save hd5
BMMC_rna.write_h5ad("{}/BMMC-raw-pair-RNA-counts.h5ad".format(output_path))

# save hd5
BMMC_atac.write_h5ad("{}/BMMC-raw-pair-ATAC-peaks.h5ad".format(output_path))

# %%
# Save RNA
## makr dir
# !mkdir BMMC-raw-pair-RNA-counts.mtx
## save X to mtx
io.mmwrite('BMMC-raw-pair-RNA-counts.mtx/matrix', BMMC_rna.X.T)
## save barcodes
with open('BMMC-raw-pair-RNA-counts.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_rna.obs_names:
        f.write(item + '\n')      
## save features
with open('BMMC-raw-pair-RNA-counts.mtx/features.tsv', 'w') as f:
    for item in BMMC_rna.var_names:
        f.write(item + '\n')
## gzip file
# !gzip BMMC-raw-pair-RNA-counts.mtx/*
## save metadata
BMMC_rna.obs.to_csv('BMMC-raw-pair-RNA-counts.mtx/metadata.csv')

# %%
# Save ATAC
## makr dir
# !mkdir BMMC-raw-pair-ATAC-peaks.mtx
## save X to mtx
io.mmwrite('BMMC-raw-pair-ATAC-peaks.mtx/matrix', BMMC_atac.X.T)
## save barcodes
with open('BMMC-raw-pair-ATAC-peaks.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_atac.obs_names:
        f.write(item + '\n')      
## save features
with open('BMMC-raw-pair-ATAC-peaks.mtx/features.tsv', 'w') as f:
    for item in BMMC_atac.var_names:
        f.write(item + '\n')
## gzip file
# !gzip BMMC-raw-pair-ATAC-peaks.mtx/*
## save metadata
BMMC_atac.obs.to_csv('BMMC-raw-pair-ATAC-peaks.mtx/metadata.csv')

# %%
# Save ATAC
BMMC_atac.X = BMMC_atac.layers['processed'].copy()

## makr dir
# !mkdir BMMC-raw-pair-ATAC-peaks-binarized.mtx
## save X to mtx
io.mmwrite('BMMC-raw-pair-ATAC-peaks-binarized.mtx/matrix', BMMC_atac.X.T)
## save barcodes
with open('BMMC-raw-pair-ATAC-peaks-binarized.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_atac.obs_names:
        f.write(item + '\n')      
## save features
with open('BMMC-raw-pair-ATAC-peaks-binarized.mtx/features.tsv', 'w') as f:
    for item in BMMC_atac.var_names:
        f.write(item + '\n')
## gzip file
# !gzip BMMC-raw-pair-ATAC-peaks-binarized.mtx/*
## save metadata
BMMC_atac.obs.to_csv('BMMC-raw-pair-ATAC-peaks-binarized.mtx/metadata.csv')

# %% [markdown]
# ## CITE-seq

# %%
BMMC_multi_h5ad = "/Data/wangsg/BM/pipeline/data/BMMC/RawData/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"

# %%
import scanpy as sc
from scipy import io
BMMC_multi = sc.read_h5ad(BMMC_multi_h5ad)

print(BMMC_multi.var.feature_types)
print(BMMC_multi.var.feature_types.value_counts())

BMMC_rna = BMMC_multi[:, BMMC_multi.var["feature_types"] == "GEX"]
BMMC_adt = BMMC_multi[:, BMMC_multi.var["feature_types"] == "ADT"]

print(BMMC_rna.var.feature_types)
print(BMMC_adt.var.feature_types)

# %%
print(BMMC_rna.X.todense())
BMMC_rna.layers['processed'] = BMMC_rna.X.copy()
BMMC_rna.X = BMMC_rna.layers['counts'].copy()
print(BMMC_rna.X.todense())

# %%
print(BMMC_adt.X.todense())
BMMC_adt.layers['processed'] = BMMC_adt.X.copy()
BMMC_adt.X = BMMC_adt.layers['counts'].copy()
print(BMMC_adt.X.todense())

# %%
output_path = "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq"

# save hd5
BMMC_rna.write_h5ad("{}/BMMC-raw-pair-RNA-counts.h5ad".format(output_path))
# save hd5
BMMC_adt.write_h5ad("{}/BMMC-raw-pair-ADT-counts.h5ad".format(output_path))

# %%
# Save RNA
## makr dir
# !mkdir RNA+ADT/CITE-seq/BMMC-raw-pair-RNA-counts.mtx
## save X to mtx
io.mmwrite('RNA+ADT/CITE-seq/BMMC-raw-pair-RNA-counts.mtx/matrix', BMMC_rna.X.T)
## save barcodes
with open('RNA+ADT/CITE-seq/BMMC-raw-pair-RNA-counts.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_rna.obs_names:
        f.write(item + '\n')      
## save features
with open('RNA+ADT/CITE-seq/BMMC-raw-pair-RNA-counts.mtx/features.tsv', 'w') as f:
    for item in BMMC_rna.var_names:
        f.write(item + '\n')
## gzip file
# !gzip RNA+ADT/CITE-seq/BMMC-raw-pair-RNA-counts.mtx/*
## save metadata
BMMC_rna.obs.to_csv('RNA+ADT/CITE-seq/BMMC-raw-pair-RNA-counts.mtx/metadata.csv')

# %%
# Save ADT
## makr dir
# !mkdir RNA+ADT/CITE-seq/BMMC-raw-pair-ADT-counts.mtx
## save X to mtx
io.mmwrite('RNA+ADT/CITE-seq/BMMC-raw-pair-ADT-counts.mtx/matrix', BMMC_adt.X.T)
## save barcodes
with open('RNA+ADT/CITE-seq/BMMC-raw-pair-ADT-counts.mtx/barcodes.tsv', 'w') as f:
    for item in BMMC_adt.obs_names:
        f.write(item + '\n')      
## save features
with open('RNA+ADT/CITE-seq/BMMC-raw-pair-ADT-counts.mtx/features.tsv', 'w') as f:
    for item in BMMC_adt.var_names:
        f.write(item + '\n')
## gzip file
# !gzip RNA+ADT/CITE-seq/BMMC-raw-pair-ADT-counts.mtx/*
## save metadata
BMMC_adt.obs.to_csv('RNA+ADT/CITE-seq/BMMC-raw-pair-ADT-counts.mtx/metadata.csv')


# %% [markdown]
# # Function: create pkl data of cobolt

# %%
def cobolt_pkl_data(input_path,
                   output_path,
                   dataset):
    from cobolt.utils import SingleData, MultiomicDataset
    from cobolt.model import Cobolt
    
    # Load Data
    rna = SingleData.from_file(path=input_path,
                               dataset_name=dataset["data_name"],
                               feature_name="GeneExpr",
                               count_file=dataset["gene_expression"],
                               barcode_file=dataset["gene_barcodes"],
                               feature_file=dataset["gene_names"])
    atac = SingleData.from_file(path=input_path,
                                dataset_name=dataset["data_name"],
                                feature_name="ChromAccess",
                                count_file=dataset["atac_expression"],
                                barcode_file=dataset["atac_barcodes"],
                                feature_file=dataset["atac_names"])
    # Filter Data
    # rna.filter_features(upper_quantile=0.99, lower_quantile=0.7)
    atac.filter_features(upper_quantile=0.99, lower_quantile=0.7)

    # Merge Data
    multi_dt = MultiomicDataset.from_singledata(rna, atac)

    # Save Data
    save_pickle(multi_dt,
                "{}/{}-{}-{}-multi-filtered.pkl".format(output_path, dataset["data_name"],
                                                "cobolt", dataset["task_type"]))


# %%
input_path = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
output_path = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset = json.load(open("/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/BMMC.json"))

# %%
cobolt_pkl_data(input_path, output_path, dataset)


# %% [markdown]
# # Function: create h5 data of multivi¶

# %%
def multivi_h5ad_data(input_path,
                     output_path,
                     dataset):
    # Load Data
    rna = sc.read_h5ad("{}/{}-{}-{}-RNA-counts.h5ad".format(input_path, 
                                                               dataset["data_name"], 
                                                               "raw", 
                                                               dataset["task_type"]))

    atac = sc.read_h5ad("{}/{}-{}-{}-ATAC-peaks.h5ad".format(input_path, 
                                                               dataset["data_name"], 
                                                               "raw", 
                                                               dataset["task_type"]))


    # Filter Data
    ## RNA
    sc.pp.filter_genes(rna, min_cells=int(rna.shape[0] * 0.001))
    # sc.pp.filter_cells(rna, min_genes=3)
    ## ATAC
    sc.pp.filter_genes(atac, min_cells=int(rna.shape[0] * 0.001))
    sc.pp.filter_genes(atac, max_cells=int(rna.shape[0] * 0.1))
    # sc.pp.filter_cells(atac, min_genes=3)

    # Merge Data
    multi = anndata.concat([rna, atac], axis=1, join="outer")
    multi.var_names_make_unique()

    # Add Metadata
    multi.var['modality'] = np.repeat(["Gene Expression", "Peaks"], [rna.shape[1], atac.shape[1]], axis=0)

    # Save Data
    del multi.raw
    multi.write_h5ad("{}/{}-{}-{}-multi-filtered.h5ad".format(output_path, dataset["data_name"],
                                                "multivi", dataset["task_type"]))

# %%
input_path = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
output_path = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset = json.load(open("/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/BMMC.json"))

# %%
multivi_h5ad_data(input_path, output_path, dataset)


# %% [markdown]
# # Function: create pkl of scMVP

# %%
def scMVP_pkl_data(input_path,
                   output_path, 
                   dataset, 
                   annotation=None, 
                   annotation_column=None):
    
    from scMVP.dataset import LoadData, GeneExpressionDataset, CellMeasurement

    # Load Data
    ## prepare dataset
    input_path = "{}/".format(input_path)
    output_path = "{}/".format(output_path)
    dataset_sub = {
        k: dataset.get(k, None)
        for k in ("gene_names", "gene_barcodes", "gene_expression",
                  "atac_names", "atac_barcodes", "atac_expression")
    }
    ## multi
    if annotation:
        cell_embeddings = pd.read_csv(input_path + annotation,
                                      sep="\t",
                                      index_col=None).iloc[:,
                                                           annotation_column]
        scmvp_multi_data = LoadData(dataset=dataset_sub,
                                   data_path=input_path,
                                   dense=False,
                                   gzipped=False,
                                   atac_threshold=0.001,
                                   cell_threshold=1,
                                   cell_meta=cell_embeddings)
    else:
        scmvp_multi_data = LoadData(dataset=dataset_sub,
                                   data_path=input_path,
                                   dense=False,
                                   gzipped=False,
                                   atac_threshold=0.001,
                                   cell_threshold=1)

    # Save Data
    save_pickle(scmvp_multi_data,
                "{}/{}-{}-{}-multi-filtered.pkl".format(output_path, dataset["data_name"],
                                                "scMVP", dataset["task_type"]))


# %%
input_path = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
output_path = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset = json.load(open("/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/BMMC.json"))

# %%
scMVP_pkl_data(input_path, output_path, dataset)

# %%
annotation=None
annotation_column=None

# %%
from scMVP.dataset import LoadData, GeneExpressionDataset, CellMeasurement

# Load Data
## prepare dataset
input_path = "{}/".format(input_path)
output_path = "{}/".format(output_path)
dataset_sub = {
    k: dataset.get(k, None)
    for k in ("gene_names", "gene_barcodes", "gene_expression",
              "atac_names", "atac_barcodes", "atac_expression")
}
## multi
if annotation:
    cell_embeddings = pd.read_csv(input_path + annotation,
                                  sep="\t",
                                  index_col=None).iloc[:,
                                                       annotation_column]
    scmvp_multi_data = LoadData(dataset=dataset_sub,
                               data_path=input_path,
                               dense=False,
                               gzipped=False,
                               atac_threshold=0.001,
                               cell_threshold=1,
                               cell_meta=cell_embeddings)
else:
    scmvp_multi_data = LoadData(dataset=dataset_sub,
                               data_path=input_path,
                               dense=False,
                               gzipped=False,
                               atac_threshold=0.001,
                               cell_threshold=1)

# # Save Data
# save_pickle(scmvp_multi_data,
#             "{}/{}-{}-{}-multi-filtered.pkl".format(output_path, dataset["data_name"],
#                                             "scMVP", dataset["task_type"]))

# %%

# %%

# %%

# %%

# %%
