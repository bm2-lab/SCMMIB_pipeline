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
#     display_name: python (metrics)
#     language: python
#     name: metrics
# ---

# %%
import os
import sys 
import json
import pickle

import torch
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from scipy import io

# %% [markdown]
# ## lymph_node

# %%
lymph_node_PATH_1 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset11_Human_Lymph_Node_A1"
lymph_node_PATH_2 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset12_Human_Lymph_Node_D1"

# %%
lymph_node_rna_1 = sc.read_h5ad("{}/adata_RNA.h5ad".format(lymph_node_PATH_1))
lymph_node_rna_2 = sc.read_h5ad("{}/adata_RNA.h5ad".format(lymph_node_PATH_2))

# %%
print(lymph_node_rna_1)
print(lymph_node_rna_2)

# %%
adata = lymph_node_rna_1
adata.var['gene_name'] = adata.var_names
adata.var_names = adata.var['gene_ids']
lymph_node_rna_1 = adata

adata = lymph_node_rna_2
adata.var['gene_name'] = adata.var_names
adata.var_names = adata.var['gene_ids']
lymph_node_rna_2 = adata


# %%
lymph_node_rna_1.obs_names = [name + "_1" for name in lymph_node_rna_1.obs_names]
lymph_node_rna_2.obs_names = [name + "_2" for name in lymph_node_rna_2.obs_names]

# %%
print(lymph_node_rna_1)
print(lymph_node_rna_2)

# %%
lymph_node_rna_1.obs['batch'] = 'A1'
lymph_node_rna_2.obs['batch'] = 'D1'

# %%
lymph_node_rna = ad.concat([lymph_node_rna_1, lymph_node_rna_2], axis=0)
lymph_node_rna.var = lymph_node_rna_1.var

# %%
lymph_node_rna

# %%
lymph_node_rna.var

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node"

lymph_node_rna.write_h5ad("{}/lymph_node-CITE_seq_RNA-counts.h5ad".format(output_path))

# %%

# %%
lymph_node_adt_1 = sc.read_h5ad("{}/adata_ADT.h5ad".format(lymph_node_PATH_1))
lymph_node_adt_2 = sc.read_h5ad("{}/adata_ADT.h5ad".format(lymph_node_PATH_2))

# %%
print(lymph_node_adt_1.shape)
print(lymph_node_adt_2.shape)

# %%
lymph_node_adt_1.X.toarray()

# %%
adata = lymph_node_adt_1
adata.var['adt_name'] = adata.var_names
adata.var_names = adata.var['gene_ids']
lymph_node_adt_1 = adata

adata = lymph_node_adt_2
adata.var['adt_name'] = adata.var_names
adata.var_names = adata.var['gene_ids']
lymph_node_adt_2 = adata

# %%
lymph_node_adt_1.obs_names = [name + "_1" for name in lymph_node_adt_1.obs_names]
lymph_node_adt_2.obs_names = [name + "_2" for name in lymph_node_adt_2.obs_names]

# %%
lymph_node_adt_1.obs['batch'] = 'A1'
lymph_node_adt_2.obs['batch'] = 'D1'

# %%
lymph_node_adt = ad.concat([lymph_node_adt_1, lymph_node_adt_2], axis=0)

# %%
lymph_node_adt.var = lymph_node_adt_1.var
lymph_node_adt.var_names = lymph_node_adt.var['adt_name']

# %%
gene_tab = pd.read_csv("/home/wsg/BM/pipeline/config/gencode.v41.gtf.uniq.txt", delimiter=' ', header=None)
gene_tab['gene_ids'] = gene_tab[0].str.split('.').str[0]
gene_tab

adt_new_var = pd.merge(lymph_node_adt.var, gene_tab, left_on='gene_ids', right_on=1, how='left')
adt_new_var = adt_new_var.drop_duplicates('gene_ids')
adt_new_var = adt_new_var.drop(0, axis=1)
adt_new_var = adt_new_var.drop(1, axis=1)
adt_new_var.astype(str)
adt_new_var = adt_new_var.set_index('adt_name')

lymph_node_adt.var = adt_new_var
lymph_node_adt = lymph_node_adt[:, ~lymph_node_adt.var['gene_ids_y'].isna()]
lymph_node_adt.var = lymph_node_adt.var.set_index('gene_ids_y')

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node"

lymph_node_adt.write_h5ad("{}/lymph_node-CITE_seq_ADT-counts.h5ad".format(output_path))

# %%

# %%

# %%

# %% [markdown]
# ### human_lymph_node_A1

# %%
lymph_node_PATH = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/human_lymph_node"

# %%
lymph_node_rna = sc.read_h5ad("{}/adata_RNA.h5ad".format(lymph_node_PATH))
lymph_node_adt = sc.read_h5ad("{}/adata_ADT.h5ad".format(lymph_node_PATH_1))

# %%
metadata = pd.read_csv("{}/annotation.csv".format(lymph_node_PATH))
metadata.columns = ["barcode", "cell_type"]
metadata.index = metadata["barcode"]

# %%
lymph_node_rna.obs = metadata

# %%
lymph_node_adt.obs = metadata

# %%
lymph_node_rna.var

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1"

lymph_node_rna.write_h5ad("{}/lymph_node-CITE_seq_RNA-counts.h5ad".format(output_path))
lymph_node_adt.write_h5ad("{}/lymph_node-CITE_seq_ADT-counts.h5ad".format(output_path))

# %%

# %%

# %% [markdown]
# ## thymus

# %%
thymus_PATH_1 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset3_Mouse_Thymus1"
thymus_PATH_2 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset4_Mouse_Thymus2"
thymus_PATH_3 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset5_Mouse_Thymus3"
thymus_PATH_4 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset6_Mouse_Thymus4"

# %%
thymus_rna_1 = sc.read_h5ad("{}/adata_RNA.h5ad".format(thymus_PATH_1))
thymus_rna_2 = sc.read_h5ad("{}/adata_RNA.h5ad".format(thymus_PATH_2))
thymus_rna_3 = sc.read_h5ad("{}/adata_RNA.h5ad".format(thymus_PATH_3))
thymus_rna_4 = sc.read_h5ad("{}/adata_RNA.h5ad".format(thymus_PATH_4))

# %%
print(thymus_rna_1)
print(thymus_rna_2)
print(thymus_rna_3)
print(thymus_rna_4)

# %%
thymus_rna_1.obs_names = [name + "_1" for name in thymus_rna_1.obs_names]
thymus_rna_2.obs_names = [name + "_2" for name in thymus_rna_2.obs_names]
thymus_rna_2.obs_names = [name + "_3" for name in thymus_rna_3.obs_names]
thymus_rna_4.obs_names = [name + "_4" for name in thymus_rna_4.obs_names]

# %%
thymus_rna_1.obs['batch'] = 'thymus_1'
thymus_rna_2.obs['batch'] = 'thymus_2'
thymus_rna_3.obs['batch'] = 'thymus_3'
thymus_rna_4.obs['batch'] = 'thymus_4'

# %%
all_vars = pd.unique(thymus_rna_1.var_names.tolist() + thymus_rna_2.var_names.tolist() + 
                     thymus_rna_3.var_names.tolist() + thymus_rna_4.var_names.tolist())
len(all_vars)


# %%
def align_anndata_vars(adata, all_vars):
    from scipy.sparse import lil_matrix
    aligned_data = lil_matrix((adata.n_obs, len(all_vars)))
    
    var_index_map = {var: idx for idx, var in enumerate(all_vars)}
    
    for var_name in adata.var_names:
        aligned_data[:, var_index_map[var_name]] = adata[:, var_name].X
    
    return ad.AnnData(X=aligned_data, obs=adata.obs, var=pd.DataFrame(index=all_vars))

thymus_rna_1_aligned = align_anndata_vars(thymus_rna_1, all_vars)
thymus_rna_2_aligned = align_anndata_vars(thymus_rna_2, all_vars)
thymus_rna_3_aligned = align_anndata_vars(thymus_rna_3, all_vars)
thymus_rna_4_aligned = align_anndata_vars(thymus_rna_4, all_vars)

# %%
thymus_rna_1_aligned.obsm['spatial'] = thymus_rna_1.obsm['spatial']
thymus_rna_2_aligned.obsm['spatial'] = thymus_rna_2.obsm['spatial']
thymus_rna_3_aligned.obsm['spatial'] = thymus_rna_3.obsm['spatial']
thymus_rna_4_aligned.obsm['spatial'] = thymus_rna_4.obsm['spatial']

# %%
thymus_rna = ad.concat([thymus_rna_1_aligned, thymus_rna_2_aligned, thymus_rna_3_aligned, thymus_rna_4_aligned], axis=0)

# %%
thymus_rna.obs

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/thymus"
thymus_rna.write_h5ad("{}/thymus-CITE_seq-raw-RNA-counts.h5ad".format(output_path))

# %%

# %%
thymus_adt_1 = sc.read_h5ad("{}/adata_ADT.h5ad".format(thymus_PATH_1))
thymus_adt_2 = sc.read_h5ad("{}/adata_ADT.h5ad".format(thymus_PATH_2))
thymus_adt_3 = sc.read_h5ad("{}/adata_ADT.h5ad".format(thymus_PATH_3))
thymus_adt_4 = sc.read_h5ad("{}/adata_ADT.h5ad".format(thymus_PATH_4))

# %%
print(thymus_adt_1)
print(thymus_adt_2)
print(thymus_adt_3)
print(thymus_adt_4)

# %%
thymus_adt_1.obs_names = [name + "_1" for name in thymus_adt_1.obs_names]
thymus_adt_2.obs_names = [name + "_2" for name in thymus_adt_2.obs_names]
thymus_adt_3.obs_names = [name + "_3" for name in thymus_adt_3.obs_names]
thymus_adt_4.obs_names = [name + "_4" for name in thymus_adt_4.obs_names]

# %%
thymus_adt_1.obs['batch'] = 'thymus_1'
thymus_adt_2.obs['batch'] = 'thymus_2'
thymus_adt_3.obs['batch'] = 'thymus_3'
thymus_adt_4.obs['batch'] = 'thymus_4'


# %%
def rename(adata):    
    import re
    new_var_names = []
    for name in adata.var_names:
        parts = re.split('[-_]', name)
        if '-'.join(parts) == 'Ms-Hu-CD11b':
            parts = re.split('[-_]', 'mouse-human-CD11b')
        elif parts[0] == 'Mouse':
            parts[0] = 'mouse'
        elif parts[0] == 'Rat':
            parts[0] = 'rat'
        elif parts[0] == 'Human':
            parts[0] = 'human'
        new_name = '-'.join(parts)
        new_var_names.append(new_name)
    adata.var_names = new_var_names
    return adata


# %%
print(thymus_adt_1.shape)
print(thymus_adt_2.shape)
print(thymus_adt_3.shape)
print(thymus_adt_4.shape)

# %%
print(thymus_adt_1.var_names)
print(thymus_adt_2.var_names)
print(thymus_adt_3.var_names)
print(thymus_adt_4.var_names)

# %%
thymus_adt_1 = rename(thymus_adt_1)
thymus_adt_2 = rename(thymus_adt_2)
thymus_adt_3 = rename(thymus_adt_3)
thymus_adt_4 = rename(thymus_adt_4)

# %%
print(thymus_adt_1.var_names)
print(thymus_adt_2.var_names)
print(thymus_adt_3.var_names)
print(thymus_adt_4.var_names)

# %%
all_vars = pd.unique(thymus_adt_1.var_names.tolist() + thymus_adt_2.var_names.tolist() + 
                     thymus_adt_3.var_names.tolist() + thymus_adt_4.var_names.tolist())
len(all_vars)


# %%
def align_anndata_vars(adata, all_vars):
    from scipy.sparse import lil_matrix
    aligned_data = lil_matrix((adata.n_obs, len(all_vars)))
    
    var_index_map = {var: idx for idx, var in enumerate(all_vars)}
    
    for var_name in adata.var_names:
        aligned_data[:, var_index_map[var_name]] = adata[:, var_name].X
    
    return ad.AnnData(X=aligned_data, obs=adata.obs, var=pd.DataFrame(index=all_vars))

thymus_adt_1_aligned = align_anndata_vars(thymus_adt_1, all_vars)
thymus_adt_2_aligned = align_anndata_vars(thymus_adt_2, all_vars)
thymus_adt_3_aligned = align_anndata_vars(thymus_adt_3, all_vars)
thymus_adt_4_aligned = align_anndata_vars(thymus_adt_4, all_vars)

# %%
thymus_adt_1_aligned.obsm['spatial'] = thymus_adt_1.obsm['spatial']
thymus_adt_2_aligned.obsm['spatial'] = thymus_adt_2.obsm['spatial']
thymus_adt_3_aligned.obsm['spatial'] = thymus_adt_3.obsm['spatial']
thymus_adt_4_aligned.obsm['spatial'] = thymus_adt_4.obsm['spatial']

# %%
thymus_adt = ad.concat([thymus_adt_1_aligned, thymus_adt_2_aligned, thymus_adt_3_aligned, thymus_adt_4_aligned], axis=0)

# %%
thymus_adt.var.index = thymus_adt.var.index.str.replace('human-', '', regex=True)
thymus_adt.var.index = thymus_adt.var.index.str.replace('mouse-', '', regex=True)
thymus_adt.var.index = thymus_adt.var.index.str.replace('rat-', '', regex=True)
thymus_adt.var['adt_name'] = thymus_adt.var.index.str.split('-').str[0]

# %%
gene_tab = pd.read_csv('/home/wsg/BM/pipeline/config/gene_adt_tab.tsv', sep=' ', header=0)
gene_tab['gene_name'] = gene_tab['Symbol'].str.capitalize()

# %%
adt_new_var = pd.merge(thymus_adt.var, gene_tab, left_on='adt_name', right_on='ADT', how='left')
adt_new_var.astype(str)

# %%
thymus_adt.var = adt_new_var
thymus_adt = thymus_adt[:, ~thymus_adt.var['gene_name'].isna()]
thymus_adt.var = thymus_adt.var.set_index('gene_name')

# %%
thymus_adt.obs

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/thymus"
thymus_adt.write_h5ad("{}/thymus-CITE_seq-raw-ADT-counts.h5ad".format(output_path))

# %%

# %%

# %% [markdown]
# ## spleen

# %%
spleen_PATH_1 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset1_Mouse_Spleen1"
spleen_PATH_2 = "/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/Dataset2_Mouse_Spleen2"

# %%
spleen_rna_1 = sc.read_h5ad("{}/adata_RNA.h5ad".format(spleen_PATH_1))
spleen_rna_2 = sc.read_h5ad("{}/adata_RNA.h5ad".format(spleen_PATH_2))

# %%
print(spleen_rna_1)
print(spleen_rna_2)

# %%
adata = spleen_rna_1
adata.var['gene_name'] = adata.var_names
adata.var_names = adata.var['gene_ids']
spleen_rna_1 = adata

adata = spleen_rna_2
adata.var['gene_name'] = adata.var_names
adata.var_names = adata.var['gene_ids']
spleen_rna_2 = adata


# %%
spleen_rna_1.obs_names = [name + "_1" for name in spleen_rna_1.obs_names]
spleen_rna_2.obs_names = [name + "_2" for name in spleen_rna_2.obs_names]

# %%
spleen_rna_1.obs['batch'] = 'spleen_1'
spleen_rna_2.obs['batch'] = 'spleen_2'

# %%
spleen_rna = ad.concat([spleen_rna_1, spleen_rna_2], axis=0)
spleen_rna.var = spleen_rna_1.var

# %%
spleen_rna.var = spleen_rna.var.astype(str)

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen"

spleen_rna.write_h5ad("{}/spleen-CITE_seq-raw-RNA-counts.h5ad".format(output_path))

# %%

# %%

# %%
spleen_adt_1.obs['batch'] = 'spleen_1'
spleen_adt_2.obs['batch'] = 'spleen_2'

# %%
spleen_adt_1 = sc.read_h5ad("{}/adata_ADT.h5ad".format(spleen_PATH_1))
spleen_adt_2 = sc.read_h5ad("{}/adata_ADT.h5ad".format(spleen_PATH_2))

# %%
print(spleen_adt_1)
print(spleen_adt_2)

# %%
spleen_adt_1.obs_names = [name + "_1" for name in spleen_adt_1.obs_names]
spleen_adt_2.obs_names = [name + "_2" for name in spleen_adt_2.obs_names]

# %%
spleen_adt = ad.concat([spleen_adt_1, spleen_adt_2], axis=0)
spleen_adt.var = spleen_adt_1.var

# %%
spleen_adt.var['adt_name'] = spleen_adt.var.index

# %%
gene_tab = pd.read_csv('/home/wsg/BM/pipeline/config/gene_adt_tab.tsv', sep=' ', header=0)
gene_tab['gene_name'] = gene_tab['Symbol'].str.capitalize()

# %%
adt_new_var = pd.merge(spleen_adt.var, gene_tab, left_on='adt_name', right_on='ADT', how='left')
adt_new_var.astype(str)

spleen_adt.var = adt_new_var
spleen_adt = spleen_adt[:, ~spleen_adt.var['gene_name'].isna()]

adt_new_var = pd.merge(spleen_adt.var, spleen_rna.var, left_on='gene_name', right_on='gene_name', how='left')
adt_new_var.astype(str)

spleen_adt.var = adt_new_var
spleen_adt.var = spleen_adt.var.set_index('gene_ids')

# %%

# %%
output_path = "/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen"
spleen_adt.write_h5ad("{}/spleen-CITE_seq-raw-ADT-counts.h5ad".format(output_path))

# %%

# %%

# %%
