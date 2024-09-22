import STAligner

import os
os.environ['R_HOME'] = "/home/shaliu_fu/miniconda3/envs/joint_bench/lib/R"
os.environ['R_USER'] = "/home/shaliu_fu/miniconda3/envs/joint_bench/lib/python3.8/site-packages/rpy2"
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg

import torch
used_device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')


section_ids = ['Dataset11_Human_Lymph_Node_A1', 'Dataset12_Human_Lymph_Node_D1']
section_id= section_ids[0]

adata = sc.read_h5ad(os.path.join("/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/" + section_id + "/adata_RNA.h5ad"))

Batch_list = []
adj_list = []

section_ids = ['Dataset11_Human_Lymph_Node_A1', 'Dataset12_Human_Lymph_Node_D1']
# for section_id in section_ids:
for index in range(len(section_ids)):
    section_id = section_ids[index]
    print(section_id)
    adata = sc.read_h5ad(os.path.join("/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/Data_SpatialGlue/" + section_id + "/adata_RNA.h5ad"))
    adata.layers['count'] = adata.X 

    # make spot name unique
    adata.obs_names = [x + '_' + str(index+1) for x in adata.obs_names]
    adata.obs['batch'] = index+1

    # make sure var name unique
    adata.var_names_make_unique()

    STAligner.Cal_Spatial_Net(adata, rad_cutoff=1.5)

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000) #ensure enough common HVGs in the combined matrix
    adata = adata[:, adata.var['highly_variable']]

    adj_list.append(adata.uns['adj'])
    Batch_list.append(adata)

adata_concat = ad.concat(Batch_list, label="slice_name", keys=section_ids)
adata_concat.obs["batch_name"] = adata_concat.obs["slice_name"].astype('category')

adata_concat = STAligner.train_STAligner_subgraph(adata_concat, verbose=True, knn_neigh = 100, n_epochs = 600, iter_comb = [(0,1)],
                                                        Batch_list=Batch_list, device=used_device)

sc.pp.neighbors(adata_concat, use_rep='STAligner', random_state=666)
sc.tl.louvain(adata_concat, random_state=666, key_added="louvain", resolution=0.2)

out_tab = pd.DataFrame(adata_concat.obsm['STAligner'],index=adata_concat.obs_names)
out_tab.to_csv("../output/staligner/lymph_node_latent.csv")