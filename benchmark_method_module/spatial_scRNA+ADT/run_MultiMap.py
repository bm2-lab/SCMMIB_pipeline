import argparse
import subprocess
import os
# 创建命令行解析器
parser = argparse.ArgumentParser(description='noimputation_mosaic_RNA_ADT')
# 添加布尔类型参数
parser.add_argument('--h5ad_path', type=str, help='h5ad_path')
parser.add_argument('--index', type=str, help='index')
parser.add_argument('--gpu_index',  type=str, help='gpu_index')

args = parser.parse_args()

h5ad_path = args.h5ad_path
index = args.index
gpu_index = int(args.gpu_index)


method = 'multimap'
target_folder = h5ad_path.replace('/home/wsg/BM/data','/home/sirm/project/SCMMIB') + '/'+method
if not os.path.isdir(target_folder):
    cmd = f'mkdir -p {target_folder}'
    subprocess.run(cmd,shell=True)
latent_embed_path =os.path.join(target_folder,f'RUN_{index}',f'{method}_latent.csv')







import scanpy as sc
import anndata
import MultiMAP
import pandas as pd
import muon
sc.settings.set_figure_params(dpi=80)
from glob import glob

RNA_data_path = glob(f'{h5ad_path}/*RNA-counts.h5ad')[0]
ADT_data_path = glob(f'{h5ad_path}/*ADT-counts.h5ad')[0]
# read data
#file_fold = '/home/shaliu_fu/joint_bench/demo_notebook/input/spatial/human_lymph_node/' #please replace 'file_fold' with the download path
adata_omics1 = sc.read_h5ad(RNA_data_path)
adata_omics2 = sc.read_h5ad(ADT_data_path)
# adata_omics1 = sc.read_h5ad(f'{file_fold}/adata_RNA.h5ad')
# adata_omics2 = sc.read_h5ad(f'{file_fold}/adata_ADT.h5ad')

adata_omics1.var_names_make_unique()
adata_omics2.var_names_make_unique()
rna = adata_omics1
adt = adata_omics2
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)


rna_pca = rna.copy()
sc.pp.pca(rna_pca,n_comps=10)
rna.obsm['X_pca'] = rna_pca.obsm['X_pca'].copy()

adt_pca = adt.copy()
muon.prot.pp.clr(adt_pca)

sc.pp.pca(adt_pca,n_comps=10) #维度必须同RNA一样，并且小于蛋白数量
adt.obsm['X_pca'] = adt_pca.obsm['X_pca'].copy()
adata1 = MultiMAP.matrix.MultiMAP([rna.obsm['X_pca'],adt.obsm['X_pca']])[2] # rna和adt共同映射
adata2 = MultiMAP.matrix.MultiMAP([adata1,rna.obsm['spatial']])[2] # 然后将2维空间和spatial共同映射


n = len(rna.obs.index)
tmp1 = adata2[:n,:]
tmp2 = adata2[n:2*n,:]
tmp3 = adata2[2*n:,:]
import numpy as np
latent = pd.DataFrame(np.c_[tmp1,tmp2,tmp3],index=rna.obs.index)

latent.to_csv(latent_embed_path)
