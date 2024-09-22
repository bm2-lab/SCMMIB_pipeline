import warnings
warnings.filterwarnings('ignore')
import torch
import scanpy as sc
import numpy as np
import pandas as pd
import muon
import multigrate as mtg
import json
import sys
import anndata as ad
import os

##### 这里替换原来的输入
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# print(input_path ,output_path,config)






RNA_data = ad.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
ATAC_data = ad.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

if 'batch' in RNA_data.obs.columns:
    data_names = RNA_data.obs['batch'].unique()
    if 's1d1' in data_names:
        paired_dataset= "s1d1"
        unpaired_dataset= "s3d10"
    else:
        paired_dataset= "s3d6"
        unpaired_dataset= "s2d1"
        
        
        

preprocess_RNA=True
preprocess_ADT=True
preprocess=True
imputation=True

if 'data_size' in RNA_data.obs.columns:
    data_names = RNA_data.obs['data_size'].unique()
    if 'c5k' in data_names:
        unpaired_dataset = "c5k"
        paired_dataset = data_names[data_names != "c5k"][0]
    elif 'c5k_1' in data_names:
        unpaired_dataset = "c5k_1"
        paired_dataset = "c5k_2"
    else:
        data_size = data_names[0].split("_")[0]
        unpaired_dataset = data_size + "_1"
        paired_dataset = data_size + "_2"
    
    RNA_data.obs['batch'] =  RNA_data.obs['data_size']
    ATAC_data.obs['batch'] = ATAC_data.obs['data_size']
    


####



def RNA_process(rna,preprocess=True):
    if preprocess:
        rna.layers['counts'] = rna.X.copy()
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
        sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor='seurat', batch_key='batch')
        rna_hvg = rna[:, rna.var.highly_variable].copy()
    else:
        rna_hvg = rna[:, rna.var.highly_variable].copy()
        rna_hvg = rna
    return rna_hvg

def ATAC_process(atac,preprocess=True):
    if preprocess:
        atac.layers['counts'] = atac.X.copy()
        muon.atac.pp.tfidf(atac, scale_factor=1e4)
        atac.layers['tf-idf'] = atac.X.copy()
        atac.X = atac.layers['counts'].copy()
        sc.pp.normalize_total(atac, target_sum=1e4)
        sc.pp.log1p(atac)
        sc.pp.highly_variable_genes(atac, n_top_genes=30000, flavor='seurat', batch_key='batch')
        atac.layers['log-norm'] = atac.X.copy()
        atac_hvf = atac[:, atac.var.highly_variable].copy()
    else:
        atac_hvf = atac[:, atac.var.highly_variable].copy()
        atac_hvf = atac
    return atac_hvf

RNA_data = RNA_process(RNA_data,preprocess)
ATAC_data = ATAC_process(ATAC_data) # ATAC没有preprocess
hvg_number = len(RNA_data.var.index)
paired_rna = RNA_data[RNA_data.obs['batch']==paired_dataset].copy()
paired_atac = ATAC_data[ATAC_data.obs['batch']==paired_dataset].copy()
unpaired_rna = RNA_data[RNA_data.obs['batch']==unpaired_dataset].copy()
unpaired_atac = ATAC_data[ATAC_data.obs['batch']==unpaired_dataset].copy()
unpaired_atac.obs.index = ["{}_atac".format(i) for i in unpaired_atac.obs.index.values]
unpaired_rna.obs.index = ["{}_rna".format(i) for i in unpaired_rna.obs.index.values]
print(RNA_data.obs['batch'].value_counts(),paired_rna,unpaired_rna, unpaired_atac)
if preprocess:
    adata = mtg.data.organize_multiome_anndatas(
        adatas = [[paired_rna,unpaired_rna,None], [paired_atac,None, unpaired_atac]],           # a list of anndata objects per modality, RNA-seq always goes first
        layers = [['counts','counts',None], ['log-norm',None,'log-norm']],    # if need to use data from .layers, if None use .X
    )
else:
    adata = mtg.data.organize_multiome_anndatas(
        adatas = [[paired_rna,unpaired_rna,None], [paired_atac,None, unpaired_atac]],           # a list of anndata objects per modality, RNA-seq always goes first
        layers = [[None,None,None], ['log-norm',None,'log-norm']],    # if need to use data from .layers, if None use .X
    )
mtg.model.MultiVAE.setup_anndata(
    adata,
    rna_indices_end=hvg_number, # how many features in the rna-seq modality
    categorical_covariate_keys=["batch"]
)
model = mtg.model.MultiVAE(
    adata,
    losses=['nb', 'mse'],
)

########################### 设置GPU
# https://multigrate.readthedocs.io/en/latest/notebooks/paired_integration_multiome.html
from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
print(gpu_index)
print(memory_free)
torch.cuda.set_device(gpu_index)
torch.set_num_threads(5)
########################### 

model.train()
model.get_latent_representation()
latent = pd.DataFrame(data=adata.obsm['latent'],
                          index=adata.obs_names)
# latent = latent.round(3)
latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-multigrate-mosaic-latent.csv"))

# 统计GPU
pid= os.getpid()        
gpu_memory = pd.Series(dtype='str')

devices = Device.all()
for device in devices:
    processes = device.processes()    
    if pid in processes.keys():
        p=processes[pid]
        gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()


gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-multigrate-mosaic-gpu_memory.csv'), header=["gpu_memory"])