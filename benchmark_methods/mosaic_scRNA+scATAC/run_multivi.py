
import os
import sys 
import json
import time
import pickle

import scvi
import torch
import random
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import anndata as ad

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
torch.set_float32_matmul_precision("high")



def multivi_module(input_path,
                 output_path,
                 config,
                 impute=True):


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
        

    del RNA_data.raw
    del ATAC_data.raw

    RNA_data.var['modality'] =  "Gene Expression"
    ATAC_data.var['modality'] =  "Peaks"


    tmp = RNA_data.obs
    adata_RNA_ATAC = sc.concat([RNA_data,ATAC_data],axis=1)
    adata_RNA_ATAC.obs =  adata_RNA_ATAC.obs.join(tmp)
    adata = adata_RNA_ATAC
    adata_rna = adata[adata.obs['batch']==unpaired_dataset, adata.var.modality == "Gene Expression"].copy()
    adata_paired = adata[adata.obs['batch']==paired_dataset].copy()
    adata_atac = adata[adata.obs['batch']==unpaired_dataset, adata.var.modality == "Peaks"].copy()
    adata_rna.obs.index = adata_rna.obs.index.map(lambda x:x+'_rna')
    adata_atac.obs.index = adata_atac.obs.index.map(lambda x:x+'_atac')
    adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired, adata_rna, adata_atac)
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
    print(adata_mvi.shape)
    sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
    print(adata_mvi.shape)
    scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality")
    model = scvi.model.MULTIVI(
        adata_mvi,
        n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
        n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
    )
    model.view_anndata_setup()

    from nvitop import Device
    devices = Device.all()  # or Device.cuda.all()
    memory_free = [device.memory_free() for device in devices]
    gpu_index = memory_free.index(max(memory_free))
    print(gpu_index)
    print(memory_free)
    torch.cuda.set_device(gpu_index)
    torch.set_num_threads(5)


    model.train(use_gpu = gpu_index)
    latent = model.get_latent_representation()
    index_change = adata_mvi.obs.index.map(lambda x :"_".join(x.split('_')[:-1]))
    latent = pd.DataFrame(latent,index=index_change)
    # latent = latent.round(3)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-multivi-mosaic-latent.csv"))
    if impute:
        imputation_file_path = os.path.join(output_path, config["output_prefix"] + "-multivi-mosaic-imputation.csv")
        imputation_file_path_rna = imputation_file_path.replace('.csv', '_rna.csv')
        imputation_file_path_atac = imputation_file_path.replace('.csv', '_atac.csv')
        imputed_expression = model.get_normalized_expression()
        imputed_expression.index = imputed_expression.index.map(lambda x :"_".join(x.split('_')[:-1]))
        rna_index = imputed_expression.index.map(lambda x:x.split('_')[-1]=='atac')
        imputed_rna = imputed_expression[rna_index ]
        # imputed_rna = imputed_rna.round(3)
        imputed_rna.to_csv(imputation_file_path_rna)
        imputed_accessibility = model.get_accessibility_estimates()
        imputed_accessibility.index = imputed_accessibility.index.map(lambda x :"_".join(x.split('_')[:-1]))
        atac_index = imputed_accessibility.index.map(lambda x:x.split('_')[-1]=='rna')
        imputed_atac = imputed_accessibility[atac_index ]
        # imputed_atac = imputed_atac.round(3)
        imputed_atac.to_csv(imputation_file_path_atac)

    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-multivi-mosaic-gpu_memory.csv'), header=["gpu_memory"])

multivi_module(input_path = sys.argv[1],
             output_path = sys.argv[2],
             config = json.load(open(sys.argv[3])),
             impute = sys.argv[4]
            )