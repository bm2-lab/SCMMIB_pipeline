import time
start = time.time()
import scvi
import os,sys,json
import anndata
# from scvi.model.MULTIVI import setup_anndata, synthetic_iid, transfer_anndata_setup
import pandas as pd
from scipy.sparse import csr_matrix
import scipy.io as sp_io
from anndata import AnnData
import numpy as np
import scipy.sparse as sparse
import scanpy as sc
import torch
import matplotlib.pyplot as plt



from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
print(gpu_index)
print(memory_free)
torch.cuda.set_device(gpu_index)
torch.set_num_threads(5)


def totalVI_module(input_path,
                    output_path,
                    config,
                    impute=True):

    RNA_data = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    ADT_data = sc.read_h5ad(os.path.join(input_path, config['adt_h5ad_filename']))

    if 'batch' in RNA_data.obs.columns:
        data_names = RNA_data.obs['batch'].unique()
        if 's1d1' in data_names:
            paired_dataset= "s1d1"
            unpaired_dataset= "s3d10"
        else:
            paired_dataset= "s3d6"
            unpaired_dataset= "s2d1"

    preprocess_RNA=True
    # imputation=True

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
        ADT_data.obs['batch'] = ADT_data.obs['data_size']




    if preprocess_RNA:
        sc.pp.normalize_total(RNA_data, target_sum=1e4)
        sc.pp.log1p(RNA_data)
    adata = RNA_data
    tmp_ADT = pd.DataFrame(ADT_data.X.toarray(),index = ADT_data.obs.index,columns = ADT_data.var.index)
    adata.obsm["protein_expression"] = tmp_ADT
    batch = adata.obs.batch.values.ravel()
    new_indexs = []
    for indexs,batch1 in zip(list(adata.obs.index),list(adata.obs.batch)):
        if batch1==paired_dataset:
            new_indexs.append(indexs)
        elif batch1==unpaired_dataset:
            new_indexs.append(indexs+'_rna')
        else:
            raise(ValueError('bad batch index'))  
    adata.obs.index = new_indexs 
    adata.obsm["protein_expression"].loc[batch == unpaired_dataset] = np.zeros_like(
        adata.obsm["protein_expression"].loc[batch == unpaired_dataset]
    )
    adata.obsm["protein_expression"].index = adata.obs.index
    sc.pp.highly_variable_genes(
        adata, batch_key="batch", flavor="seurat", n_top_genes=3000, subset=True
    )
    scvi.model.TOTALVI.setup_anndata(
        adata, batch_key="batch", protein_expression_obsm_key="protein_expression"
    )
    model = scvi.model.TOTALVI(adata, latent_distribution="normal", n_layers_decoder=2)
    model.train(use_gpu=gpu_index)
    latent = model.get_latent_representation()
    latent = pd.DataFrame(latent,index = adata.obs.index)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-totalVI-mosaic-latent.csv"))
    if imputation:
        rna, protein = model.get_normalized_expression(
        transform_batch=unpaired_dataset, n_samples=25, return_mean=True
        )
        imputation_file_path = os.path.join(output_path, config["output_prefix"] + "-totalVI-mosaic-imputation.csv")
        imputation_file_path_adt = imputation_file_path.replace('.csv', '_adt.csv')
        imputed_adt = protein
        adt_index = imputed_adt.index.map(lambda x:x.split('_')[-1]=='rna')
        imputed_adt = imputed_adt[adt_index ]
        imputed_adt.to_csv(imputation_file_path_adt)

    # 统计GPU
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()


    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-totalVI-mosaic-gpu_memory.csv'), header=["gpu_memory"])


totalVI_module(input_path = sys.argv[1],
             output_path = sys.argv[2],
             config = json.load(open(sys.argv[3])),
             impute = sys.argv[4]
            )