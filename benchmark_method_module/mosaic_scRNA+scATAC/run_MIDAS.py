#!/home/shaliu_fu/miniconda3/envs/midas_env/bin/python

import warnings
warnings.filterwarnings('ignore')

from scmidas.datasets import GenDataFromPath

from scmidas.models import MIDAS
from scmidas.datasets import GetDataInfo
import scanpy as sc
import pandas as pd
import json
import sys
import os
import numpy as np


from nvitop import Device
import torch
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)
torch.set_num_threads(10) 


def midas_module(input_path,
                 output_path,
                 config,
                 impute=True):
    torch.set_num_threads(10) # set cpu threads limit
    RNA_data = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    ATAC_data = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))
    # mosaic data separation.
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

    paired_rna = RNA_data[RNA_data.obs['batch']==paired_dataset].copy()
    paired_atac = ATAC_data[ATAC_data.obs['batch']==paired_dataset].copy()
    unpaired_rna = RNA_data[RNA_data.obs['batch']==unpaired_dataset].copy()
    unpaired_atac = ATAC_data[ATAC_data.obs['batch']==unpaired_dataset].copy()
    
    atac_out0=pd.DataFrame(paired_atac.X.A,dtype=int)
    atac_out0.index = paired_atac.obs.index.values
    atac_out0.columns = paired_atac.var.index.values

    rna_out0=pd.DataFrame(paired_rna.X.A,dtype=int)
    rna_out0.index = paired_rna.obs.index.values
    rna_out0.columns = paired_rna.var.index.values

    rna_out1=pd.DataFrame(unpaired_rna.X.A,dtype=int)
    rna_out1.index = [f"{i}_rna" for i in unpaired_rna.obs.index.values]
    rna_out1.columns = unpaired_rna.var.index.values

    atac_out2=pd.DataFrame(unpaired_atac.X.A,dtype=int)
    atac_out2.index = [f"{i}_atac" for i in unpaired_atac.obs.index.values]
    atac_out2.columns = unpaired_atac.var.index.values

    # MIDAS data preparation.
    os.system("mkdir -p {}/input/subset_0/mat/".format(output_path))
    os.system("mkdir -p {}/input/subset_1/mat/".format(output_path))
    os.system("mkdir -p {}/input/subset_2/mat/".format(output_path))
    atac_out0.to_csv("{}/input/subset_0/mat/atac.csv".format(output_path))
    rna_out0.to_csv("{}/input/subset_0/mat/rna.csv".format(output_path))
    rna_out1.to_csv("{}/input/subset_1/mat/rna.csv".format(output_path))
    atac_out2.to_csv("{}/input/subset_2/mat/atac.csv".format(output_path))

    data_path = [
        {"rna": "{}/input/subset_0/mat/rna.csv".format(output_path),
        "atac": "{}/input/subset_0/mat/atac.csv".format(output_path)},
        {"rna": "{}/input/subset_1/mat/rna.csv".format(output_path)},
        {"atac": "{}/input/subset_2/mat/atac.csv".format(output_path)},
    ]
    save_dir = "{}/input/".format(output_path)
    remove_old = False
    GenDataFromPath(data_path, save_dir, remove_old) 

    # load and run MIDAS
    data = [GetDataInfo(save_dir)]  # load midas data
    model = MIDAS(data)
    model.init_model()
    model.train(n_epoch=500, save_path="{}/train/".format(output_path)) 
    model.predict(save_dir="{}/predict/".format(output_path), impute=True,mod_latent=True, batch_correct=True,joint_latent=True)

    modality_emb = model.read_preds(impute=True,mod_latent=True, batch_correct=True,joint_latent=True)

    joint_latent = modality_emb["z"]["joint"][:, :32] # biological information
    # joint_batch = modality_emb["z"]["joint"][:, 32:] # batch information
    atac_impute_all = modality_emb['x_impt']['atac']
    rna_impute_all = modality_emb['x_impt']['rna']

    latent_out = pd.DataFrame(joint_latent,index=np.append(rna_out0.index.values,[rna_out1.index.values,atac_out2.index.values]))
    latent_out.to_csv("{}/MIDAS_{}_latent.csv".format(output_path, config['output_prefix']))

    if impute:
        atac_tab = pd.DataFrame(atac_impute_all, index=np.append(rna_out0.index.values,[rna_out1.index.values,atac_out2.index.values]), columns=atac_out0.columns.values)
        atac_tab2 = atac_tab.loc[rna_out1.index.values,:] # atac impute from rna_out1

        rna_tab = pd.DataFrame(rna_impute_all, index=np.append(rna_out0.index.values,[rna_out1.index.values,atac_out2.index.values]), columns=rna_out0.columns.values)
        rna_tab2 = rna_tab.loc[atac_out2.index.values,:]  # rna impute from atac_out2

        atac_tab2.to_csv("{}/MIDAS_{}_imputation_atac.csv.gz".format(output_path, config['output_prefix']), compression='gzip')
        rna_tab2.to_csv("{}/MIDAS_{}_imputation_rna.csv.gz".format(output_path, config['output_prefix']), compression='gzip')

    ## 统计GPU
    pid= os.getpid()
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-MIDAS-mosaic-gpu_memory.csv'), header=["gpu_memory"])


midas_module(input_path = sys.argv[1],
             output_path = sys.argv[2],
             config = json.load(open(sys.argv[3]))
            )