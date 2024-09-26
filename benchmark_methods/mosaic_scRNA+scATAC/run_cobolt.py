import os
import sys
import json
import pickle

import torch
import numpy as np
import pandas as pd

from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt
import anndata as ad


# input_path = sys.argv[1]
# output_path = sys.argv[2]
# config = json.load(open(sys.argv[3]))

# # print(input_path ,output_path,config)



def cobolt_module(input_path,
                 output_path,
                 config):


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
        

    ATAC_data.obs.index = ["{}_atac".format(i) if ATAC_data.obs.loc[i, 'batch'] == unpaired_dataset else "{}".format(i) for i in ATAC_data.obs.index.values]
    RNA_data.obs.index = ["{}_rna".format(i) if RNA_data.obs.loc[i, 'batch'] == unpaired_dataset else "{}".format(i) for i in RNA_data.obs.index.values]


    # Notice that in order to recognize the joint cells, Cobolt requires the barcode and the name of the dataset to be the same. 
    feature = RNA_data.var.index
    barcode = RNA_data.obs.index
    count = RNA_data.X
    rna = SingleData("GeneExpr", "test", feature, count, barcode)


    feature = ATAC_data.var.index
    barcode = ATAC_data.obs.index
    count = ATAC_data.X
    atac = SingleData("ChromAccess", "test", feature, count, barcode)


    rna.filter_features(upper_quantile=0.99, lower_quantile=0.7)
    atac.filter_features(upper_quantile=0.99, lower_quantile=0.7)


    multi_dt = MultiomicDataset.from_singledata(rna, atac)

    from nvitop import Device
    devices = Device.all()  # or Device.cuda.all()
    memory_free = [device.memory_free() for device in devices]
    gpu_index = memory_free.index(max(memory_free))
    print(gpu_index)
    print(memory_free)
    torch.cuda.set_device(gpu_index)
    torch.set_num_threads(5)
    ########################### 

    model = Cobolt(dataset=multi_dt, lr=0.001, n_latent=10)
    model.train(num_epochs=100) 
    model.calc_all_latent()
    latent = model.get_all_latent()


    out_embed = pd.DataFrame(latent[0])
    out_barcode = pd.DataFrame(latent[1])[0]
    out_embed.index = out_barcode
    out_embed.index  = out_embed.index.map(lambda x : x.split('~')[1])
    latent = out_embed

    # latent = latent.round(3)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-cobolt-mosaic-latent.csv"))

    # 统计GPU
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-cobolt-mosaic-gpu_memory.csv'), header=["gpu_memory"])

cobolt_module(input_path = sys.argv[1],
             output_path = sys.argv[2],
             config = json.load(open(sys.argv[3]))
            )