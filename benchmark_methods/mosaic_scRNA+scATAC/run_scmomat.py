
import sys, os,json

import numpy as np
# from umap import UMAP
import time
import torch
import matplotlib.pyplot as plt
import pandas as pd  
import scipy.sparse as sp

import scmomat 
import anndata as ad

plt.rcParams["font.size"] = 10


def scmomat_module(input_path,
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
            
    preprocess=True

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

    genes = np.array(RNA_data.var.index)
    regions = np.array(ATAC_data.var.index)
    feats_name = {"rna": genes, "atac": regions}
    count_rna_ref = RNA_data[RNA_data.obs['batch']==paired_dataset].X.toarray()
    count_rna_gap = RNA_data[RNA_data.obs['batch']==unpaired_dataset].X.toarray()
    count_atac_ref = ATAC_data[ATAC_data.obs['batch']==paired_dataset].X.toarray()
    count_atac_gap = ATAC_data[ATAC_data.obs['batch']==unpaired_dataset].X.toarray()
    cell_barcode_ref = RNA_data[RNA_data.obs['batch']==paired_dataset].obs.index
    cell_barcode_gap = RNA_data[RNA_data.obs['batch']==unpaired_dataset].obs.index
    cell_barcode_rna = cell_barcode_gap.map(lambda x : x+'_rna')
    cell_barcode_atac = cell_barcode_gap.map(lambda x : x+'_atac')


    if preprocess:
        count_rna_ref = scmomat.preprocess(count_rna_ref, modality = "RNA", log = False)
        count_rna_gap = scmomat.preprocess(count_rna_gap, modality = "RNA", log = False)
        count_atac_ref = scmomat.preprocess(count_atac_ref, modality = "ATAC")
        count_atac_gap = scmomat.preprocess(count_atac_gap, modality = "ATAC")
    else:
        count_atac_ref = scmomat.preprocess(count_atac_ref, modality = "ATAC")
        count_atac_gap = scmomat.preprocess(count_atac_gap, modality = "ATAC")

    counts_atacs = [count_atac_ref, None,count_atac_gap]
    counts_rnas = [count_rna_ref, count_rna_gap, None]

    # CREATE THE COUNTS OBJECT
    counts = {"feats_name": feats_name, "nbatches": 3, "rna":counts_rnas, "atac": counts_atacs}
    #------------------------------------------------------------------------------------------------------------------------------------
    # NOTE: Number of latent dimensions, key hyper-parameter, 20~30 works for most of the cases.
    K = 30
    #------------------------------------------------------------------------------------------------------------------------------------
    # NOTE: Here we list other parameters in the function for illustration purpose, most of these parameters are set as default value.
    # weight on regularization term, default value
    lamb = 0.001 
    # number of total iterations, default value
    T = 4000
    # print the result after each ``interval'' iterations, default value
    interval = 1000
    # batch size for each iteraction, default value
    batch_size = 0.1
    # learning rate, default value
    lr = 1e-2
    # random seed, default value
    seed = 0
    from nvitop import Device
    devices = Device.all()  # or Device.cuda.all()
    memory_free = [device.memory_free() for device in devices]
    gpu_index = memory_free.index(max(memory_free))
    print(gpu_index)
    print(memory_free)
    torch.cuda.set_device(gpu_index)
    torch.set_num_threads(5)

    # running device, can be CPU or GPU
    device = torch.device("cuda:{}".format(gpu_index) if torch.cuda.is_available() else "cpu")

    #------------------------------------------------------------------------------------------------------------------------------------

    model = scmomat.scmomat_model(counts = counts, K = K, batch_size = batch_size, interval = interval, lr = lr, lamb = lamb, seed = seed, device = device)
    losses = model.train_func(T = T)

    zs = model.extract_cell_factors()

    latents = np.concatenate(zs, axis = 0)
    latent_df = pd.DataFrame(latents)
    latent_df.index = list(cell_barcode_ref)+list(cell_barcode_rna)+list(cell_barcode_atac)
    latent = latent_df

    #latent = latent.round(3)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-scmomat-mosaic-latent.csv"))

    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()


    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-scmomat-mosaic-gpu_memory.csv'), header=["gpu_memory"])



scmomat_module(input_path = sys.argv[1],
             output_path = sys.argv[2],
             config = json.load(open(sys.argv[3]))
            )