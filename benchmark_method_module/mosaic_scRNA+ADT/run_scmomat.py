
import sys,os,json
import numpy as np
# from umap import UMAP
import time
import torch
import matplotlib.pyplot as plt
import pandas as pd  
import scipy.sparse as sp
import anndata as ad
import scmomat 

plt.rcParams["font.size"] = 10

##### 这里替换原来的输入
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# print(input_path ,output_path,config)




RNA_data = ad.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
ADT_data = ad.read_h5ad(os.path.join(input_path, config['adt_h5ad_filename']))


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
    ADT_data.obs['batch'] = ADT_data.obs['data_size']

####


genes = np.array(RNA_data.var.index)
proteins = np.array(ADT_data.var.index)


paired_rna = RNA_data[RNA_data.obs['batch']==paired_dataset].copy()
paired_ADT = ADT_data[ADT_data.obs['batch']==paired_dataset].copy()
unpaired_rna = RNA_data[RNA_data.obs['batch']==unpaired_dataset].copy()
unpaired_ADT = ADT_data[ADT_data.obs['batch']==unpaired_dataset].copy()


unpaired_ADT.obs.index = ["{}_adt".format(i) for i in unpaired_ADT.obs.index.values]
unpaired_rna.obs.index = ["{}_rna".format(i) for i in unpaired_rna.obs.index.values]


# 得到cell barcode
cell_barcode_ref = paired_rna.obs.index
cell_barcode_rna = unpaired_rna.obs.index
cell_barcode_adt = unpaired_ADT.obs.index
# 得到array
paired_rna = paired_rna.X.toarray()
paired_ADT = paired_ADT.X.toarray()
unpaired_rna = unpaired_rna.X.toarray()
unpaired_ADT = unpaired_ADT.X.toarray()



if preprocess_RNA and preprocess_ADT:
    paired_rna = scmomat.preprocess(paired_rna, modality = "RNA", log = False)
    unpaired_rna = scmomat.preprocess(unpaired_rna, modality = "RNA", log = False)
    paired_ADT = scmomat.preprocess(paired_ADT, modality = "protein")
    unpaired_ADT = scmomat.preprocess(unpaired_ADT, modality = "protein")
elif preprocess_RNA and not preprocess_ADT:
    paired_rna = scmomat.preprocess(paired_rna, modality = "RNA", log = False)
    unpaired_rna = scmomat.preprocess(unpaired_rna, modality = "RNA", log = False)
elif not preprocess_RNA and preprocess_ADT:
    paired_ADT = scmomat.preprocess(paired_ADT, modality = "protein")
    unpaired_ADT = scmomat.preprocess(unpaired_ADT, modality = "protein")
else:
    pass


# batch是按照ref gap的顺序排的，因此后续对进去cell barcode的时候也要这样
counts_rnas = [paired_rna, unpaired_rna, None]
counts_adts = [paired_ADT,None, unpaired_ADT]


feats_name = {"rna": genes, "adt": proteins}
# CREATE THE COUNTS OBJECT
counts = {"feats_name": feats_name, "nbatches": 3, "rna":counts_rnas, "adt": counts_adts}

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
# running device, can be CPU or GPU

from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
print(gpu_index)
print(memory_free)
torch.cuda.set_device(gpu_index)
torch.set_num_threads(5)




device = torch.device("cuda:{}".format(gpu_index) if torch.cuda.is_available() else "cpu")
#------------------------------------------------------------------------------------------------------------------------------------

model = scmomat.scmomat_model(counts = counts, K = K, batch_size = batch_size, interval = interval, lr = lr, lamb = lamb, seed = seed,device=device)
losses = model.train_func(T = T)


zs = model.extract_cell_factors()

# zs是批次划分的latent，这里为了输出需要合并
# 模型不输入barcode，直接输出就行。
latents = np.concatenate(zs, axis = 0)
latent_df = pd.DataFrame(latents)
latent_df.index = list(cell_barcode_ref)+list(cell_barcode_rna)+list(cell_barcode_adt)
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