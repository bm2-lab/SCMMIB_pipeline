

import sys
sys.path.append('/home/sirm/project/software/scVAEIT-main/Reproducibility materials/')
from scVAEIT.VAEIT import scVAEIT

import anndata as ad

import os
import sys
os.environ["OMP_NUM_THREADS"] = "11"
os.environ["OPENBLAS_NUM_THREADS"] = "8" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "11" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "8" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "11" # export NUMEXPR_NUM_THREADS=6
os.environ["NUMBA_CACHE_DIR"]='/tmp/numba_cache'
import numpy as np
import pandas as pd
import scipy as sp
import scipy.sparse
import h5py
import json
import tensorflow as tf
import matplotlib.pyplot as plt
import scanpy as sc


def RNA_process(rna,preprocess=True):
    if preprocess:
        X = rna.X.toarray().astype(np.float32)
        rna.X = np.log(X/np.sum(X, axis=1, keepdims=True)*1e4+1.) 
    else:
        X = rna.X.toarray().astype(np.float32)
        rna.X = X
    return rna
def ADT_process(adt,preprocess=True):
    if preprocess:
        X = adt.X.toarray().astype(np.float32)
        adt.X = np.log(X/np.sum(X, axis=1, keepdims=True)*1e4+1.) 
    else:
        X = adt.X.toarray().astype(np.float32)
        adt.X = X
    return adt

##### 这里替换原来的输入
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# print(input_path ,output_path,config)




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











RNA_data=RNA_process(RNA_data,preprocess_RNA)
ADT_data=ADT_process(ADT_data,preprocess_ADT)

gene_names = np.array(RNA_data.var.index).astype(str)
ADT_names = np.array(ADT_data.var.index).astype(str)

all_feature_names = list(RNA_data.var.index) + list(ADT_data.var.index)

paired_rna = RNA_data[RNA_data.obs['batch']==paired_dataset].copy()
paired_ADT = ADT_data[ADT_data.obs['batch']==paired_dataset].copy()
unpaired_rna = RNA_data[RNA_data.obs['batch']==unpaired_dataset].copy()
unpaired_ADT = ADT_data[ADT_data.obs['batch']==unpaired_dataset].copy()

unpaired_ADT.obs.index = ["{}_adt".format(i) for i in unpaired_ADT.obs.index.values]
unpaired_rna.obs.index = ["{}_rna".format(i) for i in unpaired_rna.obs.index.values]

zero_rna = np.zeros_like(unpaired_rna.X)
zero_adt = np.zeros_like(unpaired_ADT.X)

X = np.r_[paired_rna.X,unpaired_rna.X,zero_rna]
Y = np.r_[paired_ADT.X,zero_adt,unpaired_ADT.X]

cell_barcode = np.array(list(paired_rna.obs.index)+list(unpaired_rna.obs.index)+list(unpaired_ADT.obs.index))

len_ref = len(paired_rna.obs.index)
len_gap = len(unpaired_rna.obs.index)


batches = [0]*len_ref+[1]*len_gap+[2]*len_gap
batches = np.array(batches)
batches_tmp = np.zeros((len(batches),2))
batches_tmp[:,1] = batches
batches = batches_tmp

dim_input_arr = np.array([len(gene_names),len(ADT_names)])
data = np.c_[X, Y]
masks_raw = - np.ones_like(data, dtype=np.float32)

masks_raw[:len_ref,:]=0. # rna
masks_raw[len_ref:len_ref+len_gap,:len(gene_names)]=0.    # rna有 adt无
masks_raw[len_ref+len_gap:,len(gene_names):]=0 

from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))


physical_devices = tf.config.list_physical_devices('GPU')
try:
    tf.config.experimental.set_visible_devices(physical_devices[gpu_index], 'GPU')
    #tf.config.experimental.set_memory_growth(physical_devices[0], True)
except:
    # Invalid device or cannot modify virtual devices once initialized.
    pass

masks = masks_raw.copy()
masks = tf.convert_to_tensor(masks, dtype=tf.float32)
tmp_checkpoint_path = output_path



path_root = tmp_checkpoint_path
configs = {
    # Dimension of input features for [RNA, ADT, peaks]
    'dim_input_arr': dim_input_arr,

    # Blocks for [RNA, ADT]
    'dim_block': np.array([len(gene_names),len(ADT_names)]), # input dimension of blocks
    'dist_block':['NB','NB'], # distributions of blocks
    'dim_block_enc':np.array([128, 64]), # dimension of first layer of the encoder
    'dim_block_dec':np.array([128, 64]), # dimension of first layer of the decoder
    'dim_block_embed':np.array([128, 64]), # mask embedding dimension
    'block_names':np.array(['rna', 'adt'] ),
    'uni_block_names':np.array(['rna','adt']),
    # Internal network structure
    'dimensions':[256], # dimension of latent layers of encoder; the reversed is used for decoder
    'dim_latent':32, # the latent dimension bewteen the encoder and decoder
    'beta_kl':1.,
    
    # Weights
    'beta_unobs':.9, # weight for masked out observation; weight for observerd values will be 1-beta_unobs.
    'beta_modal':np.array([0.05,0.95]), # weights for 2 modalities, which can be adjusted based on loss in the first few epochs.
    'beta_reverse':0.,
    # Masking probability
    "p_feat" : 0.5, # probablity of randomly masking out an entry
    "p_modal" : np.ones(2)/2,
}


batches_cate = batches


model = scVAEIT(configs, data, masks, batches_cate)

epoch = 300
model.train(
        num_epoch=epoch, batch_size=256, save_every_epoch=50,
        verbose=True, checkpoint_dir=path_root+'/checkpoint/')

masks = masks.numpy().copy() 



if imputation:
    imputation_file_path = os.path.join(output_path, config["output_prefix"] + "-scVAEIT-mosaic-imputation.csv")
    imputation_file_path_rna = imputation_file_path.replace('.csv', '_rna.csv')
    imputation_file_path_adt = imputation_file_path.replace('.csv', '_adt.csv')
    denoised_data = model.get_denoised_data() # 这个是全部的数据，它在inputation的时候顺便还改变了其它的值
    denoised_data = pd.DataFrame(denoised_data,index=cell_barcode,columns=all_feature_names)
    imputed_expression = denoised_data.iloc[:,:len(gene_names)]
    rna_index = imputed_expression.index.map(lambda x:x.split('_')[-1]=='adt')
    imputed_rna = imputed_expression[rna_index ]
    # imputed_rna = imputed_rna.round(3)
    imputed_rna.to_csv(imputation_file_path_rna)
    imputed_adt = denoised_data.iloc[:,len(gene_names):]
    adt_index = imputed_adt.index.map(lambda x:x.split('_')[-1]=='rna')
    imputed_adt = imputed_adt[adt_index ]
    # imputed_adt = imputed_adt.round(3)
    imputed_adt.to_csv(imputation_file_path_adt)
    
model.update_z()
latent_df = pd.DataFrame(model.adata.X)

latent_df.index = list(cell_barcode)
latent = latent_df
# latent = latent.round(3)
latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-scVAEIT-mosaic-latent.csv"))
## 统计GPU
pid= os.getpid()        
gpu_memory = pd.Series(dtype='str')

devices = Device.all()
for device in devices:
    processes = device.processes()    
    if pid in processes.keys():
        p=processes[pid]
        gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()


gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-scVAEIT-mosaic-gpu_memory.csv'), header=["gpu_memory"])