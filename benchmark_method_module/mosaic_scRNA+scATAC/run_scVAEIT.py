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
import json
import h5py

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
def ATAC_process(atac,preprocess=True):
    if preprocess:
        X = atac.X.toarray().astype(np.float32)
        X[X>0.] = 1.
        atac.X = X
    else:
        X = atac.X.toarray().astype(np.float32)
        atac.X = X
    return atac
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





RNA_data = RNA_process(RNA_data,preprocess)
ATAC_data = ATAC_process(ATAC_data) # ATAC没有preprocess
good_peak_header = [f'chr{i}' for i in range(1,23)]
ATAC_data.var['tag'] = ATAC_data.var.index.map(lambda x : x.split('-')[0])
ATAC_data = ATAC_data[:,ATAC_data.var['tag'].isin(good_peak_header)]
ATAC_data.var['tag_i'] = ATAC_data.var['tag'].map(lambda x: int(x.replace("chr", "")))
tmp = ATAC_data.var.copy()
resorted_feature_index = list(tmp.sort_values(by='tag_i').index)
ATAC_data = ATAC_data[:,resorted_feature_index]
gene_names = np.array(RNA_data.var.index).astype(str)
peak_names = np.array(ATAC_data.var.index).astype(str)
all_feature_names = list(RNA_data.var.index) + list(ATAC_data.var.index)
paired_rna = RNA_data[RNA_data.obs['batch']==paired_dataset].copy()
paired_atac = ATAC_data[ATAC_data.obs['batch']==paired_dataset].copy()
unpaired_rna = RNA_data[RNA_data.obs['batch']==unpaired_dataset].copy()
unpaired_atac = ATAC_data[ATAC_data.obs['batch']==unpaired_dataset].copy()
unpaired_atac.obs.index = ["{}_atac".format(i) for i in unpaired_atac.obs.index.values]
unpaired_rna.obs.index = ["{}_rna".format(i) for i in unpaired_rna.obs.index.values]
zero_rna = np.zeros_like(unpaired_rna.X)
zero_atac = np.zeros_like(unpaired_atac.X)
X = np.r_[paired_rna.X,unpaired_rna.X,zero_rna]
Y = np.r_[paired_atac.X,zero_atac,unpaired_atac.X]
cell_barcode = np.array(list(paired_rna.obs.index)+list(unpaired_rna.obs.index)+list(unpaired_atac.obs.index))
len_ref = len(paired_rna.obs.index)
len_gap = len(unpaired_rna.obs.index)
batches = [0]*len_ref+[1]*len_gap+[2]*len_gap
batches = np.array(batches)
batches_tmp = np.zeros((len(batches),2))
batches_tmp[:,1] = batches
batches = batches_tmp
chunk_atac = np.array([
    np.sum(np.char.startswith(peak_names, 'chr%d-'%i)) for i in range(1,23)
    ], dtype=np.int32)
dim_input_arr = np.array([len(gene_names),len(peak_names)])
data = np.c_[X, Y]
masks_raw = - np.ones_like(data, dtype=np.float32)
masks_raw[:len_ref,:]=0. # rna
masks_raw[len_ref:len_ref+len_gap,:len(gene_names)]=0.    # rna有 atac无
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
    # Dimension of input features for [RNA,  peaks]
    'dim_input_arr': dim_input_arr,

    # Blocks for [RNA, ATAC]
    'dim_block': np.append([len(gene_names)], chunk_atac), # input dimension of blocks
    'dist_block':['NB']+ ['Bernoulli' for _ in chunk_atac], # distributions of blocks
    'dim_block_enc':np.array([256]+ [16 for _ in chunk_atac]), # dimension of first layer of the encoder
    'dim_block_dec':np.array([256]+ [16 for _ in chunk_atac]), # dimension of first layer of the decoder
    'dim_block_embed':np.array([16] + [1 for _ in range(len(chunk_atac))])*2, # mask embedding dimension
    'block_names':np.array(['rna'] + ['atac' for _ in range(len(chunk_atac))]),
    'uni_block_names':np.array(['rna','atac']),
    # Internal network structure
    'dimensions':[256], # dimension of latent layers of encoder; the reversed is used for decoder
    'dim_latent':32, # the latent dimension bewteen the encoder and decoder
    'beta_kl':1.,
    
    # Weights
    'beta_unobs':.9, # weight for masked out observation; weight for observerd values will be 1-beta_unobs.
    'beta_modal':np.array([0.95,0.05]), # weights for 2 modalities, which can be adjusted based on loss in the first few epochs.
    'beta_reverse':0.,
    # Masking probability
    "p_feat" : 0.5, # probablity of randomly masking out an entry
    "p_modal" : np.ones(2)/2,
}

batches_cate = batches


model = scVAEIT(configs, data, masks, batches_cate)
# 只以5个epoch为例子，真正运行的话300个epoch要跑三四个小时
# 内存吃的很多
epoch = 300
model.train(
        num_epoch=epoch, batch_size=256, save_every_epoch=50,
        verbose=True, checkpoint_dir=path_root+'/checkpoint/')
masks = masks.numpy().copy() 
if imputation:
    imputation_file_path = os.path.join(output_path, config["output_prefix"] + "-scVAEIT-mosaic-imputation.csv")
    imputation_file_path_rna = imputation_file_path.replace('.csv', '_rna.csv')
    imputation_file_path_atac = imputation_file_path.replace('.csv', '_atac.csv')
    denoised_data = model.get_denoised_data() # 这个是全部的数据，它在inputation的时候顺便还改变了其它的值
    denoised_data = pd.DataFrame(denoised_data,index=cell_barcode,columns=all_feature_names)
    imputed_expression = denoised_data.iloc[:,:len(gene_names)]
    rna_index = imputed_expression.index.map(lambda x:x.split('_')[-1]=='atac')
    imputed_rna = imputed_expression[rna_index ]
    #imputed_rna = imputed_rna.round(3)
    imputed_rna.to_csv(imputation_file_path_rna)
    imputed_atac = denoised_data.iloc[:,len(gene_names):]
    atac_index = imputed_atac.index.map(lambda x:x.split('_')[-1]=='rna')
    imputed_atac = imputed_atac[atac_index ]
    #imputed_atac = imputed_atac.round(3)
    imputed_atac.to_csv(imputation_file_path_atac)
model.update_z()
latent_df = pd.DataFrame(model.adata.X)
latent_df.index = list(cell_barcode)
latent = latent_df
#latent = latent.round(3)
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
        
