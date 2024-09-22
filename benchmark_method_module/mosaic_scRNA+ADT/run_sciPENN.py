import os
import json
import sys
import numpy as np
import scanpy as sc
import pandas as pd
import torch
from sciPENN.sciPENN_API import sciPENN_API


from nvitop import Device
devices = Device.all()  # or Device.cuda.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
print(gpu_index)
print(memory_free)
torch.cuda.set_device(gpu_index)
torch.set_num_threads(5)


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


    
paired_rna = RNA_data[RNA_data.obs['batch']==paired_dataset].copy()
paired_ADT = ADT_data[ADT_data.obs['batch']==paired_dataset].copy()
unpaired_rna = RNA_data[RNA_data.obs['batch']==unpaired_dataset].copy()
unpaired_ADT = ADT_data[ADT_data.obs['batch']==unpaired_dataset].copy()
unpaired_rna.obs.index = ["{}_rna".format(i) for i in unpaired_rna.obs.index.values]
unpaired_ADT.obs.index = ["{}_adt".format(i) for i in unpaired_ADT.obs.index.values]


if preprocess_RNA:
    sciPENN = sciPENN_API(gene_trainsets = [paired_rna],
                      protein_trainsets = [paired_ADT],
                      gene_test=unpaired_rna,
                      train_batchkeys=['batch'],
                      test_batchkey='batch',
                    #   type_key='cell_type'
                      )
else:
    sciPENN = sciPENN_API(gene_trainsets = [paired_rna],
                      protein_trainsets = [paired_ADT],
                      gene_test=unpaired_rna,
                      train_batchkeys=['batch'],
                      test_batchkey='batch',
                      cell_normalize = False, 
                      log_normalize = False, 
                      gene_normalize = False,
                    #   type_key='cell_type'
                      )


tmp_checkpoint_path = os.path.dirname(os.path.join(output_path, config["output_prefix"] + "-sciPENN-mosaic-latent.csv"))
path_root = tmp_checkpoint_path

print(path_root)

sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, 
             decay_step = 0.1, lr = 10**(-3), load = False)
# sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, 
#              decay_step = 0.1, lr = 10**(-3), weights_dir = './checkpoint/', load = False)


embedding = sciPENN.embed()
latent = pd.DataFrame(embedding.X,index=embedding.obs.index)
latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-sciPENN-mosaic-latent.csv"))

if imputation:
    predicted_test = sciPENN.predict()
    imputation_file_path = os.path.join(output_path, config["output_prefix"] + "-sciPENN-mosaic-imputation.csv")
    imputation_file_path_adt = imputation_file_path.replace('.csv', '_adt.csv')
    imputed_adt = pd.DataFrame(predicted_test.X,index=predicted_test.obs.index,columns=paired_ADT.var.index)
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


gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-sciPENN-mosaic-gpu_memory.csv'), header=["gpu_memory"])


