# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: python (benchmark)
#     language: python
#     name: benchmark
# ---

# %%
import os
import sys 
import json

import h5py
import torch
import numpy as np
import pandas as pd
import scanpy as sc

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)

# %%
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# %%
# input_path = "/home/wsg/BM/data/COVID19/RNA+ADT/c1k"
# output_path = "/home/wsg/BM/results/task/scRNA+ADT/accuracy/COVID19/RawData/rep_1/run_scMDC"
# config = json.load(open("/home/wsg/BM/data/COVID19/RNA+ADT/c1k/c1k.json"))

# %%
# Make Dir
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Load data
rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
adt = sc.read_h5ad(os.path.join(input_path, config['adt_h5ad_filename']))
rna.var_names_make_unique()
adt.var_names_make_unique()

# Load Metadata
metadata = pd.read_csv(os.path.join(input_path, config['metadata']), header=0)
metadata.index=metadata[config['barcode_key']].values

# %%
os.chdir(output_path)
h5_file = os.path.join(output_path, config['output_prefix'] + "-scMDC-merged-matrix.h5")
h5f = h5py.File(h5_file, 'w')
h5f.create_dataset('X1', data=rna.X.todense(), dtype='f8', compression='gzip')
h5f.create_dataset('X2', data=adt.X.todense(), dtype='f8', compression='gzip')
h5f.create_dataset('Y', data=np.zeros((rna.X.shape[0])), compression='gzip')
h5f.close()


# %%
# # cp the scMDC source code to output_path
scMDC_path = os.path.join(config['utils_path'], 'scMDC')
os.system("cp -r " + scMDC_path + " " + output_path + "/scMDC_src")

# %%
conda_env = '/home/wsg/software/miniconda3/envs/benchmark/bin/python'

# %%
bash_cmd = conda_env + ' scMDC_src/src/run_scMDC.py --n_clusters 15 --data_file ' + h5_file + ' --ae_weight_file AE_weights.pth.tar \
--save_dir ' + output_path + ' --embedding_file --filter1 --filter2 --f1 2000 --f2 2000 \
-el 256 128 64 -dl1 64 128 256 -dl2 64 128 256 --phi1 0.005 --phi2 0.005 --sigma2 2.5 --tau .1'
with open("run_scMDC.sh", "w") as sh_file:
    sh_file.write(bash_cmd)

# %%
os.system("bash run_scMDC.sh")

# %%
latent = pd.read_csv('1_embedding.csv', header=None) 
latent.columns = pre_res = ['feat_' + str(colname) for colname in latent.columns]
latent.index = rna.obs_names

latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMDC-multi-latent.csv"))

# %%
os.remove(h5_file)
os.remove("1_embedding.csv")
os.remove("AE_weights.pth.tar")
os.system("rm -r scMDC_src")

# %%
