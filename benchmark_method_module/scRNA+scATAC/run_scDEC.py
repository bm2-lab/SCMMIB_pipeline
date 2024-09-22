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
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import io

# https://github.com/kimmo1019/scDEC

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_scDEC"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# %%
# Make Dir
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Load data
rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

# Load Metadata
metadata = pd.read_csv(os.path.join(input_path, config['metadata']), header=0)
metadata.index=metadata[config['barcode_key']].values

# %%
rna.var['feature_types'] = 'Gene Expression'
atac.var['feature_types'] = 'Peaks'
multi = ad.concat([rna, atac], axis=1, join="outer")

# %%
# # cp the scDEC source code to output_path
scMDC_path = os.path.join(config['utils_path'], 'scDEC')
os.system("cp -r " + scMDC_path + " " + output_path + "/scDEC_src")

# %%
# Save data
data_path = os.path.join(output_path, "scDEC_src/datasets/data")
os.makedirs(data_path)
## save X to mtx
io.mmwrite(os.path.join(data_path, "matrix.mtx"), multi.X.T)
## save barcodes
with open(os.path.join(data_path, "barcodes.tsv"), 'w') as f:
    for item in multi.obs_names:
        f.write(item + '\n')      
## save features
multi.var.to_csv(os.path.join(data_path, "features.tsv"), sep="\t", header=False)

# %%
os.chdir(os.path.join(output_path, "scDEC_src"))

# %%
conda_env = '/home/wsg/software/miniconda3/envs/scDEC/bin/python'

# %%
bash_cmd1 = conda_env + ' main_clustering.py --data data --K 4 --dx 16 --dy 50 --alpha 10 --beta 10 --ratio 0.2 --train True --mode 3 --no_label'
bash_cmd2 = conda_env + ' eval.py --data data --timestamp joint_bench --train True --no_label --epoch 30000'
with open("run_scDEC.sh", "w") as sh_file:
    sh_file.write(bash_cmd1 + '\n')
    sh_file.write(bash_cmd2 + '\n')

# %%
os.system("bash run_scDEC.sh")

# %%
os.system("mv results/data/*/scDEC* " + output_path)

# %%
os.chdir(output_path)

# %%
## save latent
latent = pd.read_csv(os.path.join(output_path, "scDEC_embedding.csv"), index_col=0, header=0)
latent.index = rna.obs_names
latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-scDEC-multi-latent.csv"))

# %%
## save umap
umap = pd.read_csv(os.path.join(output_path, "scDEC_UMAP.csv"), index_col=0, header=0)
umap.index = rna.obs_names
umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-scDEC-multi-umap.csv"))

# %%
os.remove("scDEC_UMAP.csv")
os.remove("scDEC_cluster.txt")
os.remove("scDEC_embedding.csv")
os.system("rm -r scDEC_src")

# %%
