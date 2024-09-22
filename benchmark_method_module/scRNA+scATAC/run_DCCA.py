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
#     display_name: python (DCCA)
#     language: python
#     name: dcca
# ---

# %%
import os
import sys 
import json
import glob

import numpy as np
import pandas as pd
import scanpy as sc
import pyranges as pr
from scanpy import AnnData

import torch
import torch.utils.data as data_utils
from sklearn.metrics import precision_recall_curve, auc

from DCCA.MVAE_cycleVAE import DCCA
from DCCA.utilities import read_dataset, normalize, parameter_setting, save_checkpoint, load_checkpoint

# https://github.com/cmzuo11/DCCA/wiki/Analysis-of-PBMC_3K-dataset-from-10X-Genomics-by-DCCA-model

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def train_with_argas( args ):
    torch.set_num_threads(5)
    
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    model_file      =  os.path.join( args.outdir, 'model_DCCA.pth.tar' )

    # Load data
    adata = sc.read_h5ad(os.path.join(input_path, args.config['rna_h5ad_filename']))
    if not isinstance(adata.X, np.ndarray):
        adata.X = adata.X.toarray().astype('float32')

    adata1 = sc.read_h5ad(os.path.join(input_path, args.config['atac_h5ad_filename']))
    if not isinstance(adata1.X, np.ndarray):
        adata1.X = adata1.X.toarray().astype('float32')
    
    # Feature Selection: RNA
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
    hvg = adata.var[adata.var.highly_variable]
    
    # Feature Selection: ATAC
    gtf = pr.read_gtf(config['gtf_file'])
    tmp_bed = os.path.join(output_path, "tmp.bed")
    if config["specie"] == "human":
        chr_normal = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    elif config["specie"] == "mouse":
        chr_normal = ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY"]
    gtf_gene_index = (gtf.Feature == "gene") & (gtf.Chromosome.isin(chr_normal))  & (gtf.gene_name.isin(hvg.index.values))
    gtf_gene = gtf[gtf_gene_index]

    peaks = adata1.var.index.values
    peaks_split = [i.split("-") for i in peaks]
    with open(tmp_bed,"w") as fo:
        for line in peaks_split:
            fo.write("{}\n".format("\t".join(line)))

    peaks_bed = pr.read_bed(tmp_bed)
    gtf_gene2 = gtf_gene.extend({"5":100000}) # 5‘ 100kb
    peaks_select = peaks_bed.intersect(gtf_gene2)
    peak_name=peaks_select.to_csv(sep="-")
    peak_name2=peak_name.split("\n")[1:]
    
    adata = adata[:, adata.var.features.isin(hvg.index.values)]
    adata1 = adata1[:, adata1.var.features.isin(peak_name2)]

    adata  = normalize(adata, filter_min_counts=False, size_factors=True, 
                       normalize_input=False, logtrans_input=True ) 
    adata1 = normalize(adata1, filter_min_counts=False, size_factors=False, 
                       normalize_input=False, logtrans_input=False ) 

    train_index, test_index = list(range( adata.n_obs )), list(range( adata.n_obs ))

    Nsample1, Nfeature1   =  np.shape( adata.X )
    Nsample2, Nfeature2   =  np.shape( adata1.X )

    train           = data_utils.TensorDataset( torch.from_numpy( adata[train_index].X ),
                                                torch.from_numpy( adata.raw[train_index].X ), 
                                                torch.from_numpy( adata.obs['size_factors'][train_index].values ),
                                                torch.from_numpy( adata1[train_index].X ),
                                                torch.from_numpy( adata1.raw[train_index].X ), 
                                                torch.from_numpy( adata1.obs['size_factors'][train_index].values ))

    train_loader    = data_utils.DataLoader( train, batch_size = args.batch_size, shuffle = True )

    test            = data_utils.TensorDataset( torch.from_numpy( adata[test_index].X ),
                                                torch.from_numpy( adata.raw[test_index].X ), 
                                                torch.from_numpy( adata.obs['size_factors'][test_index].values ),
                                                torch.from_numpy( adata1[test_index].X ),
                                                torch.from_numpy( adata1.raw[test_index].X ), 
                                                torch.from_numpy( adata1.obs['size_factors'][test_index].values ) )

    test_loader     = data_utils.DataLoader( test, batch_size = len(test_index), shuffle = False )

    total           = data_utils.TensorDataset( torch.from_numpy( adata.X  ),
                                                torch.from_numpy( adata.raw.X ),
                                                torch.from_numpy( adata.obs['size_factors'].values ),
                                                torch.from_numpy( adata1.X  ),
                                                torch.from_numpy( adata1.raw.X ),
                                                torch.from_numpy( adata1.obs['size_factors'].values ) )
    total_loader    = data_utils.DataLoader( total, batch_size = (len(train_index)+ len(test_index)) , shuffle = False )

    label_ground_truth=np.ones(len(adata.obs_names))

    model =  DCCA( layer_e_1 = [Nfeature1, 128], hidden1_1 = 128, Zdim_1 = 4, layer_d_1 = [4, 128],
                   hidden2_1 = 128, layer_e_2 = [Nfeature2, 1500, 128], hidden1_2 = 128, Zdim_2 = 4,
                   layer_d_2 = [4], hidden2_2 = 4, args = args, ground_truth = label_ground_truth,
                   ground_truth1 = label_ground_truth, Type_1 = "NB", Type_2 = "Bernoulli", cycle = 1, 
                   attention_loss = "Eucli" )

    if args.use_cuda:
        model.cuda()

    NMI_score1, ARI_score1, NMI_score2, ARI_score2  =  model.fit_model(train_loader, test_loader, total_loader, "RNA" )

    save_checkpoint(model, model_file ) 

    latent_z1, latent_z2, norm_x1, _, norm_x2, _ = model( total_loader )

    if latent_z1 is not None:
        pd.DataFrame( latent_z1, index= adata.obs_names ).to_csv( os.path.join( args.outdir, args.config["output_prefix"] + "-DCCA-RNA-latent_1.csv" ) ) 

    if latent_z2 is not None:
        pd.DataFrame( latent_z2, index= adata1.obs_names ).to_csv( os.path.join( args.outdir, args.config["output_prefix"] + "-DCCA-ATAC-latent_2.csv" ) ) 

    latent = np.hstack((latent_z1, latent_z2))
    pd.DataFrame( latent, index= adata.obs_names ).to_csv( os.path.join( args.outdir, args.config["output_prefix"] + "-DCCA-multi-latent.csv" ) ) 
    
    os.remove(os.path.join(output_path, "model_DCCA.pth.tar"))
    
    # imputation
    sys.path.append(args.config['utils_path'])
    from knn_smooth import knn_smoothing

    # metadata
    metadata = pd.read_csv(os.path.join(input_path, args.config['metadata']))
    metadata.index = metadata[args.config['barcode_key']]

    
    # imputation
    ## RNA
    imputation_rna = pd.DataFrame(norm_x1)
    imputation_rna.columns = adata.var.index.values
    imputation_rna.index = adata.obs.index.values
    imputation_rna.to_csv(os.path.join(output_path, config["output_prefix"] + "-DCCA-imputation-rna.csv.gz"), compression='gzip')

    ## ATAC
    imputation_atac = pd.DataFrame(norm_x2)
    imputation_atac.columns = adata1.var.index.values
    imputation_atac.index = adata1.obs.index.values
    imputation_atac.to_csv(os.path.join(output_path, config["output_prefix"] + "-DCCA-imputation-atac.csv.gz"), compression='gzip')

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    if len(gpu_memory):
        gpu_memory.to_csv(os.path.join(output_path, args.config['output_prefix'] + '-DCCA-gpu_memory.csv'), header=["gpu_memory"])



# %%

# %%
input_path = sys.argv[1]
output_path = sys.argv[2]
config = json.load(open(sys.argv[3]))

# %%
parser = parameter_setting()
args   = parser.parse_args("")

args.workdir  =  input_path
args.outdir   =  output_path
args.config  =  config

args.batch_size     = 64
args.epoch_per_test = 10
args.use_cuda       =  args.use_cuda and torch.cuda.is_available()

args.sf1        =  5
args.sf2        =  1
args.cluster1   =  args.cluster2   =  4
args.lr1        =  0.01
args.flr1       =  0.001
args.lr2        =  0.005
args.flr2       =  0.0005

# %%
train_with_argas(args)

# %%

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/data/test/run_DCCA"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
# input_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/p10"
# output_path = "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/BMMC/p10/rep_2/run_DCCA"
# config = json.load(open("/home/wsg/BM/data/BMMC/RNA+ATAC/p10/p10.json"))

# %%
# parser = parameter_setting()
# args   = parser.parse_args("")

# args.workdir  =  input_path
# args.outdir   =  output_path
# args.config  =  config

# args.batch_size     = 64
# args.epoch_per_test = 10
# args.use_cuda       =  args.use_cuda and torch.cuda.is_available()

# args.sf1        =  5
# args.sf2        =  1
# args.cluster1   =  args.cluster2   =  4
# args.lr1        =  0.01
# args.flr1       =  0.001
# args.lr2        =  0.005
# args.flr2       =  0.0005

# %%
# train_with_argas(args)

# %%

# %%

# %%
