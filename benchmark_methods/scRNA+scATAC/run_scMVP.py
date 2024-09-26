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

import muon
import torch
import anndata
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
import matplotlib.pyplot as plt

from scipy import io
import scipy.io as sp_io
from scipy.sparse import csr_matrix, issparse
from sklearn.metrics import precision_recall_curve, auc

from scMVP.dataset import LoadData, GeneExpressionDataset, CellMeasurement
from scMVP.models import VAE_Attention, Multi_VAE_Attention, VAE_Peak_SelfAttention
from scMVP.inference import UnsupervisedTrainer
from scMVP.inference import MultiPosterior, MultiTrainer


# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def scmvp_module(input_path,
                 output_path,
                 config):
    torch.set_num_threads(5)

    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    if rna.X.shape[0] <= 5000:
        train_epochs = [10, 5, 10]
    elif rna.X.shape[0] > 5000 and rna.X.shape[0] <= 20000:
        train_epochs = [15, 10, 15]
    elif rna.X.shape[0] > 20000:
        train_epochs = [15, 15, 15]

    print(train_epochs[0])
    print(train_epochs[1])
    print(train_epochs[2])

    # ATAC preprocess
    atac.layers['raw'] = atac.X.copy()
    muon.atac.pp.tfidf(atac, scale_factor=1e4)  # run tf-idf
    atac.layers['tf-idf'] = atac.X.copy()
    atac.X = atac.layers['raw'].copy()
    sc.pp.normalize_total(atac, target_sum=1e4)  # run log-norm
    sc.pp.log1p(atac)

    ## Feature selection
    # hvp=min(100000,atac.X.shape[1])
    sc.pp.highly_variable_genes(atac, n_top_genes=30000)
    atac.X = atac.layers['tf-idf'].copy()
    atac = atac[:, atac.var.highly_variable].copy()

    ## make ATAC object 
    scmvp_atac_dataset = GeneExpressionDataset()
    cell_attributes_dict = {"barcodes": atac.obs_names.to_numpy()}
    scmvp_atac_dataset.populate_from_data(
        X=atac.X.toarray(),  # notice the normalization
        batch_indices=None,
        gene_names=atac.var_names.to_numpy(),
        cell_attributes_dict=cell_attributes_dict,
        Ys=[]
    )

    # RNA preprocess
    rna.layers['raw'] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)  # run log-norm
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000)

    # Feature selection
    rna.X = rna.layers['raw'].copy()
    rna_hvg = rna[:, rna.var.highly_variable].copy()

    ## make Y for RNA object
    scmvp_rna_dataset = GeneExpressionDataset()
    Ys = []
    measurement = CellMeasurement(
        name="atac_expression",
        data=scmvp_atac_dataset.X,
        columns_attr_name="atac_names",
        columns=scmvp_atac_dataset.gene_names,
    )
    Ys.append(measurement)

    # make RNA object
    cell_attributes_dict = {"barcodes": rna_hvg.obs_names.to_numpy()}
    scmvp_rna_dataset.populate_from_data(
        X=rna_hvg.X.toarray(),
        batch_indices=None,
        gene_names=rna_hvg.var_names.to_numpy(),
        cell_attributes_dict=cell_attributes_dict,
        Ys=Ys,
    )

    # Pre Train
    lr = 5e-3
    use_batches = False
    use_cuda = True
    n_centroids = 15
    n_alfa = 1.0

    ## ATAC
    pre_atac_vae = VAE_Peak_SelfAttention(scmvp_atac_dataset.nb_genes,
                                          n_latent=20,
                                          n_batch=0,
                                          n_layers=1,
                                          log_variational=True,
                                          reconstruction_loss="nb")
    pre_atac_trainer = UnsupervisedTrainer(
        pre_atac_vae,
        scmvp_atac_dataset,
        train_size=0.9,
        use_cuda=use_cuda,
        frequency=5,
    )
    pre_atac_trainer.train(n_epochs=train_epochs[0], lr=lr)
    torch.save(pre_atac_trainer.model.state_dict(),
               '{}/pre_atac_trainer.pkl'.format(output_path))
    pre_atac_trainer.model.eval()

    ## ATAC pretrainer_posterior:
    full = pre_atac_trainer.create_posterior(pre_atac_trainer.model,
                                             scmvp_atac_dataset,
                                             indices=np.arange(len(scmvp_atac_dataset)))
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()

    ## RNA
    pre_vae = VAE_Attention(scmvp_rna_dataset.nb_genes,
                            n_latent=20,
                            n_batch=0,
                            n_layers=1,
                            log_variational=True,
                            reconstruction_loss="nb")
    pre_trainer = UnsupervisedTrainer(
        pre_vae,
        scmvp_rna_dataset,
        train_size=0.9,
        use_cuda=use_cuda,
        frequency=5,
    )
    pre_trainer.train(n_epochs=train_epochs[1], lr=lr)
    torch.save(pre_trainer.model.state_dict(),
               '{}/pre_trainer.pkl'.format(output_path))
    pre_trainer.model.eval()

    ## RNA pretrainer_posterior:
    full = pre_trainer.create_posterior(pre_trainer.model,
                                        scmvp_rna_dataset,
                                        indices=np.arange(len(scmvp_rna_dataset)))
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()
    imputed_values = full.sequential().imputation()

    # Reload pretrainer for multiVAE
    pre_trainer = UnsupervisedTrainer(
        pre_vae,
        scmvp_rna_dataset,
        train_size=0.9,
        use_cuda=use_cuda,
        frequency=5,
    )
    pre_trainer.model.load_state_dict(
        torch.load('%s/pre_trainer.pkl' % output_path))

    pre_atac_trainer = UnsupervisedTrainer(
        pre_atac_vae,
        scmvp_atac_dataset,
        train_size=0.9,
        use_cuda=use_cuda,
        frequency=5,
    )
    pre_atac_trainer.model.load_state_dict(
        torch.load('%s/pre_atac_trainer.pkl' % output_path))

    # Visulize Data
    prior_adata = anndata.AnnData(X=scmvp_rna_dataset.X)
    prior_adata.obsm["X_multi_vi"] = latent
    prior_adata.obs['cell_type'] = torch.tensor(labels.reshape(-1, 1))
    sc.pp.neighbors(prior_adata, use_rep="X_multi_vi", n_neighbors=30)
    sc.tl.umap(prior_adata, min_dist=0.3)
    sc.tl.louvain(prior_adata)

    n_centroids = len(np.unique(prior_adata.obs['louvain'].tolist()))

    # joint RNA and ATAC embedding
    multi_vae = Multi_VAE_Attention(
        scmvp_rna_dataset.nb_genes,
        len(scmvp_rna_dataset.atac_names),
        n_batch=0,
        n_latent=20,
        n_centroids=n_centroids,
        n_alfa=n_alfa,
        mode="mm-vae")  # should provide ATAC num, alfa, mode and loss type
    trainer = MultiTrainer(
        multi_vae,
        scmvp_rna_dataset,
        train_size=0.9,
        use_cuda=use_cuda,
        frequency=5,
    )

    trainer.model.init_gmm_params_with_louvain(
        latent,
        np.array(prior_adata.obs['louvain'].tolist()).astype(int))

    trainer.model.RNA_encoder.load_state_dict(
        pre_trainer.model.z_encoder.state_dict())
    for param in trainer.model.RNA_encoder.parameters():
        param.requires_grad = True
    trainer.model.ATAC_encoder.load_state_dict(
        pre_atac_trainer.model.z_encoder.state_dict())
    for param in trainer.model.ATAC_encoder.parameters():
        param.requires_grad = True
    trainer.train(n_epochs=train_epochs[2], lr=lr)

    torch.save(trainer.model.state_dict(),
               '%s/multi_vae_trainer.pkl' % output_path)
    trainer.model.eval()

    # LoadData
    scmvp_multi_data = scmvp_rna_dataset

    full = trainer.create_posterior(trainer.model,
                                    scmvp_multi_data,
                                    indices=np.arange(len(scmvp_multi_data)),
                                    type_class=MultiPosterior)
    latent, latent_rna, latent_atac, cluster_gamma, cluster_index, batch_indices, labels = full.sequential(
    ).get_latent()
    batch_indices = batch_indices.ravel()
    imputed_values = full.sequential().imputation()

    # Visulize Data
    prior_adata = anndata.AnnData(X=latent)
    prior_adata.obsm["X_multi_vi"] = latent
    prior_adata.obs['cell_type'] = torch.tensor(labels.reshape(-1, 1))
    sc.pp.neighbors(prior_adata, use_rep="X_multi_vi", n_neighbors=30)
    sc.tl.umap(prior_adata, min_dist=0.3)
    sc.tl.louvain(prior_adata)

    # Save Results
    ## save UMAP
    umap = pd.DataFrame(prior_adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], 
                        index=scmvp_rna_dataset.barcodes)
    umap.insert(2, "cluster", prior_adata.obs['louvain'].values)
    umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMVP-multi-umap.csv"))

    ## save latent
    latent = pd.DataFrame(prior_adata.obsm["X_multi_vi"], index=scmvp_rna_dataset.barcodes)
    latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMVP-multi-latent.csv"))

    # save imputation
    ## RNA
    imputation_rna = pd.DataFrame(imputed_values[0]) # imputation
    imputation_rna.columns = rna_hvg.var.index.values
    imputation_rna.index = rna_hvg.obs.index.values
    imputation_rna.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMVP-imputation-rna.csv.gz"), compression='gzip')

    ## ATAC
    imputation_atac = pd.DataFrame(imputed_values[1]) # imputation
    imputation_atac.columns = atac.var.index.values
    imputation_atac.index = atac.obs.index.values
    imputation_atac.to_csv(os.path.join(output_path, config["output_prefix"] + "-scMVP-imputation-atac.csv.gz"), compression='gzip')

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-scMVP-gpu_memory.csv'), header=["gpu_memory"])



# %%
scmvp_module(input_path = sys.argv[1],
             output_path = sys.argv[2],
             config = json.load(open(sys.argv[3]))
            )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/data/test/run_scMVP"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
# scmvp_module(input_path,
#              output_path,
#              config)

# %%

# %%
