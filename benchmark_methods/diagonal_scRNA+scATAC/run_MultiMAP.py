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
#     display_name: python (MultiMAP)
#     language: python
#     name: multimap
# ---

# %%
import os
import sys 
import json

import anndata
import MultiMAP
import scanpy as sc
import pandas as pd

# https://nbviewer.org/github/Teichlab/MultiMAP/blob/master/examples/tutorial.ipynb

# %%
def MultiMAP_module(input_path,
                   output_path,
                   config
                  ):
    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))
    gam = sc.read_h5ad(os.path.join(input_path, config['gam_h5ad_filename']))

    # Load Metadata
    metadata = pd.read_csv(os.path.join(input_path, config['metadata']), header=0)
    metadata.index=metadata[config['barcode_key']].values

    # preprocess
    ## ATAC
    MultiMAP.TFIDF_LSI(atac)
    gam.obsm['X_lsi'] = atac.obsm['X_lsi'].copy()
    ## RNA
    rna_pca = rna.copy()
    sc.pp.scale(rna_pca)
    sc.pp.pca(rna_pca)
    rna.obsm['X_pca'] = rna_pca.obsm['X_pca'].copy()
    # Integration
    adata = MultiMAP.Integration([rna, gam], ['X_pca', 'X_lsi'])

    # sc.pl.embedding(adata, 'X_multimap', color=['batch','cell_type'])

    latent = pd.DataFrame(adata.obsm['X_multimap'])
    latent.index = adata.obs_names

    # save latent
    ## rna
    rna_latent = latent.iloc[:int(latent.shape[0]/2),]
    rna_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-MultiMAP-RNA-latent.csv"))
    ## atac
    atac_latent = latent.iloc[int(latent.shape[0]/2):,]
    atac_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-MultiMAP-ATAC-latent.csv"))


# %%
MultiMAP_module(input_path = sys.argv[1], 
                output_path = sys.argv[2],
                config = json.load(open(sys.argv[3]))
               )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/test"
# output_path = "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_MultiMAP"
# config = json.load(open("/home/wsg/BM/data/test/c1k.json"))

# %%
# MultiMAP_module(input_path,
#                 output_path,
#                 config)

# %%
