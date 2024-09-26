import argparse
import subprocess
import os


import scanpy as sc
import anndata
import MultiMAP
import pandas as pd
import muon
import numpy as np
sc.settings.set_figure_params(dpi=80)
from glob import glob
import sys
import json


def MultiMAP_module(input_path,
                   output_path,
                   config
                  ):
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
    adata_omics1 = rna
    adata_omics2 = adt
    # adata_omics1 = sc.read_h5ad(f'{file_fold}/adata_RNA.h5ad')
    # adata_omics2 = sc.read_h5ad(f'{file_fold}/adata_ADT.h5ad')

    adata_omics1.var_names_make_unique()
    adata_omics2.var_names_make_unique()
    rna = adata_omics1
    adt = adata_omics2
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)


    rna_pca = rna.copy()
    sc.pp.pca(rna_pca,n_comps=10)
    rna.obsm['X_pca'] = rna_pca.obsm['X_pca'].copy()

    adt_pca = adt.copy()
    muon.prot.pp.clr(adt_pca)

    sc.pp.pca(adt_pca,n_comps=10)
    adt.obsm['X_pca'] = adt_pca.obsm['X_pca'].copy()
    adata1 = MultiMAP.matrix.MultiMAP([rna.obsm['X_pca'],adt.obsm['X_pca']])[2] 
    adata2 = MultiMAP.matrix.MultiMAP([adata1,rna.obsm['spatial']])[2] 

    n = len(rna.obs.index)
    tmp1 = adata2[:n,:]
    tmp2 = adata2[n:2*n,:]
    tmp3 = adata2[2*n:,:]

    latent = pd.DataFrame(np.c_[tmp1,tmp2,tmp3],index=rna.obs.index)
    latent_embed_path = os.path.join(output_path, config["output_prefix"] + "-MultiMAP-multi-latent.csv")
    latent.to_csv(latent_embed_path)

MultiMAP_module(input_path = sys.argv[1], 
                output_path = sys.argv[2],
                config = json.load(open(sys.argv[3]))
               )
