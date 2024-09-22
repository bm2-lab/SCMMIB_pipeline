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
#     display_name: python (scglue)
#     language: python
#     name: scglue
# ---

# %%
import anndata as ad
import sys
import json
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
# from scMVP.dataset import LoadData
import numpy as np
import pandas as pd
import scipy.io as sp_io
from scipy.sparse import csr_matrix
import itertools
from itertools import chain
import seaborn as sns
import os
import torch
import pyranges as pr

# %%
from nvitop import Device
devices = Device.all()
memory_free = [device.memory_free() for device in devices]
gpu_index = memory_free.index(max(memory_free))
torch.cuda.set_device(gpu_index)


# %%
def glue_module(config,
                input_path,
                output_path,
                max_cpu=30):
    torch.set_num_threads(5)

    # Make Dir
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Load data
    rna = sc.read_h5ad(os.path.join(input_path, config['rna_h5ad_filename']))
    atac = sc.read_h5ad(os.path.join(input_path, config['atac_h5ad_filename']))

    # Get the GTF infor
    gtf = pr.read_gtf(config['gtf_file'])
    if config["specie"] == "human":
        chr_normal = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    elif config["specie"] == "mouse":
        chr_normal = ["chr" + str(i) for i in range(1, 20)] + ["chrX", "chrY"]
    gtf_gene_index = (gtf.Feature == "gene") & (gtf.Chromosome.isin(chr_normal))
    gtf_gene = gtf[gtf_gene_index]

    # filter gene
    rna = rna[:, rna.var.features.isin(gtf_gene.gene_name.tolist())]

    # RNA Preprocess
    rna.layers["counts"] = rna.X.copy()
    #     sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor="seurat_v3", batch_key='batch')
    sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor="seurat_v3")

    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")

    ## annotation
    scglue.data.get_gene_annotation(
        rna,
        gtf=config['gtf_file'],
        gtf_by='gene_name',
    )
    rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

    # ATAC Preprocess
    scglue.data.lsi(atac, n_components=100, n_iter=15)

    # Get genome location
    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    atac.var.head()

    # get graph
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)

    # check
    scglue.graph.check_graph(guidance, [rna, atac])

    # filter hvf
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )

    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"
    )

    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()

    # train model
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": "glue"}
    )

    # Save Results
    glue.save("{}/glue.dill".format(output_path))

    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    combined = ad.concat([rna, atac])
    sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
    sc.tl.umap(combined)

    feature_embeddings = glue.encode_graph(guidance_hvf)
    feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
    rna_cells = rna.obsm["X_glue"].shape[0]
    atac_cells = atac.obsm["X_glue"].shape[0]

    sc.tl.louvain(combined)

    # latent
    ## RNA
    rna_latent = pd.DataFrame(data=combined.obsm["X_glue"][0:rna_cells, :],
                      index=rna.obs.barcode)
    rna_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-GLUE-RNA-latent.csv"))
    ## ATAC
    atac_latent = pd.DataFrame(data=combined.obsm["X_glue"][rna_cells:rna_cells + atac_cells, :],
                      index=atac.obs.barcode)
    atac_latent.to_csv(os.path.join(output_path, config["output_prefix"] + "-GLUE-ATAC-latent.csv"))

    # UMAP
    ## RNA
    rna_umap = pd.DataFrame(combined.obsm["X_umap"][0:rna_cells, :], columns=["UMAP1", "UMAP2"], index=rna.obs.barcode)
    rna_umap.insert(2, "cluster", combined.obs['louvain'].values[0:rna_cells])
    rna_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-GLUE-RNA-umap.csv"))

    ## ATAC
    atac_umap = pd.DataFrame(combined.obsm["X_umap"][rna_cells:rna_cells + atac_cells, :], columns=["UMAP1", "UMAP2"],
                             index=atac.obs.barcode)
    atac_umap.insert(2, "cluster", combined.obs['louvain'].values[rna_cells:rna_cells + atac_cells])
    atac_umap.to_csv(os.path.join(output_path, config["output_prefix"] + "-GLUE-ATAC-umap.csv"))

    # gpu_memory
    pid= os.getpid()        
    gpu_memory = pd.Series(dtype='str')

    devices = Device.all()
    for device in devices:
        processes = device.processes()    
        if pid in processes.keys():
            p=processes[pid]
            gpu_memory['device ' + str(device.index)] = p.gpu_memory_human()

    gpu_memory.to_csv(os.path.join(output_path, config['output_prefix'] + '-GLUE-gpu_memory.csv'), header=["gpu_memory"])


# %%
glue_module(input_path = sys.argv[1], 
            output_path = sys.argv[2],
            config = json.load(open(sys.argv[3]))
           )

# %%

# %%

# %%

# %%
# input_path = "/home/wsg/BM/data/BMMC/RNA+ATAC/p10"
# output_path = "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/BMMC/p10_GAM/rep_1/run_GLUE"
# config = json.load(open("/home/wsg/BM/data/BMMC/RNA+ATAC/p10/p10.json"))

# %%
# glue_module(input_path = input_path, 
#             output_path = output_path,
#             config = config)

# %%
