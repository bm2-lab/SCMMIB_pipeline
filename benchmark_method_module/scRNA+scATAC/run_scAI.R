# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: R (r4)
#     language: R
#     name: r4
# ---

# %% vscode={"languageId": "r"}
args <- commandArgs(T) 

# %% vscode={"languageId": "r"}
library(here)
library(rjson)
library(Matrix)

library(scAI)
# https://htmlpreview.github.io/?https://github.com/sqjin/scAI/blob/master/examples/walkthrough_Kidneydataset.html

# %%
scAI_module<-function(input_path,
                      output_path, 
                      config){
    # Make Dir
    if (!dir.exists(output_path)){
        dir.create(output_path, 
                   recursive = TRUE)
    }

    # Load Data
    rna <- readRDS(here(input_path, config["rna_rds_filename"]))
    atac <- readRDS(here(input_path, config["atac_rds_filename"]))

    metadata <- read.csv(here(input_path, config["metadata"]), row.names = config["barcode_key"])

    # Merge Data
    X = list(RNA = rna, ATAC = atac)

    # Create a scAI object
    scAI_outs <- create_scAIobject(raw.data = X)
    scAI_outs <- addpData(scAI_outs, pdata = metadata, pdata.name = colnames(metadata))

    # Data preprocess
    # Perform quality control to remove low-quality cells and genes, and normalize the data.
    scAI_outs <- preprocessing(scAI_outs, minFeatures = 0, minCells = 0)

    # Run scAI model
    scAI_outs <- run_scAI(scAI_outs, K = 20, 
                          nrun = 1, do.fast = T)

    # Identify cell clusters
    scAI_outs <- identifyClusters(scAI_outs, resolution = 0.9) 
    scAI_outs <- getAggregatedData(scAI_outs, group = scAI_outs@identity)

    # Visualize cells onto the low-dimensional space
    scAI_outs <- reducedDims(scAI_outs, method = "umap")
    # cellVisualization(scAI_outs, scAI_outs@embed$umap, color.by = "cluster")

    ## save UMAP
    write.table(scAI_outs@embed$umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "scAI-multi-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## save latent
    write.table(t(scAI_outs@fit$H),
                file = here(output_path, paste(config["output_prefix"], 
                                               "scAI-multi-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %% vscode={"languageId": "r"}
scAI_module(input_path = args[1],
            output_path = args[2],
            config = unlist(fromJSON(file=args[3]))
           )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/test"
# output_path <- "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/HSPC/p10/rep_1/run_scAI"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/test/c1k.json"))

# %%
# scAI_module(input_path,
#             output_path,
#             config)

# %%
