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

# %%
args <- commandArgs(T) 

# %% vscode={"languageId": "python"}
library(here)
library(rjson)
library(dplyr)
library(tidyverse)

library(Matrix)
library(stringr)
library(cowplot)
library(ggplot2)
library(data.table)

library(MOFA2)
library(Seurat)
library(SeuratDisk)
library(magrittr)

# https://biofam.github.io/MOFA2/MEFISTO.html
# https://raw.githack.com/bioFAM/MEFISTO_tutorials/master/MEFISTO_spatial.html

# %%
MEFISTO_module <- function(input_path,
                           output_path,
                           config
                           ){
    # Make Dir
    if (!dir.exists(output_path)){
        dir.create(output_path, 
                   recursive = TRUE)
    }

    # Load Data
    rna <- readRDS(here(input_path, config["rna_rds_filename"]))
    adt <- readRDS(here(input_path, config["adt_rds_filename"]))

    metadata <- read.csv(here(input_path, config["metadata"]), row.names = config["barcode_key"])

    multi <- CreateSeuratObject(counts = rna, meta.data = metadata)
    multi[["ADT"]] <- CreateAssayObject(counts = adt)

    DefaultAssay(multi) <- "RNA"
    multi  <- SCTransform(multi , verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    multi <- FindVariableFeatures(multi, selection.method = "vst", nfeatures = 3000) 

    DefaultAssay(multi) <- "ADT"
    multi <- NormalizeData(multi, normalization.method = 'CLR', margin = 2) %>% ScaleData() 
    multi<- FindVariableFeatures(multi)

    mofa <- create_mofa(multi, assays = c("SCT","ADT"))

    coords <- t(multi@meta.data[,c("X", "Y")])
    rownames(coords) <- c("coord_1", "coord_2")

    mofa <- set_covariates(mofa, coords)

    data_opts <- get_default_data_options(mofa)

    model_opts <- get_default_model_options(mofa)
    model_opts$num_factors <- 4

    train_opts <- get_default_training_options(mofa)
    train_opts$maxiter <- 1

    mefisto_opts <- get_default_mefisto_options(mofa)

    mofa <- prepare_mofa(mofa, model_options = model_opts,
                       mefisto_options = mefisto_opts,
                       training_options = train_opts,
                       data_options = data_opts)

    mofa <- run_mofa(mofa, outfile=here(output_path, "mefisto_model.hdf5"), use_basilisk = TRUE)

    mofa <-load_model(here(output_path, "mefisto_model.hdf5"), remove_inactive_factors = FALSE)
    factors <- 1:get_dimensions(mofa)[["K"]]

    mofa <- run_umap(mofa, factors = factors, n_neighbors = 15, min_dist = 0.30)

    ## save latent
    write.table(mofa@expectations$Z$group1,
                file = here(output_path, paste(config["output_prefix"], 
                                               "MEFISTO-multi-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %%
MEFISTO_module(input_path = args[1],
               output_path = args[2],
               config = unlist(fromJSON(file=args[3]))
              )

# %%

# %% vscode={"languageId": "r"}

# %% vscode={"languageId": "r"}

# %%
# input_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1"
# output_path <- "/home/wsg/BM/results/task/spatial_scRNA+ADT/accuracy/lymph_node_A1/rep_1/run_MEFISTO"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/lymph_node.json"))

# %%
# MEFISTO_module(input_path,
#                output_path,
#                config)

# %%
