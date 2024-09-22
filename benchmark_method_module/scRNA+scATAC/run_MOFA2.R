# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: Rmd,R:percent
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

# %%
library(here)
library(rjson)
library(Matrix)

library(MOFA2)
library(Seurat)
library(Signac)
library(tidyverse)

# https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html

# %%
mofa_module <- function(input_path, 
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

    multi_data <- CreateSeuratObject(counts = rna, assay = "RNA", meta.data = metadata)

    chrom_assay <- CreateChromatinAssay(counts = atac, sep = c("-", "-"))
    multi_data[["ATAC"]] <- chrom_assay

    # Data preprocess
    # Normalize
    ## RNA
    DefaultAssay(multi_data) <- "RNA"
    multi_data <- NormalizeData(multi_data, normalization.method = "LogNormalize", assay = "RNA")
    multi_data <- ScaleData(multi_data, do.center = TRUE, do.scale = FALSE, assay = "RNA")
    ## ATAC
    DefaultAssay(multi_data) <- "ATAC"
    multi_data <- NormalizeData(multi_data)
    multi_data <- ScaleData(multi_data)
    multi_data <- RunTFIDF(multi_data, assay = "ATAC")

    # Feature selection
    ## RNA
    DefaultAssay(multi_data) <- "RNA"
    multi_data <- FindVariableFeatures(multi_data, selection.method = "vst", nfeatures = 3000, assay = "RNA") 
    ## ATAC
    DefaultAssay(multi_data) <- "ATAC"
    # FindTopFeatures(multi_data, min.cutoff = 2000, assay = "ATAC")
    multi_data <- FindTopFeatures(multi_data, min.cutoff = 'q80', assay = "ATAC")

    # Merge Data
    mofa <- create_mofa(multi_data, assays = c("RNA","ATAC"))
    # plot_data_overview(mofa)

    # Define options
    ## Define data options
    data_opts <- get_default_data_options(mofa)
    ## Define model options
    model_opts <- get_default_model_options(mofa)
    model_opts$num_factors <- 15
    ## Define train options
    train_opts <- get_default_training_options(mofa)
    train_opts$gpu_mode=TRUE

    # stochastic：SGD
    # train_opts$stochastic=TRUE
    # stochastic_opts <- get_default_stochastic_options(mofa)
    # stochastic_opts$batch_size<-0.25

    # Build and train Data
    mofa <- prepare_mofa(object = mofa,
                         data_options = data_opts,
                         model_options = model_opts,
                         training_options= train_opts
                         # stochastic_options = stochastic_opts
                        )
    mofa <- run_mofa(mofa, 
                     outfile =here(output_path, "mofa_model.hdf5") ,
                     use_basilisk =  TRUE)

    model <- load_model(here(output_path, "mofa_model.hdf5"),
                       remove_inactive_factors = FALSE)

    # Visualize Data
    model <- run_umap(model,
                      n_neighbors = 20, 
                      min_dist = 0.30
    )

    ## plot_dimred(model, method = "UMAP")

    # Cluster Data
    mofa_cluster <- cluster_samples(mofa, k=20, 
                                    factors = "all")
    mofa_cluster <- as.data.frame(mofa_cluster$cluster)
    colnames(mofa_cluster) <- "cluster"

    # Save Results
    ## prepare UMAP
    umap <- cbind(model@dim_red$UMAP, mofa_cluster)
    umap <- umap[, -1]

    ## save UMAP
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "MOFA2-multi-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## save latent
    write.table(model@expectations$Z$group1,
                file = here(output_path, paste(config["output_prefix"], 
                                               "MOFA2-multi-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %%
mofa_module(input_path = args[1],
            output_path = args[2],
            config = unlist(fromJSON(file=args[3]))
           )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/test"
# output_path <- "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/HSPC/p10/rep_1/run_MOFA2"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/test/c1k.json"))

# %%
# mofa_module(input_path,
#             output_path,
#             config)

# %%
