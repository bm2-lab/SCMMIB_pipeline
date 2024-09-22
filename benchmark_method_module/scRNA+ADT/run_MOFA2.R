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

# %%
mofa_module <- function(input_path, 
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

    # SCTransform on RNA
    mofa_data <- CreateSeuratObject(counts = rna)
    DefaultAssay(mofa_data) <- "RNA"
    mofa_data <- SCTransform(mofa_data, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    mofa_data <- FindVariableFeatures(mofa_data, selection.method = "vst", nfeatures = 3000) 

    # Add ADT
    mofa_data[["ADT"]] <- CreateAssayObject(counts = adt)

    ## ADT
    DefaultAssay(mofa_data) <- 'ADT'
    # we will use all ADT features for dimensional reduction
    VariableFeatures(mofa_data) <- rownames(mofa_data[["ADT"]])
    mofa_data <- NormalizeData(mofa_data, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca')

    # Filter Data
    #     mofa_data <- subset(x = mofa_data, 
    #                         subset = nCount_RNA < 25000 & nCount_RNA > 10)

    # Merge Data
    mofa <- create_mofa(mofa_data, assays = c("SCT","ADT"))
    plot_data_overview(mofa)

    #     mofa <- readRDS(here(input_path, 
    #                          paste(config["data_name"], 
    #                                "MOFA2", 
    #                                config["task_type"],
    #                                "multi", "filtered.rds", sep = "-")))

    # Define options
    ## Define data options
    data_opts <- get_default_data_options(mofa)
    ## Define model options
    model_opts <- get_default_model_options(mofa)
    model_opts$num_factors <- 15
    ## Define train options
    train_opts <- get_default_training_options(mofa)
    train_opts$gpu_mode=TRUE

    # stochastic：SGD，适合>10000个细胞和GPU。
    #     train_opts$stochastic=TRUE
    #     stochastic_opts <- get_default_stochastic_options(mofa)
    #     stochastic_opts$batch_size<-0.25

    # Build and train Data
    mofa <- prepare_mofa(object = mofa,
                         data_options = data_opts,
                         model_options = model_opts,
                         training_options= train_opts
    ##                      stochastic_options = stochastic_opts
                        )
    mofa <- run_mofa(mofa, 
                     outfile =here(output_path, "mofa_model.hdf5") ,
                     use_basilisk =  TRUE)

    model <- load_model(here(output_path, "mofa_model.hdf5"),
                       remove_inactive_factors = FALSE)

    # Visualize Data
    set.seed(42)
    model <- run_umap(model,
                      n_neighbors = 20, 
                      min_dist = 0.30
    )

    # Cluster Data
    mofa_cluster <- cluster_samples(mofa, k=20, 
                                    factors = "all")
    mofa_cluster <- as.data.frame(mofa_cluster$cluster)
    colnames(mofa_cluster) <- "cluster"

    # Save Results
    ## prepare UMAP
    umap <- cbind(model@dim_red$UMAP, mofa_cluster)
    umap <- umap[, -1]

    ## plot_dimred(model, method = "UMAP")

    # Save Results
    ## save UMAP
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "MOFA2-multi-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
    # save latent
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
# input_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1"
# output_path <- "/home/wsg/BM/results/task/scRNA+ADT/accuracy/SPATIAL/lymph_node_A1/rep_1/run_MOFA2"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node_A1/lymph_node.json"))

# %%
# mofa_module(input_path,
#             output_path,
#             config)
