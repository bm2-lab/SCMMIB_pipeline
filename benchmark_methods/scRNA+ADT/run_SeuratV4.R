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

# %%
library(here)
library(rjson)
library(Matrix)
library(tidyverse)
library(Seurat)

# %%
SeuratV4_module <- function(input_path,
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

    # Integrate Data
    multi <- CreateSeuratObject(counts = rna)
    DefaultAssay(multi) <- "RNA"
    multi[["ADT"]] <- CreateAssayObject(counts = adt)
    ## RNA
    DefaultAssay(multi) <- 'RNA'
    multi <- NormalizeData(multi) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
    ## ADT
    DefaultAssay(multi) <- 'ADT'
    # we will use all ADT features for dimensional reduction
    VariableFeatures(multi) <- rownames(multi[["ADT"]])
    multi <- NormalizeData(multi, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca')

    # Cluster
    multi <- FindMultiModalNeighbors(
      multi, reduction.list = list("pca", "apca"), 
      dims.list = list(1:30, 1:min(18, length(rownames(multi[["ADT"]]))-1)), modality.weight.name = "RNA.weight"
    )
    multi <- RunUMAP(multi, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    multi <- FindClusters(multi, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

    # Save Results
    ## UMAP
    umap <- as.data.frame(multi@reductions$wnn.umap@cell.embeddings)
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap[, "cluster"] <- multi@meta.data$seurat_clusters
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4-multi-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## graph
    write.table(multi@graphs$wsnn,
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4-multi-graph.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %%
SeuratV4_module(input_path = args[1],
                output_path = args[2],
                config = unlist(fromJSON(file=args[3]))
               )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData"
# output_path <- "/home/wsg/BM/results/task/scRNA+ADT/accuracy/10x_NSCLC/RawData/rep_1/run_SeuratV4"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData/RawData.json"))

# %%
# SeuratV4_module(input_path,
#                 output_path, 
#                 config)

# %%
