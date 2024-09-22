# -*- coding: utf-8 -*-
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
#     display_name: R (SeuratV3)
#     language: R
#     name: seuratv3
# ---

# %%
args <- commandArgs(T) 

# %%
library(here)
library(rjson)

library(rliger)
library(Seurat)
library(stringr)

# http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/SNAREseq_walkthrough.html
# http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html


# %%
LIGER_UINMF_module <- function(input_path,
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
    atac <- readRDS(here(input_path, config["atac_rds_filename"]))

    metadata <- read.csv(here(input_path, config["metadata"]), row.names = config["barcode_key"])

    # Get GAM of ATAC
    activity.matrix <- readRDS(here(input_path, config["gam_rds_filename"]))


    # rename cell barcodes
    colnames(rna) <- paste0(colnames(rna), "_RNA")
    colnames(activity.matrix) <- paste0(colnames(activity.matrix), "_ATAC")

    # split the GAM by intersection
    shared_atac <- activity.matrix[rownames(activity.matrix) %in% rownames(rna),]
    unshared_atac <- activity.matrix[!(rownames(activity.matrix) %in% rownames(rna)),]

    # unshared atac
    ## select unshared gene-centric features
    liger <- createLiger(list(rna = rna, atac = unshared_atac))
    liger <- normalize(liger)
    liger <- selectGenes(liger, unshared = TRUE, unshared.datasets = list(2), unshared.thresh= 0.2)
    #     length(liger@var.unshared.features$atac)
    liger <- scaleNotCenter(liger)
    #     dim(liger@scale.unshared.data$atac)
    ## get preprocessed unshared atac features
    unshared_feats <- liger@scale.unshared.data$atac

    # shared atac
    # 1. Create a LIGER object and normalize the shared data
    liger <- createLiger(list(rna = rna, atac = activity.matrix))
    liger <- normalize(liger)
    # 2. use the RNA dataset to select variable shared features
    liger <- selectGenes(liger, var.thresh = 0.1, datasets.use =1 , unshared = TRUE,  unshared.datasets = list(2), unshared.thresh= 0.2)
    # liger <- selectGenes(liger, var.thresh = 0.1, datasets.use =1)
    # 3. Scale
    liger <- scaleNotCenter(liger)

    # 4. Add the unshared features that have been properly selected as genes
    unshared_feats_names <- rownames(unshared_feats)
    liger@var.unshared.features[[2]] = unshared_feats_names
    liger@scale.unshared.data[[2]] = unshared_feats

    # Joint Matrix Factorization
    ## 1. To factorize the datasets and include the unshared datasets, set the use.unshared parameter to TRUE.
    liger <- optimizeALS(liger, k=30, use.unshared = TRUE, max_iters =30, thresh=1e-10)

    # Quantile Normalization and Joint Clustering
    liger <- quantile_norm(liger)
    liger <- louvainCluster(liger)

    # Visualizations and Save results
    liger <- runUMAP(liger)
    umap_plots <-plotByDatasetAndCluster(liger, axis.labels = c("UMAP1","UMAP2"), return.plots = TRUE)
    #     umap_plots[[2]]
    ## UMAP
    umap <- as.data.frame(liger@tsne.coords)
    umap[, "cluster"] <- liger@clusters

    umap_RNA <- umap[which(sub(".*_", "", rownames(umap)) == 'RNA'), ]
    rownames(umap_RNA) <- sub("_.*", "", rownames(umap_RNA))
    colnames(umap_RNA) <- c("UMAP1", "UMAP2", "cluster")
    write.table(umap_RNA, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_UINMF-RNA-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)


    umap_ATAC <- umap[which(sub(".*_", "", rownames(umap)) == 'ATAC'), ]
    rownames(umap_ATAC) <- sub("_.*", "", rownames(umap_ATAC))
    colnames(umap_ATAC) <- c("UMAP1", "UMAP2", "cluster")
    write.table(umap_ATAC, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_UINMF-ATAC-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## latent
    rna_latent <- as.data.frame(liger@H$rna)
    rownames(rna_latent) <- sub("_.*", "", rownames(rna_latent))
    write.table(rna_latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_UINMF-RNA-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    atac_latent <- as.data.frame(liger@H$atac)
    rownames(atac_latent) <- sub("_.*", "", rownames(atac_latent))
    write.table(atac_latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_UINMF-ATAC-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)


    # write.table(liger@H.norm,
    #             file = here(output_path, paste(config["output_prefix"], 
    #                                            "LIGER_UINMF-multi-latent.csv", sep = "-")), 
    #             row.names =T, col.names = T, sep=',', quote=F)

}

# %%
LIGER_UINMF_module(input_path = args[1],
                   output_path = args[2],
                   config = unlist(fromJSON(file=args[3]))
                  )

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/10x_mouse_brain/2p5/rep_1/run_LIGER_UINMF"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/2p5.json"))

# %%
# LIGER_UINMF_module(input_path,
#                    output_path,
#                    config)

# %%
