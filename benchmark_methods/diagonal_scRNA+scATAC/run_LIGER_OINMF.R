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

# http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/online_iNMF_tutorial.html


# %%
LIGER_OINMF_module <- function(input_path,
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

    # RNA processing
    rna = createLiger(list(rna = rna))
    rna = normalize(rna)
    rna = selectGenes(rna, var.thresh = 0.1)

    # ATAC processing
    atac = createLiger(list(atac = activity.matrix))
    atac = normalize(atac)

    # filter rna name not in atac
    rna@var.genes = rna@var.genes[rna@var.genes %in% rownames(activity.matrix)]
    rna_name_notin_atac <- rownames(rna@raw.data$rna)[!(rownames(rna@raw.data$rna) %in% rownames(atac@raw.data$atac))]
    rna@var.genes <- rna@var.genes[!(rna@var.genes %in% rna_name_notin_atac)]

    rna.vargenes = rna@var.genes
    rna = scaleNotCenter(rna)

    atac@var.genes = rna@var.genes
    atac = scaleNotCenter(atac)

    # rna = online_iNMF(rna, k = 20, max.epochs = 1)
    rna = online_iNMF(rna, k = 20, max.epochs = 1, miniBatch_size = 300) # 细胞数量低于5000时，需要降低miniBatch_size
    rna = quantile_norm(rna)
    rna = runUMAP(rna)
    #     plotByDatasetAndCluster(rna, axis.labels = c("UMAP1","UMAP2"))

    # online_iNMF
    # liger = online_iNMF(rna, X_new = list(atac = atac), k = 20, max.epochs = 1)
    liger = online_iNMF(rna, X_new = list(atac = atac), k = 20, max.epochs = 1, miniBatch_size = 300) # 细胞数量低于5000时，需要降低miniBatch_size

    liger = quantile_norm(liger)
    liger = runUMAP(liger)
    # plotByDatasetAndCluster(liger, axis.labels = c("UMAP1","UMAP2"))

    # Save Results
    ## UMAP
    umap <- as.data.frame(liger@tsne.coords)
    umap[, "cluster"] <- liger@clusters

    umap_RNA <- umap[which(sub(".*_", "", rownames(umap)) == 'RNA'), ]
    rownames(umap_RNA) <- sub("_.*", "", rownames(umap_RNA))
    colnames(umap_RNA) <- c("UMAP1", "UMAP2", "cluster")
    write.table(umap_RNA, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_OINMF-RNA-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    umap_ATAC <- umap[which(sub(".*_", "", rownames(umap)) == 'ATAC'), ]
    rownames(umap_ATAC) <- sub("_.*", "", rownames(umap_ATAC))
    colnames(umap_ATAC) <- c("UMAP1", "UMAP2", "cluster")
    write.table(umap_ATAC, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_OINMF-ATAC-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## latent
    rna_latent <- as.data.frame(liger@H$rna)
    rownames(rna_latent) <- sub("_.*", "", rownames(rna_latent))
    write.table(rna_latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_OINMF-RNA-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    atac_latent <- as.data.frame(liger@H$atac)
    rownames(atac_latent) <- sub("_.*", "", rownames(atac_latent))
    write.table(atac_latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_OINMF-ATAC-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    # write.table(liger@H.norm,
    #             file = here(output_path, paste(config["output_prefix"], 
    #                                            "LIGER_OINMF-multi-latent.csv", sep = "-")), 
    #             row.names =T, col.names = T, sep=',', quote=F)

}

# %%
LIGER_OINMF_module(input_path = args[1],
                   output_path = args[2],
                   config = unlist(fromJSON(file=args[3]))
                  )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/test"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_LIGER_OINMF"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/test/c1k.json"))

# %%
# LIGER_OINMF_module(input_path,
#                    output_path,
#                    config)
