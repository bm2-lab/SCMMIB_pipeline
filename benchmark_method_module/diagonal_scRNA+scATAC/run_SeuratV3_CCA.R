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

library(Seurat)
library(Matrix)
library(ggplot2)
# https://satijalab.org/seurat/archive/v3.2/atacseq_integration_vignette

# %%
SeuratV3_module <- function(input_path,
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

    # RNA Preprocess
    SeuratV3.rna <- CreateSeuratObject(counts = rna, project = "SeuratV3_RNA")
    # SeuratV3.rna[["percent.mt"]] <- PercentageFeatureSet(SeuratV3.rna, pattern = "^MT-")
    SeuratV3.rna$tech <- "rna"

    SeuratV3.rna <- NormalizeData(SeuratV3.rna, normalization.method = "LogNormalize", scale.factor = 10000)
    SeuratV3.rna <- FindVariableFeatures(SeuratV3.rna, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(SeuratV3.rna)
    SeuratV3.rna <- ScaleData(SeuratV3.rna, features = all.genes)

    SeuratV3.rna <- RunPCA(SeuratV3.rna, features = VariableFeatures(object = SeuratV3.rna))
    SeuratV3.rna <- FindNeighbors(SeuratV3.rna, dims = 1:20)
    SeuratV3.rna <- FindClusters(SeuratV3.rna, resolution = 0.2)
    ElbowPlot(SeuratV3.rna)
    SeuratV3.rna <- RunUMAP(SeuratV3.rna, dims = 1:10)
    # SeuratV3.rna <- RunTSNE(SeuratV3.rna, dims = 1:10)

    # ATAC Preprocess
    activity.matrix <- readRDS(here(input_path, config['gam_rds_filename']))

    SeuratV3.atac <- CreateSeuratObject(counts = atac, assay = "ATAC", project = "SeuratV3_ATAC")
    SeuratV3.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
    SeuratV3.atac$tech <- "atac"

    DefaultAssay(SeuratV3.atac) <- "ACTIVITY" 
    SeuratV3.atac <- FindVariableFeatures(SeuratV3.atac)
    SeuratV3.atac <- NormalizeData(SeuratV3.atac)
    SeuratV3.atac <- ScaleData(SeuratV3.atac)

    DefaultAssay(SeuratV3.atac) <- "ATAC"
    VariableFeatures(SeuratV3.atac) <- names(which(Matrix::rowSums(SeuratV3.atac) > 25))
    SeuratV3.atac <- RunLSI(SeuratV3.atac, n = 50, scale.max = NULL)
    SeuratV3.atac <- RunUMAP(SeuratV3.atac, reduction = "lsi", dims = 1:15)

    transfer.anchors <- FindTransferAnchors(reference = SeuratV3.rna, query = SeuratV3.atac, features = VariableFeatures(object = SeuratV3.rna), 
                                            reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

    genes.use <- VariableFeatures(SeuratV3.rna)
    refdata <- GetAssayData(SeuratV3.rna, assay = "RNA", slot = "data")[genes.use, ]

    # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
    # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = SeuratV3.atac[["lsi"]], dims = 1:10)

    # this line adds the imputed data matrix to the SeuratV3.atac object
    SeuratV3.atac[["RNA"]] <- imputation
    coembed <- merge(x = SeuratV3.rna, y = SeuratV3.atac)

    # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
    # configs
    coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
    coembed <- RunUMAP(coembed, dims = 1:30)

    obj.list <- SplitObject(coembed, split.by = "tech")

    obj.list[["rna"]] <- RenameCells(
        obj.list[["rna"]],
        new.names = sub("_.*", "", colnames(x = obj.list[["rna"]]))
    )

    obj.list[["atac"]] <- RenameCells(
        obj.list[["atac"]],
        new.names = sub("_.*", "", colnames(x = obj.list[["rna"]]))
    )

    # RNA
    umap <- as.data.frame(obj.list[["rna"]]@reductions$umap@cell.embeddings)
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap[, "cluster"] <- obj.list[["rna"]]@meta.data$seurat_clusters
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV3_CCA-RNA-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(obj.list[["rna"]]@reductions$pca@cell.embeddings,
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV3_CCA-RNA-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    # ATAC
    umap <- as.data.frame(obj.list[["atac"]]@reductions$umap@cell.embeddings)
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap[, "cluster"] <- obj.list[["atac"]]@meta.data$seurat_clusters
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV3_CCA-ATAC-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(obj.list[["atac"]]@reductions$pca@cell.embeddings,
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV3_CCA-ATAC-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %%
SeuratV3_module(input_path = args[1],
                output_path = args[2],
                config = unlist(fromJSON(file=args[3]))
               )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_SeuratV3_CCA"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/RawData.json"))

# %%
# SeuratV3_module(input_path,
#                 output_path,
#                 config
#                )
