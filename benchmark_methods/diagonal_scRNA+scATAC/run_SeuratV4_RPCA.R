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
#     display_name: R (r4)
#     language: R
#     name: r4
# ---

# %%
args <- commandArgs(T) 

# %%
library(here)
library(rjson)

library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)

# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac

# %%
SeuratV4_RPCA_module <- function(input_path,
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
    SeuratV4.rna <- CreateSeuratObject(counts = rna, project = "SeuratV4_RNA")
    # SeuratV4.rna[["percent.mt"]] <- PercentageFeatureSet(SeuratV4.rna, pattern = "^MT-")
    SeuratV4.rna$tech <- "rna"

    SeuratV4.rna <- NormalizeData(SeuratV4.rna, normalization.method = "LogNormalize", scale.factor = 10000)
    SeuratV4.rna <- FindVariableFeatures(SeuratV4.rna, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(SeuratV4.rna)
    SeuratV4.rna <- ScaleData(SeuratV4.rna, features = all.genes)

    SeuratV4.rna <- RunPCA(SeuratV4.rna, features = VariableFeatures(object = SeuratV4.rna))
    SeuratV4.rna <- FindNeighbors(SeuratV4.rna, dims = 1:20)
    SeuratV4.rna <- FindClusters(SeuratV4.rna, resolution = 0.2)
    ElbowPlot(SeuratV4.rna)
    SeuratV4.rna <- RunUMAP(SeuratV4.rna, dims = 1:10)
    # SeuratV4.rna <- RunTSNE(SeuratV4.rna, dims = 1:10)

    # ATAC Preprocess
    activity.matrix <- readRDS(here(input_path, config['gam_rds_filename']))

    SeuratV4.atac <- CreateSeuratObject(counts = atac, assay = "ATAC", project = "SeuratV4_ATAC")
    SeuratV4.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
    SeuratV4.atac$tech <- "atac"

    DefaultAssay(SeuratV4.atac) <- "ACTIVITY" 
    SeuratV4.atac <- FindVariableFeatures(SeuratV4.atac)
    SeuratV4.atac <- NormalizeData(SeuratV4.atac)
    SeuratV4.atac <- ScaleData(SeuratV4.atac)
    
    DefaultAssay(SeuratV4.atac) <- "ATAC"
    SeuratV4.atac <- RunTFIDF(SeuratV4.atac)
    SeuratV4.atac <- FindTopFeatures(SeuratV4.atac, min.cutoff = 'q0')
    SeuratV4.atac <- RunSVD(SeuratV4.atac)
    SeuratV4.atac <- RunUMAP(SeuratV4.atac, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    

    transfer.anchors <- FindTransferAnchors(reference = SeuratV4.rna, query = SeuratV4.atac, features = VariableFeatures(object = SeuratV4.rna), 
                                            reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

    genes.use <- VariableFeatures(SeuratV4.rna)
    refdata <- GetAssayData(SeuratV4.rna, assay = "RNA", slot = "data")[genes.use, ]

    # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
    # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = SeuratV4.atac[["lsi"]], dims = 1:10)

    # this line adds the imputed data matrix to the SeuratV4.atac object
    SeuratV4.atac[["RNA"]] <- imputation
    coembed <- merge(x = SeuratV4.rna, y = SeuratV4.atac)

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
                                               "SeuratV4_RPCA-RNA-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(obj.list[["rna"]]@reductions$pca@cell.embeddings,
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4_RPCA-RNA-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    # ATAC
    umap <- as.data.frame(obj.list[["atac"]]@reductions$umap@cell.embeddings)
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap[, "cluster"] <- obj.list[["atac"]]@meta.data$seurat_clusters
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4_RPCA-ATAC-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(obj.list[["atac"]]@reductions$pca@cell.embeddings,
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4_RPCA-ATAC-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %%
SeuratV4_RPCA_module(input_path = args[1],
                     output_path = args[2],
                     config = unlist(fromJSON(file=args[3]))
                    )

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_SeuratV4_RPCA"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/RawData.json"))

# %%
# SeuratV4_RPCA_module(input_path,
#                      output_path,
#                      config)
