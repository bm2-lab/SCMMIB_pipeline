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
library(Matrix)
library(peakRAM)

library(bindSC)
library(Seurat)
library(Signac)
library(pbapply)
library(future)

library(tidyverse)

library(densityClust)
library(cluster)

# %%
bindSC_module <- function(input_path,
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
    multi <- CreateSeuratObject(counts = rna)
    DefaultAssay(multi) <- "RNA"
    multi <- SCTransform(multi, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    multi <- FindVariableFeatures(multi, selection.method = "vst", nfeatures = 3000)

    multi_hvgs <- VariableFeatures(multi)
    rna_normalize <- as.matrix(multi[["SCT"]]@data)

    # CLR normalization on ADT
    adt_assay <- CreateAssayObject(counts = adt)
    multi[["ADT"]] <- adt_assay

    DefaultAssay(multi) <- "ADT"
    multi <- NormalizeData(multi, normalization.method = 'CLR', margin = 2)

    # prepare input of bindSC

    ## X: adt matrix from CITE-seq data
    X <- multi@assays$ADT@data

    ## Y: gene expression matrix from RNA data
    Y <- rna_normalize[intersect(multi_hvgs, rownames(rna_normalize)),]
    # rownames(X) <- str_split(rownames(X), pattern = "\\.", simplify = T)[,1]

    # Z: initilized gene activity matrix from RNA data
    # gene_adt_tab  <- read.table(config["gene_adt_table"], header = 1)
    gene_adt_list <- intersect(rownames(rna), rownames(adt))
    score.matrix <- as.data.frame(multi@assays$SCT@data)
    Z <- na.omit(score.matrix[gene_adt_list, ])
    Z <- as.matrix(Z)

    # prepare input for BiCCA
    shared_adt <- intersect(rownames(X), rownames(Z))
    y <- Y
    x <- X[shared_adt, ]
    z0 <- Z[shared_adt, ]


    bindsc_res <- BiCCA(X = x ,
                 Y = y, 
                 Z0 =z0, 
                 X.clst = rep(0, dim(x)[2]),
                 Y.clst = rep(0, dim(y)[2]),
                 alpha = 0.1, 
                 lambda = 0.7,
                 K = min(15, length(shared_adt)-1),
                 temp.path  = "out",
                 num.iteration = 50,
                 tolerance = 0.01,
                 save = TRUE,
                 parameter.optimize = FALSE, 
                 block.size = 0)

    # Visulize
    library(umap)
    latent <- rbind(bindsc_res$u, bindsc_res$r)
    umap_plt <- umap(latent)
    umap_plt2  <- data.frame("UMAP1"=umap_plt$layout[, 1],
                             "UMAP2"=umap_plt$layout[, 2],
                             "data" = c(rep("RNA", nrow(bindsc_res$u)),
                                        rep("ADT", nrow(bindsc_res$r))))
    umap_plt2$barcodes <- rownames(umap_plt$layout)

    # # Cluster
    # DRdist <- dist(umap_plt2)
    # dclust <- densityClust(DRdist, gaussian=T)
    # dclust <- findClusters(dclust, rho = 20, delta = 2.5)
    # densityClusts <- dclust$clusters
    # umap_plt2$cluster <- densityClusts
    # # densityClusts <- as.data.frame(densityClusts)

    # Save Results
    ## save latent
    write.table(bindsc_res$u,
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-RNA-latent_1.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(bindsc_res$r,
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-ADT-latent_2.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    latent <- merge(bindsc_res$u, bindsc_res$r, by = "row.names", all = TRUE)
    rownames(latent) <- latent[, 1]
    latent <- as.matrix(latent[-1])
    write.table(latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-multi-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## save UMAP
    rna_tab <- subset(umap_plt2, data=="RNA")[,c(1,2,4)]
    rownames(rna_tab) <- rna_tab$barcodes
    adt_tab <- subset(umap_plt2, data=="ADT")[,c(1,2,4)]
    rownames(adt_tab) <- adt_tab$barcodes

    write.table(rna_tab, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-RNA-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(adt_tab, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-ADT-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    # # impute atac promoter
    # Z_impu <- impuZ(X=X, bicca = bindsc_res)
    # write.table(Z_impu,
    #             file = here(output_path, paste(config["output_prefix"], 
    #                                            "bindSC-ATAC-imputation.csv", sep = "-")), 
    #             row.names =T, col.names = T, sep=',', quote=F)
}


# %%
bindSC_module(input_path = args[1],
              output_path = args[2],
              config = unlist(fromJSON(file=args[3]))
             )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData"
# output_path <- "/home/wsg/BM/results/task/scRNA+ADT/accuracy/10x_NSCLC/RawData/rep_1/run_bindSC"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData/RawData.json"))

# %%
# bindSC_module(input_path,
#               output_path,
#               config)
