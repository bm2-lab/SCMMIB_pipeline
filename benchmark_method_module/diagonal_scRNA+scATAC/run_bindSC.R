# -*- coding: utf-8 -*-
# %%
args <- commandArgs(T) 

# %%
library(here)
library(rjson)
library(Matrix)

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
    atac <- readRDS(here(input_path, config["atac_rds_filename"]))

    metadata <- read.csv(here(input_path, config["metadata"]), row.names = config["barcode_key"])

    # Get GAM of ATAC
    activity.matrix <- readRDS(here(input_path, config["gam_rds_filename"]))

    # SCTransform on RNA
    wnn_data <- CreateSeuratObject(counts = rna)
    DefaultAssay(wnn_data) <- "RNA"
    wnn_data <- SCTransform(wnn_data, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    wnn_data <- FindVariableFeatures(wnn_data, selection.method = "vst", nfeatures = 3000)     

    wnn_hvgs <-VariableFeatures(wnn_data)
    rna_normalize=as.matrix(wnn_data[["RNA"]]@data)

    # TFIDF on ATAC
    chrom_assay <- CreateChromatinAssay(
       counts = atac,
       sep = c("-", "-")
     )

    wnn_data[["ATAC"]] <- chrom_assay

    DefaultAssay(wnn_data) <- "ATAC"
    wnn_data <- RunTFIDF(wnn_data)
    wnn_data <- FindTopFeatures(wnn_data, min.cutoff = 'q0')
    wnn_data <- RunSVD(wnn_data)
    wnn_data <- RunUMAP(wnn_data, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

    # prepare input of bindSC

    ## X: gene expression matrix from scRNA-seq data
    X <- rna_normalize[intersect(wnn_hvgs,rownames(rna_normalize)),]
    rownames(X) <- str_split(rownames(X),pattern = "\\.",simplify = T)[,1]

    ## Y: chromatin accessibility matrix from scATAC-seq data
    y <- wnn_data[["lsi"]]@cell.embeddings

    # Z0: initilized gene activity matrix from scATAC-seq data
    Z0 <- as.matrix(activity.matrix)
    Z0 <- Z0[intersect(rownames(Z0), rownames(X)),]

    X <- X[intersect(rownames(Z0), rownames(X)),]

    # dimension reductions (SVD)
    ## Since the feature dimension is high for single cell epigenetic profiles. 
    ## Here we use low-dimension representations rather than the orignal matrixs for bindSC alignment. 
    ## You should perform dimension reductions (SVD) on X and Z0 jointly.
    out <- dimReduce( dt1 =  X, dt2 = Z0,  K = 50)
    x <- out$dt1
    z0 <- out$dt2

    # Run bindSC
    K <- 5
    bindsc_res <- BiCCA(X=t(x),
               Z0=t(z0), 
               Y=t(y), 
               alpha = 0.1,
               lambda = 0.5,
               K = K, 
               #X.clst = X.clst, # 已知标签或者聚类。
               X.clst = rep(0,dim(x)[1]),
               #Y.clst = Y.clst,
               Y.clst = rep(0,dim(y)[1]),
               num.iteration = 50, 
               temp.path = output_path,
               tolerance = 0.01, 
               save = TRUE, 
               block.size = 0)

#     # Visulize
#     library(umap)
#     latent <- rbind(bindsc_res$u, bindsc_res$r)
#     umap_plt <- umap(latent)
#     umap_plt2  <- data.frame("UMAP1"=umap_plt$layout[,1],
#                             "UMAP2"=umap_plt$layout[,2],
#                             "data" = c(rep("scRNA-seq",nrow(x)),
#                                        rep("scATAC-seq",nrow(y))))

#     # Cluster
#     DRdist <- dist(umap_plt2)
#     dclust <- densityClust(DRdist,gaussian=T)
#     dclust <- findClusters(dclust, rho = 20, delta = 2.5)
#     densityClusts <- dclust$clusters
#     umap_plt2$cluster <- densityClusts
#     #densityClusts <- as.data.frame(densityClusts)

    # Save Results
    ## save latent
    write.table(bindsc_res$u,
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-RNA-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    write.table(bindsc_res$r,
                file = here(output_path, paste(config["output_prefix"], 
                                               "bindSC-ATAC-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

#     ## save UMAP
#     rna_tab <- subset(umap_plt2, data=="scRNA-seq")[,c(1,2,4)]
#     rownames(rna_tab) <- colnames(rna)
#     write.table(rna_tab, 
#                 file = here(output_path, paste(config["output_prefix"], 
#                                                "bindSC-RNA-umap.csv", sep = "-")), 
#                 row.names =T, col.names = T, sep=',', quote=F)

#     atac_tab <- subset(umap_plt2, data=="scATAC-seq")[,c(1,2,4)]
#     rownames(atac_tab) <- colnames(atac)
#     write.table(atac_tab, 
#                 file = here(output_path, paste(config["output_prefix"], 
#                                                "bindSC-ATAC-umap.csv", sep = "-")), 
#                 row.names =T, col.names = T, sep=',', quote=F)

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
# input_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/R75_A100"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/robustness/SHARE/R75_A100/run_bindSC"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/SHARE/RNA+ATAC/R75_A100/RawData.json"))

# %%
# bindSC_module(input_path,
#               output_path,
#               config)

# %%
