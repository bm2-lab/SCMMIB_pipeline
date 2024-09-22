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


# %%
#' @importFrom methods slot<-
#'
#' @rdname optimizeALS
#' @export
#' @method optimizeALS liger
#'
optimizeALS_no_filter <- function(
  object,
  k,
  lambda = 5.0,
  thresh = 1e-6,
  max.iters = 30,
  nrep = 1,
  H.init = NULL,
  W.init = NULL,
  V.init = NULL,
  use.unshared = FALSE,
  rand.seed = 1,
  print.obj = FALSE,
  verbose = TRUE,
  ...
) {
  
  if (use.unshared == FALSE){
#     object <- removeMissingObs(
#     object = object,
#     slot.use = 'scale.data',
#     use.cols = FALSE,
#     verbose = TRUE
#   )
    out <- optimizeALS(
      object = object@scale.data,
      k = k,
      lambda = lambda,
      thresh = thresh,
      max.iters = max.iters,
      nrep = nrep,
      H.init = H.init,
      W.init = W.init,
      V.init = V.init,
      use.unshared = FALSE,
      rand.seed = rand.seed,
      print.obj = print.obj,
      verbose = verbose
    )
    names(x = out$H) <- names(x = out$V) <- names(x = object@raw.data)
    for (i in 1:length(x = object@scale.data)) {
      rownames(x = out$H[[i]]) <- rownames(x = object@scale.data[[i]])
    }
    colnames(x = out$W) <- object@var.genes
    for (i in names(x = out)) {
      slot(object = object, name = i) <- out[[i]]
    }
    object@parameters$lambda <- lambda
    return(object)
  }
  if(use.unshared == TRUE){
    object <- optimize_UANLS(object = object,
                             k = k,
                             lambda = lambda,
                             thresh = thresh,
                             max.iters = max.iters,
                             nrep = nrep,
                             rand.seed = rand.seed,
                             print.obj = print.obj)
  }
}

# %%
LIGER_INMF_module <- function(input_path,
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

    # 1. Create a LIGER object and normalize the shared data
    liger <- createLiger(list(atac = activity.matrix, rna = rna))
    liger <- normalize(liger)
    # 2. use the RNA dataset to select variable shared features
    liger <- selectGenes(liger, datasets.use =2)
    # liger <- selectGenes(liger, var.thresh = 0.1, datasets.use =1)
    # 3. Scale
    liger <- scaleNotCenter(liger)

    # Joint Matrix Factorization
    ## 1. To factorize the datasets and include the unshared datasets, set the use.unshared parameter to TRUE.
    liger <- optimizeALS_no_filter(liger, k=20)

    # after optimizeALS, liger@scale.data$rna and liger@scale.data$atac will remove some zero expression cells
    # make sure liger@cell.data has the same dims of liger@scale.data
    filtered_cells_barcode <- c(rownames(liger@scale.data[[1]]), rownames(liger@scale.data[[2]]))
    liger@cell.data <- liger@cell.data[filtered_cells_barcode, ]

    # Quantile Normalization and Joint Clustering
    liger <- quantile_norm(liger)
    liger <- louvainCluster(liger, resolution = 0.2)

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
                                               "LIGER_INMF-RNA-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    umap_ATAC <- umap[which(sub(".*_", "", rownames(umap)) == 'ATAC'), ]
    rownames(umap_ATAC) <- sub("_.*", "", rownames(umap_ATAC))
    colnames(umap_ATAC) <- c("UMAP1", "UMAP2", "cluster")
    write.table(umap_ATAC, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_INMF-ATAC-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## latent
    rna_latent <- as.data.frame(liger@H$rna)
    rownames(rna_latent) <- sub("_.*", "", rownames(rna_latent))
    write.table(rna_latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_INMF-RNA-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    atac_latent <- as.data.frame(liger@H$atac)
    rownames(atac_latent) <- sub("_.*", "", rownames(atac_latent))
    write.table(atac_latent,
                file = here(output_path, paste(config["output_prefix"], 
                                               "LIGER_INMF-ATAC-latent.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    # write.table(liger@H.norm,
    #             file = here(output_path, paste(config["output_prefix"], 
    #                                            "LIGER_INMF-multi-latent.csv", sep = "-")), 
    #             row.names =T, col.names = T, sep=',', quote=F)
    
}

# %%
LIGER_INMF_module(input_path = args[1],
                  output_path = args[2],
                  config = unlist(fromJSON(file=args[3]))
           )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/test"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_LIGER_INMF"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/test/c1k.json"))

# %%
# LIGER_INMF_module(input_path,
#                   output_path,
#                   config
#                  )
