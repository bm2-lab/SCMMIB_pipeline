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

library(CiteFuse)
library(scater)
library(SingleCellExperiment)
library(DT)

# %%
CiteFuse_module <- function(input_path,
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

    # preprocess
    CITEseq_multi=list()
    CITEseq_multi[['RNA']] <- rna
    CITEseq_multi[['ADT']] <- adt

    sce_citeseq <- preprocessing(CITEseq_multi)
    sce_citeseq <- scater::logNormCounts(sce_citeseq)
    sce_citeseq <- normaliseExprs(sce_citeseq, altExp_name = "ADT", transform = "log")
    sce_citeseq <- CiteFuse(sce_citeseq) # Merge Data

    sce_citeseq <- reducedDimSNF(sce_citeseq,
                             method = "UMAP",
                             dimNames = "UMAP_joint") # Reduce Dim
    SNF_W_louvain <- igraphClustering(sce_citeseq, method = "louvain") # Cluster
    sce_citeseq$SNF_W_louvain <- as.factor(SNF_W_louvain)

    out_umap <- data.frame(sce_citeseq@int_colData$reducedDims)
    rownames(out_umap) <- rownames(sce_citeseq@colData)
    colnames(out_umap) <- c("UMAP1","UMAP2")
    out_umap$louvain_cluster <- SNF_W_louvain
    out_snf <- sce_citeseq@metadata$SNF_W

    # Save Results
    ## save latent
    write.table(out_snf,
                file = here(output_path, paste(config["output_prefix"], 
                                               "CiteFuse-multi-graph.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## save UMAP
    write.table(out_umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "CiteFuse-multi-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
    }

# %%
CiteFuse_module(input_path = args[1],
                output_path = args[2],
                config = unlist(fromJSON(file=args[3]))
               )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/COVID19/RNA+ADT/RawData"
# output_path <- "/home/wsg/BM/results/task/scRNA+ADT/accuracy/COVID19/RawData/rep_1/run_CiteFuse"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/COVID19/RNA+ADT/RawData/RawData.json"))

# %%
# CiteFuse_module(input_path,
#                 output_path,
#                 config)

# %%
