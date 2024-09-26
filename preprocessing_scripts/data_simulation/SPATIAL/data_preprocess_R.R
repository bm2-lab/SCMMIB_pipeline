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
library(tidyverse)
library(dplyr)
library(DropletUtils) 

library(Seurat)
library(Signac)
library(SeuratDisk)
library(SeuratData)

library(rhdf5)
library(anndata)

# %%
# InstallData("stxBrain")

# %%
# brain <- LoadData("stxBrain", type = "anterior1")

# %% [markdown]
# # convert spatial h5ad to rds

# %% [markdown]
# ## lymph node

# %%
input_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node"
output_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node"
# dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/SPATIAL/RNA_ADT/lymph_node/lymph_node.json"))

# %%
h5ad <- read_h5ad(file = paste0(input_path, "/lymph_node-CITE_seq_RNA-counts.h5ad"))

# %%
h5ad$obs

# %%
metadata <- h5ad$obs
metadata['barcode'] <- rownames(metadata)

# %%
spatial_axis <- as.data.frame(h5ad$obsm$spatial)
colnames(spatial_axis) = c("X", "Y")
metadata <- cbind(metadata, spatial_axis)
head(metadata)

# %%
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
process = "raw"
RNA_counts <- t(h5ad$X)
dim(RNA_counts)

# %%
RNA_counts <- as(RNA_counts, 'CsparseMatrix')
RNA_counts

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste("lymph_node-CITE_seq", process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_counts, 
        file = here(output_path, 
                    paste("lymph_node-CITE_seq", process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste("lymph_node-CITE_seq", process,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# # Convert h5Seurat to h5ad
# setwd(output_path)
# Convert(paste("lymph_node-CITE_seq", process, "RNA", "counts.h5Seurat", sep = "-"), 
#         dest = "h5ad")

# %%

# %%
h5ad <- read_h5ad(file = paste0(input_path, "/lymph_node-CITE_seq-raw-ADT-counts.h5ad"))

# %%
h5ad

# %%
metadata <- h5ad$obs
metadata['barcode'] <- rownames(metadata)

# %%
spatial_axis <- as.data.frame(h5ad$obsm$spatial)
colnames(spatial_axis) = c("X", "Y")
metadata <- cbind(metadata, spatial_axis)
head(metadata)

# %%
process = "raw"
ADT_counts <- t(h5ad$X)
ADT_counts <- as(ADT_counts, "sparseMatrix")
dim(ADT_counts)

# %%
ADT_counts <- as(ADT_counts, 'CsparseMatrix')
ADT_counts

# %%
# ADT_counts
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
# save raw adt to mtx
data_path <- here(output_path,
                 paste("lymph_node-CITE_seq", process,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_counts, path = data_path, version = "3")
# write_csv(metadata, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_counts, 
        file = here(output_path, 
                    paste("lymph_node-CITE_seq", process,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_counts, assay = "ADT", meta.data = metadata)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste("lymph_node-CITE_seq", process,
                                   "ADT", "counts.h5Seurat", sep = "-")))

# # save raw adt to h5ad
# setwd(output_path)
# Convert(paste("lymph_node-CITE_seq", process, "ADT", "counts.h5Seurat", sep = "-"), 
#         dest = "h5ad")

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
process = "p10"
# save raw rna to rds
RNA_subset_counts <- HSPC_RNA_p10$X
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))

# metadata
metadata <- HSPC_RNA_p10$obs
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# %%

# %%
# ATAC

# %%
input_path <- "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
output_path <- "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10/p10.json"))

# %%
h5Seurat@assays

# %%
HSPC_ATAC_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-multiome-p10-ATAC-peaks.h5ad"))
HSPC_ATAC_p10

# %%
process = "p10"
# save raw rna to rds
ATAC_subset_counts <- HSPC_ATAC_p10$X
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ATAC", "peaks.rds", sep = "-")))

# metadata
metadata <- HSPC_ATAC_p10$obs
# Create Seurat Object
ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ATAC", "peaks.h5Seurat", sep = "-")))

# %%

# %% [markdown]
# ## thymus

# %%
input_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/thymus"
output_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/thymus"

# %%
h5ad <- read_h5ad(file = paste0(input_path, "/thymus-CITE_seq-raw-RNA-counts.h5ad"))

# %%
h5ad$obs

# %%
metadata <- h5ad$obs
metadata['barcode'] <- rownames(metadata)

# %%
spatial_axis <- as.data.frame(h5ad$obsm$spatial)
colnames(spatial_axis) = c("X", "Y")
metadata <- cbind(metadata, spatial_axis)
head(metadata)

# %%
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
process = "raw"
RNA_counts <- t(h5ad$X)
dim(RNA_counts)

# %%
RNA_counts <- as(RNA_counts, 'CsparseMatrix')
RNA_counts

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste("thymus-CITE_seq", process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_counts, 
        file = here(output_path, 
                    paste("thymus-CITE_seq", process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste("thymus-CITE_seq", process,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# # Convert h5Seurat to h5ad
# setwd(output_path)
# Convert(paste("thymus-CITE_seq", process, "RNA", "counts.h5Seurat", sep = "-"), 
#         dest = "h5ad")

# %%

# %%
h5ad <- read_h5ad(file = paste0(input_path, "/thymus-CITE_seq-raw-ADT-counts.h5ad"))

# %%
h5ad

# %%
metadata <- h5ad$obs
metadata['barcode'] <- rownames(metadata)

# %%
spatial_axis <- as.data.frame(h5ad$obsm$spatial)
colnames(spatial_axis) = c("X", "Y")
metadata <- cbind(metadata, spatial_axis)
head(metadata)

# %%
process = "raw"
ADT_counts <- t(h5ad$X)
ADT_counts <- as(ADT_counts, "sparseMatrix")
dim(ADT_counts)

# %%
ADT_counts <- as(ADT_counts, 'CsparseMatrix')
ADT_counts

# %%
# save raw adt to mtx
data_path <- here(output_path,
                 paste("thymus-CITE_seq", process,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_counts, path = data_path, version = "3")
# write_csv(metadata, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_counts, 
        file = here(output_path, 
                    paste("thymus-CITE_seq", process,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_counts, assay = "ADT", meta.data = metadata)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste("thymus-CITE_seq", process,
                                   "ADT", "counts.h5Seurat", sep = "-")))

# # save raw adt to h5ad
# setwd(output_path)
# Convert(paste("thymus-CITE_seq", process, "ADT", "counts.h5Seurat", sep = "-"), 
#         dest = "h5ad")

# %%

# %%

# %% [markdown]
# ## spleen

# %%
input_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen"
output_path <- "/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen"

# %%
h5ad <- read_h5ad(file = paste0(input_path, "/spleen-CITE_seq-raw-RNA-counts.h5ad"))

# %%
metadata <- h5ad$obs
metadata['barcode'] <- rownames(metadata)

# %%
spatial_axis <- as.data.frame(h5ad$obsm$spatial)
colnames(spatial_axis) = c("X", "Y")
metadata <- cbind(metadata, spatial_axis)
head(metadata)

# %%
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
process = "raw"
RNA_counts <- t(h5ad$X)
dim(RNA_counts)

# %%
RNA_counts <- as(RNA_counts, 'CsparseMatrix')
RNA_counts

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste("spleen-CITE_seq", process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_counts, 
        file = here(output_path, 
                    paste("spleen-CITE_seq", process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste("spleen-CITE_seq", process,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# # Convert h5Seurat to h5ad
# setwd(output_path)
# Convert(paste("spleen-CITE_seq", process, "RNA", "counts.h5Seurat", sep = "-"), 
#         dest = "h5ad")

# %%

# %%

# %%
h5ad <- read_h5ad(file = paste0(input_path, "/spleen-CITE_seq-raw-ADT-counts.h5ad"))

# %%
metadata <- h5ad$obs
metadata['barcode'] <- rownames(metadata)

# %%
spatial_axis <- as.data.frame(h5ad$obsm$spatial)
colnames(spatial_axis) = c("X", "Y")
metadata <- cbind(metadata, spatial_axis)
head(metadata)

# %%
process = "raw"
ADT_counts <- t(h5ad$X)
ADT_counts <- as(ADT_counts, "sparseMatrix")
dim(ADT_counts)

# %%
ADT_counts <- as(ADT_counts, 'CsparseMatrix')
ADT_counts

# %%
# save raw adt to mtx
data_path <- here(output_path,
                 paste("spleen-CITE_seq", process,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_counts, path = data_path, version = "3")
# write_csv(metadata, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_counts, 
        file = here(output_path, 
                    paste("spleen-CITE_seq", process,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_counts, assay = "ADT", meta.data = metadata)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste("spleen-CITE_seq", process,
                                   "ADT", "counts.h5Seurat", sep = "-")))

# # save raw adt to h5ad
# setwd(output_path)
# Convert(paste("spleen-CITE_seq", process, "ADT", "counts.h5Seurat", sep = "-"), 
#         dest = "h5ad")

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## CITE-seq

# %%
SPOTS_cite_RNA_path = "/home/wsg/BM/data/SPOTS/GEO/GSE198353_spleen_rep_1_filtered_feature_bc_matrix.h5"

# %%
library(rhdf5)

# 打开并读取HDF5文件
file_path <- SPOTS_cite_RNA_path
dataset <- h5read(file = file_path, name = "matrix")


# %%
str(dataset)

# %%
as.matrix(dataset)

# %%
h5ls(SPOTS_cite_RNA_path)

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
input_path <- "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
output_path <- "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/HSPC/RNA+ADT/p10/p10.json"))

# %%

# %%
HSPC_RNA_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-CITE_seq-p10-RNA-counts.h5ad"))

# %%
HSPC_RNA_p10

# %%
process = "p10"
# save raw rna to rds
RNA_subset_counts <- HSPC_RNA_p10$X
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))

# metadata
metadata <- HSPC_RNA_p10$obs
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# %%

# %%
HSPC_ADT_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-CITE_seq-p10-ADT-counts.h5ad"))
HSPC_ADT_p10

# %%
process = "p10"
# save raw rna to rds
ADT_subset_counts <- HSPC_ADT_p10$X
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "counts.rds", sep = "-")))

# metadata
metadata <- HSPC_ADT_p10$obs
# Create Seurat Object
ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "counts.h5Seurat", sep = "-")))

# %%

# %%

# %%

# %%

# %%

# %%
# save raw rna to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"],"raw", dataset["task_type"], "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")


# %%
# Load ATAC
BMMC_ATAC_Dir <- "/Data/wangsg/BM/pipeline/results/BMMC/data_preprocess/BMMC-raw-pair-ATAC-peaks.mtx/"
BMMC_ATAC_counts <- Read10X(data.dir = BMMC_ATAC_Dir, gene.column = 1)

# %%
metadata <- read.csv(paste0(BMMC_ATAC_Dir, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

# set.seed(1234)
# random sample 500 cells of each donor
# bmmc_rna_500_meta <- metadata %>% group_by(batch) %>% slice_sample(n=500)

# random sample 10% cells of each donor
# bmmc_rna_10p_meta <- metadata %>% group_by(batch) %>% sample_frac(.1)

table(metadata$batch)
bmmc_atac_10p_meta <- metadata[rownames(metadata) %in% rownames(bmmc_rna_10p_meta), ]
table(bmmc_atac_10p_meta$batch)

# %%
BMMC_ATAC_counts_10p <- BMMC_ATAC_counts[ , colnames(BMMC_ATAC_counts) %in% bmmc_atac_10p_meta$barcode]

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       "raw",
                       dataset["task_type"], 
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = BMMC_ATAC_counts_10p, path = data_path, version = "3")
write_csv(bmmc_atac_10p_meta, here(data_path, "metadata.csv"))

# %%
# save raw atac to rds
saveRDS(BMMC_ATAC_counts_10p, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          "raw",
                          dataset["task_type"], 
                          "ATAC", "peaks.rds", sep = "-")))

# %%
chrom_assay <- CreateChromatinAssay(
    counts = BMMC_ATAC_counts_10p,
    sep = c("-", "-")
)
bmmc_atac_10p <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", meta.data = bmmc_atac_10p_meta)
