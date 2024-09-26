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

library(rhdf5)
library(anndata)

# %% [markdown]
# # convert HSPC mtx to rds

# %% [markdown]
# ## Multiome

# %%
input_path <- "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
output_path <- "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10/p10.json"))

# %%
HSPC_RNA_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-multiome-p10-RNA-counts.h5ad"))

# %%
HSPC_RNA_p10$X

# %%
mat <- t(HSPC_RNA_p10$X)
RNA_subset_counts <- Matrix(mat, sparse=TRUE)

# %%
RNA_subset_counts

# %%
process = "p10"
# # save raw rna to rds
# saveRDS(RNA_subset_counts, 
#         file = here(output_path, "HSPC-multiome-p10-RNA-counts.rds"))

# # save raw rna to mtx
# data_path <- here(output_path, "HSPC-multiome-p10-RNA-counts.mtx")
# write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# metadata
metadata <- HSPC_RNA_p10$obs
metadata['cell_id'] <- metadata['barcode']
metadata <- metadata[,c(6,1,2,3,4,5)]
# metadata
write_csv(metadata, here(output_path, "metadata.csv"))

# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, "HSPC-multiome-p10-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "HSPC-multiome-p10-RNA-counts.h5Seurat"), dest = "h5ad")

# %%
# ATAC

# %%
input_path <- "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
output_path <- "/home/wsg/BM/data/HSPC/RNA+ATAC/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/HSPC/RNA+ATAC/p10/p10.json"))

# %%
HSPC_ATAC_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-multiome-p10-ATAC-peaks.h5ad"))
HSPC_ATAC_p10

# %%
HSPC_ATAC_p10$X

# %%
mat <- t(HSPC_ATAC_p10$X)
ATAC_subset_counts <- Matrix(mat, sparse=TRUE)

# %%
ATAC_subset_counts

# %%
rownames(ATAC_subset_counts) <- gsub(pattern = ":", replacement = "-", x = rownames(ATAC_subset_counts))

# %%
# ATAC_subset_counts

# %%
# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, "HSPC-multiome-p10-ATAC-peaks.rds"))

# save raw atac to mtx
data_path <- here(output_path, "HSPC-multiome-p10-ATAC-peaks.mtx")
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")

# metadata
metadata <- HSPC_ATAC_p10$obs
# metadata['cell_id'] <- metadata['barcode']
# metadata <- metadata[,c(6,1,2,3,4,5)]
# metadata
# write_csv(metadata, here(output_path, "metadata.csv"))

# Create Seurat Object
ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, "HSPC-multiome-p10-ATAC-peaks.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "HSPC-multiome-p10-ATAC-peaks.h5Seurat"), dest = "h5ad")

# %% [markdown]
# ## CITE-seq

# %%
input_path <- "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
output_path <- "/home/wsg/BM/data/HSPC/RNA+ADT/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/HSPC/RNA+ADT/p10/p10.json"))

# %%
HSPC_RNA_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-CITE_seq-raw_p10-RNA-counts.h5ad"))

# %%
HSPC_RNA_p10

# %%
var_names_list = strsplit(HSPC_RNA_p10$var_names, split = "_")
var_names_symbols = sapply(var_names_list, `[`, 2)

# %%
HSPC_RNA_p10$var_names <- var_names_symbols

# %%
mat <- t(HSPC_RNA_p10$X)
RNA_subset_counts <- Matrix(mat, sparse=TRUE)

# %%
RNA_subset_counts[c("RF00017","RF00019"), ]

# %%
mat <- RNA_subset_counts

row_totals <- rowSums(mat)

row_names <- rownames(mat)
duplicated_rows <- row_names[duplicated(row_names) | duplicated(row_names, fromLast = TRUE)]

result_mat <- mat[!row_names %in% duplicated_rows, , drop = FALSE]

for (row_name in unique(duplicated_rows)) {
  duplicate_indices <- which(row_names == row_name)
  max_index <- duplicate_indices[which.max(row_totals[duplicate_indices])]
  result_mat <- rbind(result_mat, mat[max_index, , drop = FALSE])
}

result_mat <- result_mat[order(rownames(result_mat)), , drop = FALSE]

# %%
RNA_subset_counts <- result_mat

# %%
RNA_subset_counts

# %%
process = "p10"

# save raw rna to mtx
data_path <- here(output_path, "HSPC-CITE_seq-p10-RNA-counts.mtx")
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, "HSPC-CITE_seq-p10-RNA-counts.rds"))

# %%
# # metadata
# metadata <- HSPC_RNA_p10$obs
# write_csv(metadata, here(output_path, "metadata.csv"))

# # Create Seurat Object
# h5Seurat <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# # save Seurat to h5Seurat
SaveH5Seurat(h5Seurat, 
             file = here(output_path, "HSPC-CITE_seq-p10-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert("HSPC-CITE_seq-p10-RNA-counts.h5Seurat", dest = "h5ad", overwrite = TRUE)

# %%

# %%

# %%
HSPC_ADT_p10 <- read_h5ad(file = paste0(output_path, "/HSPC-CITE_seq-raw_p10-ADT-counts.h5ad"))
HSPC_ADT_p10

# %%
mat <- t(HSPC_ADT_p10$X)
ADT_subset_counts <- Matrix(mat, sparse=TRUE)

# %%
ADT_subset_counts

# %%
length(HSPC_ADT_p10$var_names)

length(intersect(HSPC_ADT_p10$var_names, HSPC_RNA_p10$var_names))

# %%
process = "p10"

# save raw rna to mtx
data_path <- here(output_path, "HSPC-CITE_seq-p10-ADT-counts.mtx")
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")

saveRDS(ADT_subset_counts, 
        file = here(output_path, "HSPC-CITE_seq-p10-ADT-counts.rds"))

# metadata
metadata <- HSPC_ADT_p10$obs
# # Create Seurat Object
h5Seurat <- CreateSeuratObject(counts = ADT_subset_counts, meta.data = metadata)

# # save Seurat to h5Seurat
SaveH5Seurat(h5Seurat, 
             file = here(output_path, "HSPC-CITE_seq-p10-ADT-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert("HSPC-CITE_seq-p10-ADT-counts.h5Seurat", dest = "h5ad", overwrite = TRUE)

# %%

# %%

# %%

# %%
atac_subset_meta <- as.data.frame(atac_subset_meta)

# subset counts matrix
ATAC_subset_counts <- ATAC_counts[, atac_subset_meta$barcode]
colnames(ATAC_subset_counts) <- make.unique(colnames(ATAC_subset_counts))

atac_subset_meta$barcode <- make.unique(atac_subset_meta$barcode)
rownames(atac_subset_meta) <- atac_subset_meta$barcode

# %%
dim(ATAC_subset_counts)

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(atac_subset_meta, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = atac_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ATAC", "peaks.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ATAC", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
