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
library(tidyr)
library(dplyr)
library(DropletUtils) 

library(Seurat)
library(Signac)
library(SeuratDisk)

library(rhdf5)
library(anndata)

# %% [markdown]
# # convert H5Seurat

# %% [markdown]
# ## metadata

# %%
input_path <- "/home/wsg/BM/data/10x_mouse_brain/RawData"
output_path <- "/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5"

# %%
brain_multi <- Read10X_h5(file = here(input_path, "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5"))

# %%
column_names <- colnames(RNA_counts)
column_names_split <- strsplit(column_names, "-")

# %%
metadata <- data.frame(barcode = column_names, 
                       index = sapply(column_names_split, tail, 1))
match_list <- c('1' = 'AD_17p9_rep4', 
                '2' = 'AD_17p9_rep5', 
                '3' = 'AD_2p5_rep2',
                '4' = 'AD_2p5_rep3', 
                '5' = 'AD_5p7_rep2', 
                '6' = 'AD_5p7_rep6',
                '7' = 'WT_13p4_rep2', 
                '8' = 'WT_13p4_rep5', 
                '9' = 'WT_2p5_rep2',
                '10' = 'WT_2p5_rep7', 
                '11' = 'WT_5p7_rep2', 
                '12' = 'WT_5p7_rep3')
metadata$sample <- match_list[as.character(metadata$index)]

# %%
metadata$separate <- metadata$sample
metadata <- metadata %>%
    separate(col = separate, into = c("type", "time", "rep"), sep = "_")
metadata

# %%
metadata_2p5 <- metadata[which(metadata$time == "2p5"),]
rownames(metadata_2p5) <- metadata_2p5$barcode
metadata_2p5

# %%
# metadata
metadata_2p5
write_csv(metadata_2p5, here(output_path, "metadata.csv"))

# %% [markdown]
# ## RNA

# %%
RNA_counts <- brain_multi$`Gene Expression`

# %%
RNA_counts

# %%
sum(metadata_2p5$barcode %in% colnames(RNA_counts))

# %%
RNA_subset_counts <- RNA_counts[, metadata_2p5$barcode]

# %%
RNA_subset_counts

# %%
process = "2p5"

# save 2p5 rna to mtx
data_path <- here(output_path, "brain-multiome-2p5-RNA-counts.mtx")
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# save 2p5 rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, "brain-multiome-2p5-RNA-counts.rds"))

# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata_2p5)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, "brain-multiome-2p5-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "brain-multiome-2p5-RNA-counts.h5Seurat"), dest = "h5ad")

# %% [markdown]
# ## ATAC

# %%
ATAC_counts <- brain_multi$Peaks

# %%
ATAC_counts_rename <- ATAC_counts
rownames(ATAC_counts_rename) <- gsub(rownames(ATAC_counts_rename), pattern = ":", replacement = "-")

# %%
ATAC_subset_counts <- ATAC_counts_rename[, metadata_2p5$barcode]

# %%
rows_to_keep <- !grepl("^GRCh38_chr21", rownames(ATAC_subset_counts))
ATAC_subset_counts <- ATAC_subset_counts[rows_to_keep, ]

# %%
ATAC_subset_counts
# tail(ATAC_subset_counts)

# %%
# save 2p5 ATAC to mtx
data_path <- here(output_path, "brain-multiome-2p5-ATAC-peaks.mtx")
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")

# save 2p5 ATAC to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, "brain-multiome-2p5-ATAC-peaks.rds"))

# Create Seurat Object
ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, meta.data = metadata_2p5)

# save Seurat to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, "brain-multiome-2p5-ATAC-peaks.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "brain-multiome-2p5-ATAC-peaks.h5Seurat"), dest = "h5ad")

# %%
