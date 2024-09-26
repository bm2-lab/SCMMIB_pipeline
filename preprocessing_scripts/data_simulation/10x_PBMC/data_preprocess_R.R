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
# # convert H5Seurat

# %% [markdown]
# ## Multiome

# %%
input_path_10k <- "/home/wsg/BM/data/10x_PBMC/RawData/PBMC_10k"
input_path_3k <- "/home/wsg/BM/data/10x_PBMC/RawData/PBMC_3k"
output_path <- "/home/wsg/BM/data/10x_PBMC/RNA+ATAC/RawData"

# %%
pbmc_10k <- Read10X_h5(file = here(input_path_10k, "pbmc_unsorted_10k_filtered_feature_bc_matrix.h5"))

# %%
pbmc_3k <- Read10X_h5(file = here(input_path_3k, "pbmc_unsorted_3k_filtered_feature_bc_matrix.h5"))

# %%

# %% [markdown]
# ## metadata

# %%
metadata_10k <- data.frame(barcode = paste0(colnames(pbmc_10k$`Gene Expression`), "-10k"),
                           batch = "pbmc_10k")
metadata_3k <- data.frame(barcode = paste0(colnames(pbmc_3k$`Gene Expression`), "-3k"),
                          batch = "pbmc_3k")
metadata <- rbind(metadata_10k, metadata_3k)
rownames(metadata) <- metadata$barcode
# metadata

# %%
table(metadata$batch)
write_csv(metadata, here(output_path, "metadata.csv"))

# %%

# %% [markdown]
# ## RNA

# %%
RNA_counts_10k <- pbmc_10k$`Gene Expression`
RNA_counts_3k <- pbmc_3k$`Gene Expression`

# %%
colnames(RNA_counts_10k) <- paste0(colnames(RNA_counts_10k), "-10k")
colnames(RNA_counts_3k) <- paste0(colnames(RNA_counts_3k), "-3k")

# %%
if(all(rownames(RNA_counts_10k) == rownames(RNA_counts_3k))){
    RNA_counts <- cbind(RNA_counts_10k, RNA_counts_3k)
}

# %%
RNA_counts

# %%
RNA_subset_counts <- RNA_counts

# %%
process = "raw"

# save raw rna to mtx
data_path <- here(output_path, "PBMC-multiome-raw-RNA-counts.mtx")
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, "PBMC-multiome-raw-RNA-counts.rds"))

# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, "PBMC-multiome-raw-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "PBMC-multiome-raw-RNA-counts.h5Seurat"), dest = "h5ad")

# %%

# %% [markdown]
# ## ATAC

# %%
ATAC_counts_10k <- pbmc_10k$Peaks
ATAC_counts_3k <- pbmc_3k$Peaks

# %%
colnames(ATAC_counts_10k) <- paste0(colnames(ATAC_counts_10k), "-10k")
colnames(ATAC_counts_3k) <- paste0(colnames(ATAC_counts_3k), "-3k")

# %%
length(intersect(rownames(ATAC_counts_10k), rownames(ATAC_counts_3k)))

# %%
# remap the peaks
library(rtracklayer)

merged_peaks_path <- paste0(output_path, "/merged.bed")
merged_peaks <- import(merged_peaks_path, format = "BED")
merged_peaks

# %%
library(GenomicRanges)
library(Matrix)

# 将两个原始的peaks矩阵的行名转换为GRanges对象
gr_10k <- GRanges(rownames(ATAC_counts_10k))
gr_3k <- GRanges(rownames(ATAC_counts_3k))

# 计算每个原始peak与合并后peaks的重叠情况
ov_10k <- findOverlaps(gr_10k, merged_peaks)
ov_3k <- findOverlaps(gr_3k, merged_peaks)


# %%
merged_peaks_name <- paste0(merged_peaks@seqnames, "-", start(merged_peaks), "-", end(merged_peaks))

# %%
head(ov_10k)
head(queryHits(ov_10k))
head(subjectHits(ov_10k))

rownames(ATAC_counts_10k)[head(queryHits(ov_10k))]
merged_peaks_name[head(subjectHits(ov_10k))]

# %%
ATAC_counts_10k_rename <- ATAC_counts_10k
rownames(ATAC_counts_10k_rename) <- merged_peaks_name[subjectHits(ov_10k)]


# %%
ATAC_counts_3k_rename <- ATAC_counts_3k
rownames(ATAC_counts_3k_rename) <- merged_peaks_name[subjectHits(ov_3k)]

# %%
length(intersect(rownames(ATAC_counts_10k_rename), rownames(ATAC_counts_3k_rename)))

# %%
ATAC_counts <- SeuratObject::RowMergeSparseMatrices(ATAC_counts_10k_rename, ATAC_counts_3k_rename)

# %%
ATAC_subset_counts <- ATAC_counts

# %%
# save raw ATAC to mtx
data_path <- here(output_path, "PBMC-multiome-raw-ATAC-peaks.mtx")
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")

# save raw ATAC to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, "PBMC-multiome-raw-ATAC-peaks.rds"))

# Create Seurat Object
ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, "PBMC-multiome-raw-ATAC-peaks.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "PBMC-multiome-raw-ATAC-peaks.h5Seurat"), dest = "h5ad")

# %%
