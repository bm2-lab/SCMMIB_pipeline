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
input_path <- "/home/wsg/BM/data/HWBC/RawData"
output_path <- "/home/wsg/BM/data/HWBC/RNA+ADT/RawData"

# %%
HWBC_multi <- LoadH5Seurat(file = paste0(input_path, "/multi.h5seurat"))

# %%
metadata <- HWBC_multi@meta.data

# %%
mat <- HWBC_multi@assays$SCT@counts

# %%
mat

# %%
# metadata
metadata$batch <- paste0(metadata$donor, "_", metadata$time)
metadata$barcode <- rownames(metadata)
# metadata
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
RNA_subset_counts <- mat

# %%
process = "raw"

# save raw rna to mtx
data_path <- here(output_path, "HWBC-CITE_seq-raw-RNA-counts.mtx")
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, "HWBC-CITE_seq-raw-RNA-counts.rds"))

# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, "HWBC-CITE_seq-raw-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "HWBC-CITE_seq-raw-RNA-counts.h5Seurat"), dest = "h5ad")

# %%
# ADT

# %%
mat <- HWBC_multi@assays$ADT@counts

# %%
mat

# %%
ADT_subset_counts <- mat

# %%
# save raw ADT to mtx
data_path <- here(output_path, "HWBC-CITE_seq-raw-ADT-counts.mtx")
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")

# save raw ADT to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, "HWBC-CITE_seq-raw-ADT-counts.rds"))

# Create Seurat Object
ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, "HWBC-CITE_seq-raw-ADT-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "HWBC-CITE_seq-raw-ADT-counts.h5Seurat"), dest = "h5ad")

# %%
