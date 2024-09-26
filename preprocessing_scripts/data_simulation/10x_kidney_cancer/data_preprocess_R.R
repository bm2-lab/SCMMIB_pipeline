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
input_path <- "/home/wsg/BM/data/10x_kidney_cancer/RawData"
output_path <- "/home/wsg/BM/data/10x_kidney_cancer/RNA+ADT/RawData"

# %%
kidney_multi <- Read10X_h5(file = here(input_path, "4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex_count_raw_feature_bc_matrix.h5"))

# %%
barcodes_json <- fromJSON(file = here(input_path, "4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Multiplex_multiplexing_analysis_cells_per_tag.json"))
print(length(barcodes_json$BC001))
print(length(barcodes_json$BC002))
sapply(barcodes_json, length)

# %%
rep1 <- data.frame(barcode = barcodes_json$BC001,
                   batch = rep("kidney_rep1", length(barcodes_json$BC001)),
                   probe_tag = rep("BC001", length(barcodes_json$BC001)))

rep2 <- data.frame(barcode = barcodes_json$BC002,
                   batch = rep("kidney_rep2", length(barcodes_json$BC002)),
                   probe_tag = rep("BC002", length(barcodes_json$BC002)))

# %%
metadata <- rbind(rep1, rep2)

# %%
table(metadata$batch)

# %%
# metadata
metadata
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
RNA_counts <- kidney_multi$`Gene Expression`

# %%
RNA_counts

# %%
sum(metadata$barcode %in% colnames(RNA_counts))

# %%
RNA_subset_counts <- RNA_counts[, metadata$barcode]

# %%
process = "raw"

# save raw rna to mtx
data_path <- here(output_path, "kidney-CITE_seq-raw-RNA-counts.mtx")
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, "kidney-CITE_seq-raw-RNA-counts.rds"))

# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, "kidney-CITE_seq-raw-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "kidney-CITE_seq-raw-RNA-counts.h5Seurat"), dest = "h5ad")

# %%
# ADT

# %%
ADT_counts <- kidney_multi$`Antibody Capture`

# %%
ADT_counts

# %%
ADT_subset_counts <- ADT_counts[, metadata$barcode]
rownames(ADT_subset_counts) <- sapply(strsplit(rownames(ADT_subset_counts),"_"), tail, 1)

# %%
ADT_subset_counts

# %%
# save raw ADT to mtx
data_path <- here(output_path, "kidney-CITE_seq-raw-ADT-counts.mtx")
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")

# save raw ADT to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, "kidney-CITE_seq-raw-ADT-counts.rds"))

# Create Seurat Object
ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, "kidney-CITE_seq-raw-ADT-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "kidney-CITE_seq-raw-ADT-counts.h5Seurat"), dest = "h5ad")

# %%
