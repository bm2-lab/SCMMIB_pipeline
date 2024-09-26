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
input_path <- "/home/wsg/BM/data/10x_NSCLC/RawData"
output_path <- "/home/wsg/BM/data/10x_NSCLC/RNA+ADT/RawData"

# %%
NSCLC_multi <- Read10X_h5(file = here(input_path, "20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5"))

# %%
NSCLC_table <- read.csv(file = here(input_path, "20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_multiplexing_analysis_assignment_confidence_table.csv"))
table(NSCLC_table$Assignment)

# %%
metadata <- NSCLC_table[,c("Barcodes", "Assignment", "Assignment_Probability")]
colnames(metadata) <- c("barcode", "donor", "donor_probability")
metadata <- metadata[-which(metadata$donor %in% c("Blanks", "Multiplet", "Unassigned")),]

# %%
# metadata
metadata
write_csv(metadata, here(output_path, "metadata.csv"))

# %%
RNA_counts <- NSCLC_multi$`Gene Expression`

# %%
RNA_counts

# %%
sum(metadata$barcode %in% colnames(RNA_counts))

# %%
RNA_subset_counts <- RNA_counts[, metadata$barcode]

# %%
process = "raw"

# save raw rna to mtx
data_path <- here(output_path, "NSCLC-CITE_seq-raw-RNA-counts.mtx")
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, "NSCLC-CITE_seq-raw-RNA-counts.rds"))

# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, "NSCLC-CITE_seq-raw-RNA-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "NSCLC-CITE_seq-raw-RNA-counts.h5Seurat"), dest = "h5ad")

# %%
# ADT

# %%
ADT_counts <- NSCLC_multi$`Antibody Capture`

# %%
ADT_counts

# %%
ADT_subset_counts <- ADT_counts[, metadata$barcode]

# %%
ADT_subset_counts

# %%
gene_adt_tab <- read.table("/home/wsg/BM/pipeline/config/gene_adt_tab.tsv", header = 1)

# %%
adt_raw <- rownames(ADT_subset_counts)
adt_raw

# %%
adt_new <- sapply(strsplit(adt_raw, "\\."), '[', 1)
adt_new

# %%
gene_adt_tab_sub <- gene_adt_tab[which(gene_adt_tab$ADT %in% adt_new), ]
gene_adt_tab_sub

# %%
adt_new

# %%
adt_gene <- c('CD3', 'CD4', 'CD8A', 'ITGAX', 'CD14', 'FCGR3A', 'CD19', 'NCAM1', 'PTPRC')
adt_gene

# %%
rownames(ADT_subset_counts) <- adt_gene
ADT_subset_counts

# %%
# save raw ADT to mtx
data_path <- here(output_path, "NSCLC-CITE_seq-raw-ADT-counts.mtx")
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")

# save raw ADT to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, "NSCLC-CITE_seq-raw-ADT-counts.rds"))

# Create Seurat Object
ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, "NSCLC-CITE_seq-raw-ADT-counts.h5Seurat"))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(here(output_path, "NSCLC-CITE_seq-raw-ADT-counts.h5Seurat"), dest = "h5ad")

# %%
