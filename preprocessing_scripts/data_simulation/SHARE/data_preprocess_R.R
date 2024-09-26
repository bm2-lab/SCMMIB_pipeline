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

# %% [markdown]
# # prepare mtx, rds, h5ad and fragments

# %% [markdown]
# ## RNA

# %%
input_path <- "/home/wsg/BM/data/SHARE/GEO"
output_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/RawData.json"))

# %%
RNA_PATH <- paste0(input_path, "/GSM4156608_skin.late.anagen.rna.counts.txt.gz")
META_PATH <- paste0(input_path, "/GSM4156597_skin_celltype.txt")

# %%
rna_counts = read.table(RNA_PATH, sep="\t", stringsAsFactors = FALSE, header = 1, row.names = 1)
metadata = read.table(META_PATH, sep="\t", stringsAsFactors = FALSE, header = 1, row.names = 1)

# %%
length(which(colnames(rna_counts) %in% metadata$rna.bc))

# %%
rna_counts_filtered <- rna_counts[, metadata$rna.bc]
rna_counts_mat <- as.matrix(rna_counts_filtered)

# %%
rna_counts_sparse_mat <- as(rna_counts_mat, "sparseMatrix")

# %%
rna_counts_sparse_mat

# %%
RNA_subset_counts = rna_counts_sparse_mat
rna_subset_meta = metadata
sample = "raw"

# %%
rna_subset_meta

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"],
                       dataset["task_type"],
                       sample,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %% [markdown]
# ## ATAC

# %%
input_path <- "/home/wsg/BM/data/SHARE/GEO"
META_PATH <- paste0(input_path, "/GSM4156597_skin_celltype.txt")
ATAC_COUNTS_PATH <- paste0(input_path, "/GSM4156597_skin.late.anagen.counts.txt.gz")
PEAKS_PATH <- paste0(input_path, "/GSM4156597_skin.late.anagen.peaks.bed.gz")
BARCODES_PATH <- paste0(input_path, "/GSM4156597_skin.late.anagen.barcodes.txt.gz")
FRAGMENTS_PATH <- paste0(input_path, "/GSM4156597_skin.late.anagen.atac.fragments.bed.gz")

# %%
metadata = read.table(META_PATH, sep="\t", stringsAsFactors = FALSE, header = 1, row.names = 1)

# %%
atac_counts = readMM(ATAC_COUNTS_PATH)

# %%
atac_counts

# %%
peaks = read.table(PEAKS_PATH, sep="\t")
peaks["features"] = paste(peaks$V1, peaks$V2, peaks$V3, sep = "-")

# %%
barcodes = read.table(BARCODES_PATH, sep="\t")

# %%
dimnames(atac_counts) <- list(peaks[,"features"], barcodes$V1)

# %%
atac_counts

# %%
head(colnames(rna_counts_filtered))
head(barcodes$V1)
head(metadata$rna.bc)

tail(colnames(rna_counts_filtered))
tail(barcodes$V1)
tail(metadata$rna.bc)

# %%
colnames(metadata) <- c("barcode", "cell_type")
rownames(metadata) <- metadata$barcode
head(metadata)

# %%
dim(atac_counts)

# %%
ATAC_subset_counts = atac_counts
rna_subset_meta = metadata
sample = "raw"

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       sample,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = rna_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample, 
                                   "ATAC", "peaks.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"ATAC", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%

# %% [markdown]
# # Data Manipulation: Load and Subset BMMC

# %%
input_path <- "/Data/wangsg/BM/pipeline/results/BMMC/data_preprocess"
output_path <- "/Data/wangsg/BM/pipeline/results/BMMC/data_preprocess/pair_10p"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/results/BMMC/data_preprocess/BMMC.json"))

# %%
# Load RNA
BMMC_RNA_Dir <- "/Data/wangsg/BM/pipeline/results/BMMC/data_preprocess/BMMC-raw-pair-RNA-counts.mtx"
BMMC_RNA_counts <- Read10X(data.dir = BMMC_RNA_Dir, gene.column = 1)

# %%
metadata <- read.csv(paste0(BMMC_RNA_Dir, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

set.seed(1234)
# random sample 500 cells of each donor
# bmmc_rna_500_meta <- metadata %>% group_by(batch) %>% slice_sample(n=500)

# random sample 10% cells of each donor
bmmc_rna_10p_meta <- metadata %>% group_by(batch) %>% sample_frac(.1)
table(bmmc_rna_10p_meta$batch)
bmmc_rna_10p_meta <- as.data.frame(bmmc_rna_10p_meta)
rownames(bmmc_rna_10p_meta) <- bmmc_rna_10p_meta$barcode

# %%
BMMC_RNA_counts_10p <- BMMC_RNA_counts[ , colnames(BMMC_RNA_counts) %in% bmmc_rna_10p_meta$barcode]

# %%
# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       "raw",
                       dataset["task_type"], 
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = BMMC_RNA_counts_10p, path = data_path, version = "3")
write_csv(bmmc_rna_10p_meta, here(data_path, "metadata.csv"))

# %%
# save raw rna to rds
saveRDS(BMMC_RNA_counts_10p, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          "raw",
                          dataset["task_type"], 
                          "RNA", "counts.rds", sep = "-")))

# %%
bmmc_rna_10p <- CreateSeuratObject(counts = BMMC_RNA_counts_10p, meta.data = bmmc_rna_10p_meta)

# %%
# save raw rna to h5Seurat
SaveH5Seurat(bmmc_rna_10p, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   "raw", 
                                   dataset["task_type"], 
                                   "RNA", "counts.h5Seurat", sep = "-")))

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

# %%
# save raw rna to h5Seurat
SaveH5Seurat(bmmc_atac_10p, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   "raw", 
                                   dataset["task_type"], 
                                   "ATAC", "peaks.h5Seurat", sep = "-")))

# %%
# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"],"raw", dataset["task_type"], "ATAC", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")


# %% [markdown]
# ## subset cells from single sample

# %%
input_path <- "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/RawData"
output_path <- "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/s3d10"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/dataset.json"))

# %%
# 设置提取样本
# s1d1  s1d2  s1d3  s2d1  s2d4  s2d5 s3d10  s3d3  s3d6  s3d7  s4d1  s4d8  s4d9 
# 6224  6740  4279  4220  6111  4895  6781  4325  1679  1771  8023  9876  4325 
sample = "s3d10"

# %%
paste0(input_path, "/BMMC-multiome-raw-RNA-counts.mtx")

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-multiome-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(RNA_Dir, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

set.seed(1234)


# %%
RNA_counts

# %%
rna_subset_meta <- metadata[metadata$batch==sample, ]
table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)

rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode


# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"],
                       dataset["task_type"],
                       sample,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load binarized ATAC
ATAC_Dir <- paste0(input_path, "/BMMC-multiome-binarized-ATAC-peaks.mtx")
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

metadata <- read.csv(paste0(ATAC_Dir, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

atac_subset_meta <- metadata[metadata$batch==sample, ]
table(atac_subset_meta$batch)
table(atac_subset_meta$cell_type)

atac_subset_meta <- as.data.frame(atac_subset_meta)
rownames(atac_subset_meta) <- atac_subset_meta$barcode


# %%
# subset atac counts matrix
ATAC_subset_counts <- ATAC_counts[ , colnames(ATAC_counts) %in% atac_subset_meta$barcode]

# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       sample,
                       "ATAC", "binarized_peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(atac_subset_meta, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "ATAC", "binarized_peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = atac_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample, 
                                   "ATAC", "binarized_peaks.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"ATAC", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ATAC
ATAC_Dir <- paste0(input_path, "/BMMC-multiome-raw-ATAC-peaks.mtx")
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

metadata <- read.csv(paste0(ATAC_Dir, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

atac_subset_meta <- metadata[metadata$batch==sample, ]
table(atac_subset_meta$batch)
table(atac_subset_meta$cell_type)

atac_subset_meta <- as.data.frame(atac_subset_meta)
rownames(atac_subset_meta) <- atac_subset_meta$barcode

# %%
# subset atac counts matrix
ATAC_subset_counts <- ATAC_counts[ , colnames(ATAC_counts) %in% atac_subset_meta$barcode]

# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       sample,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(atac_subset_meta, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = atac_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample, 
                                   "ATAC", "peaks.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"ATAC", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %% [markdown]
# ## sample N cells

# %%
input_path <- "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/RawData"
output_path <- "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/c50k"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/dataset.json"))

# 设置提取数量
Ncell = 50000
process = "c50k"

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-multiome-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)


# %%
# random sample 500 cells of each donor
# bmmc_rna_500_meta <- metadata %>% group_by(batch) %>% slice_sample(n=500)

# random sample 10% cells of each donor
# rna_subset_meta <- metadata %>% group_by(batch) %>% sample_frac(percent)

# random sample Ncells of data
rna_subset_meta <- metadata %>% slice_sample(n=Ncell)

table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)



# %%
rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load binarized ATAC
ATAC_Dir <- paste0(input_path, "/BMMC-multiome-raw-ATAC-peaks.mtx")
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

metadata <- read.csv(paste0(ATAC_Dir, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)


# %%
atac_subset_meta <- metadata[rownames(metadata) %in% rownames(rna_subset_meta), ]
table(atac_subset_meta$batch)
table(atac_subset_meta$cell_type)

atac_subset_meta <- as.data.frame(atac_subset_meta)
rownames(atac_subset_meta) <- atac_subset_meta$barcode


# %%
# subset atac counts matrix
ATAC_subset_counts <- ATAC_counts[ , colnames(ATAC_counts) %in% atac_subset_meta$barcode]

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

# %%

# %% [markdown]
# ## sample N + N cells

# %%

# %%
input_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData"
output_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c3k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/RawData.json"))

# %%
process = "c5k_c3k"

# %%
# # Load RNA
# RNA_Dir <- paste0(input_path, "/SHARE-multiome-raw-RNA-counts.mtx")
# RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# %%
dim(RNA_counts)

# %%
# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)
metadata

# %%
# set.seed(1234)
# # random sample Ncells of data
# meta_qry <- metadata %>% slice_sample(n=5000)

# %%
metadata_others <- metadata[!(metadata$barcode %in% meta_qry$barcode), ]
meta_ref <- metadata_others %>% slice_sample(n=3000)

meta_qry$data_size <- "c5k"
meta_ref$data_size <- "c3k"

rna_subset_meta <- rbind(meta_qry, meta_ref)

table(rna_subset_meta$data_size)
dim(rna_subset_meta)

# %%
max(table(rna_subset_meta$barcode))

# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# # Load ATAC
# ATAC_Dir <- paste0(input_path, "/SHARE-multiome-raw-ATAC-peaks.mtx")
# ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

# %%
atac_subset_meta <- rna_subset_meta

# subset atac counts matrix
ATAC_subset_counts <- ATAC_counts[ , colnames(ATAC_counts) %in% atac_subset_meta$barcode]

# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
# write_csv(atac_subset_meta, here(output_path, "metadata.csv"))

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

# %% [markdown]
# ## sample (downsample N) + N cells

# %%
library(scuttle)
input_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k/c5k_c5k.json"))

# %%
# 设置提取比例
proportion = 0.10
process = "c5k_R10_A10_c5k"

output_path <- paste0("/home/wsg/BM/data/SHARE/RNA+ATAC/c5k_c5k_robust/", process)

# %%
# # Load RNA
# RNA_Dir <- paste0(input_path, "/SHARE-multiome-c5k_c5k-RNA-counts.mtx")
# RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# # 添加barcode到metadata
# metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
# metadata['barcode'] <- rownames(metadata)

# %%
# metadata_data_1 <- metadata[metadata$data_size == "c5k_1", ]
# metadata_data_2 <- metadata[metadata$data_size == "c5k_2", ]

# table(metadata_data_1$data_size)
# table(metadata_data_2$data_size)

# RNA_counts_data_1  <- RNA_counts[, metadata_data_1$barcode]
# RNA_counts_data_2  <- RNA_counts[, metadata_data_2$barcode]

# %%
RNA_ds = RNA_counts_data_1
RNA_raw = RNA_counts_data_2

# %%
# RNA_counts
dim(RNA_ds)
nnzero(RNA_ds)
sum(RNA_ds)

# %%
RNA_ds_subset <- downsampleMatrix(RNA_ds, prop = proportion, bycol = T)

# %%
# RNA_subset_counts
dim(RNA_ds_subset)
nnzero(RNA_ds_subset)
sum(RNA_ds_subset)

# %%
if(!identical(rownames(RNA_ds_subset), rownames(RNA_raw))) {
  stop("The row names of the two matrices do not match.")
}

RNA_subset_counts <- cbind(RNA_raw, RNA_ds_subset)

if(!identical(colnames(RNA_subset_counts), metadata$barcode)) {
  stop("The row names of the two matrices do not match.")
}

# %%
sum(RNA_counts)
sum(RNA_raw)
sum(RNA_ds)

sum(RNA_subset_counts)
sum(RNA_raw)
sum(RNA_ds_subset)

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%
# # Load binarized ATAC
# ATAC_Dir <- paste0(input_path, "/SHARE-multiome-c5k_c5k-ATAC-peaks.mtx")
# ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

# metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
# # 添加barcode到metadata
# metadata['barcode'] <- rownames(metadata)

# %%
# metadata_data_1 <- metadata[metadata$data_size == "c5k_1", ]
# metadata_data_2 <- metadata[metadata$data_size == "c5k_2", ]

# table(metadata_data_1$data_size)
# table(metadata_data_2$data_size)

# ATAC_counts_data_1  <- ATAC_counts[, metadata_data_1$barcode]
# ATAC_counts_data_2  <- ATAC_counts[, metadata_data_2$barcode]

# %%
ATAC_raw = ATAC_counts_data_1
ATAC_ds = ATAC_counts_data_2

# %%
# ATAC_counts
dim(ATAC_ds)
nnzero(ATAC_ds)
sum(ATAC_ds)

# %%
ATAC_ds_subset <- downsampleMatrix(ATAC_ds, prop = proportion, bycol = T)

# %%
dim(ATAC_ds_subset)
nnzero(ATAC_ds_subset)
sum(ATAC_ds_subset)

# %%
if(!identical(rownames(ATAC_ds_subset), rownames(ATAC_raw))) {
  stop("The row names of the two matrices do not match.")
}

ATAC_subset_counts <- cbind(ATAC_raw, ATAC_ds_subset)

if(!identical(colnames(ATAC_subset_counts), metadata$barcode)) {
  stop("The row names of the two matrices do not match.")
}

# %%
sum(ATAC_counts)
sum(ATAC_raw)
sum(ATAC_ds)

sum(ATAC_subset_counts)
sum(ATAC_raw)
sum(ATAC_ds_subset)

# %%

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = metadata)

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

# %%

# %%

# %% [markdown]
# ## sample N cells (N > rowdata)

# %%
input_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/RawData"
output_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/c500k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/dataset.json"))

# 设置提取数量
Ncell = 500000
process = "c100k"

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-multiome-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

# %%
# random sample Ncells more than raw data
set.seed(1)
rna_subset_meta_1 <- metadata %>% slice_sample(n=50000)
set.seed(2)
rna_subset_meta_2 <- metadata %>% slice_sample(n=50000)
set.seed(3)
rna_subset_meta_3 <- metadata %>% slice_sample(n=50000)
set.seed(4)
rna_subset_meta_4 <- metadata %>% slice_sample(n=50000)
set.seed(5)
rna_subset_meta_5 <- metadata %>% slice_sample(n=50000)
set.seed(6)
rna_subset_meta_6 <- metadata %>% slice_sample(n=50000)
set.seed(7)
rna_subset_meta_7 <- metadata %>% slice_sample(n=50000)
set.seed(8)
rna_subset_meta_8 <- metadata %>% slice_sample(n=50000)
set.seed(9)
rna_subset_meta_9 <- metadata %>% slice_sample(n=50000)
set.seed(10)
rna_subset_meta_10 <- metadata %>% slice_sample(n=50000)

# rna_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2)
rna_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2, rna_subset_meta_3, 
                         rna_subset_meta_4, rna_subset_meta_5, rna_subset_meta_6, 
                         rna_subset_meta_7, rna_subset_meta_8, rna_subset_meta_9, 
                         rna_subset_meta_10)
dim(rna_subset_meta)

# %%
rna_subset_meta <- as.data.frame(rna_subset_meta)

# subset counts matrix
RNA_subset_counts <- RNA_counts[, rna_subset_meta$barcode]
colnames(RNA_subset_counts) <- make.unique(colnames(RNA_subset_counts))

rna_subset_meta$barcode <- make.unique(rna_subset_meta$barcode)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# %%
dim(RNA_subset_counts)

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load binarized ATAC
ATAC_Dir <- paste0(input_path, "/BMMC-multiome-raw-ATAC-peaks.mtx")
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

metadata <- read.csv(paste0(ATAC_Dir, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

# %%
# atac_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2)
atac_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2, rna_subset_meta_3, 
                         rna_subset_meta_4, rna_subset_meta_5, rna_subset_meta_6, 
                         rna_subset_meta_7, rna_subset_meta_8, rna_subset_meta_9, 
                         rna_subset_meta_10)
dim(atac_subset_meta)

table(atac_subset_meta$batch)
table(atac_subset_meta$cell_type)


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

# %% [markdown]
# ## split sites and donors

# %%
input_path <- "/home/wsg/BM/data/BMMC/RNA+ATAC/RawData"
output_path <- "/home/wsg/BM/data/BMMC/RNA+ATAC/s1d2_s3d10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/BMMC/RNA+ATAC/RawData/dataset.json"))

# %%
# 设置提取样本
# s1d1  s1d2  s1d3  s2d1  s2d4  s2d5 s3d10  s3d3  s3d6  s3d7  s4d1  s4d8  s4d9 
# 6224  6740  4279  4220  6111  4895  6781  4325  1679  1771  8023  9876  4325 

# sample = "s4d9"
sample = "s1d2_s3d10"

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-multiome-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# %%
# 添加barcode到metadata
metadata <- read.csv(paste0(RNA_Dir, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

set.seed(1234)

# %%
# rna_subset_meta <- metadata[metadata$batch==sample, ]

samples = strsplit(sample, split = "_")[[1]]
rna_subset_meta <- metadata[metadata$batch==samples[1] | metadata$batch==samples[2], ]

table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)

rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"],
                       dataset["task_type"],
                       sample,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
SeuratDisk::Convert(paste(dataset["data_name"], dataset["task_type"], sample, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ATAC
ATAC_Dir <- paste0(input_path, "/BMMC-multiome-raw-ATAC-peaks.mtx")
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

# %%
# subset atac counts matrix
ATAC_subset_counts <- ATAC_counts[ , colnames(ATAC_counts) %in% rna_subset_meta$barcode]

# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       sample,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = rna_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ATAC_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample, 
                                   "ATAC", "peaks.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"ATAC", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %% [markdown]
# ## dowmsample from single sample

# %%
input_path <- "/home/wsg/BM/data/BMMC/RNA+ATAC/site3/donor10"
output_path <- "/home/wsg/BM/data/BMMC/RNA+ATAC/site3/donor10_c5k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/BMMC/RNA+ATAC/site3/donor10/s3d10.json"))

# 设置提取数量
Ncell = 5000
process = "s3d10c5k"

# %%
# Load RNA
RNA_Dir <- list.files(input_path, pattern = "\\RNA-counts.mtx$", full.names = TRUE)
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# add metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

## random sample 500 cells of each donor
# bmmc_rna_500_meta <- metadata %>% group_by(batch) %>% slice_sample(n=500)

## random sample 10% cells of each donor
# rna_subset_meta <- metadata %>% group_by(batch) %>% sample_frac(percent)

# random sample Ncells of data
rna_subset_meta <- metadata %>% slice_sample(n=Ncell)

table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)

# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))

write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load binarized ATAC
ATAC_Dir <- list.files(input_path, pattern = "\\ATAC-peaks.mtx$", full.names = TRUE)
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

# add metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

atac_subset_meta <- metadata[rownames(metadata) %in% rownames(rna_subset_meta), ]

table(atac_subset_meta$batch)
table(atac_subset_meta$cell_type)

# subset atac counts matrix
ATAC_subset_counts <- ATAC_counts[ , colnames(ATAC_counts) %in% atac_subset_meta$barcode]

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

# %%

# %%

# %% [markdown]
# # Data Manipulation: BMMC CITE-seq

# %%
input_path <- ""
output_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/dataset.json"))

# %%
# Load RNA
BMMC_RNA_Dir <- paste0(input_path, "/BMMC-raw-pair-RNA-counts.mtx")
BMMC_RNA_counts <- Read10X(data.dir = BMMC_RNA_Dir, gene.column = 1)

# %%
# save raw rna to rds
saveRDS(BMMC_RNA_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          "raw",
                          dataset["task_type"], 
                          "RNA", "counts.rds", sep = "-")))

# %%
# Load ADT
BMMC_ADT_Dir <- paste0(input_path, "/BMMC-raw-pair-ADT-counts.mtx")
BMMC_ADT_counts <- Read10X(data.dir = BMMC_ADT_Dir, gene.column = 1)
BMMC_ADT_counts

# %%
# save raw adt to rds
saveRDS(BMMC_ADT_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          "raw",
                          dataset["task_type"], 
                          "ADT", "counts.rds", sep = "-")))

# %% [markdown]
# ## subset 5 percents cells of BMMC CITE-seq

# %%
input_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/RawData"
output_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/RawData/dataset.json"))


# %%
# 设置提取比例
percent = 0.1
process = "p10"

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)


set.seed(1234)

# %%
RNA_counts

# %%
# random sample 500 cells of each donor
# bmmc_rna_500_meta <- metadata %>% group_by(batch) %>% slice_sample(n=500)

# random sample 10% cells of each donor
rna_subset_meta <- metadata %>% group_by(batch) %>% sample_frac(percent)
table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)

rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ADT
ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-ADT-counts.mtx")
ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

adt_subset_meta <- metadata[rownames(metadata) %in% rownames(rna_subset_meta), ]
table(adt_subset_meta$batch)
table(adt_subset_meta$cell_type)

adt_subset_meta <- as.data.frame(adt_subset_meta)
rownames(adt_subset_meta) <- adt_subset_meta$barcode

# %%
# set.seed(1234)
# random sample 500 cells of each donor
# bmmc_rna_500_meta <- metadata %>% group_by(batch) %>% slice_sample(n=500)

# random sample 10% cells of each donor
# bmmc_rna_10p_meta <- metadata %>% group_by(batch) %>% sample_frac(.1)



# subset adt counts matrix
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(adt_subset_meta, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = adt_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "counts.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ADT", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %% [markdown]
# ## subset cells from single sample

# %%
input_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/RawData"
output_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/s2d1"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/dataset.json"))


# %%
# 设置提取样本
#  s1d1  s1d2  s1d3  s2d1  s2d4  s2d5  s3d1  s3d6  s3d7  s4d1  s4d8  s4d9 
#  5227  4978  6106 10465  5584  9122  9521 11035 11473  5456  3929  7365 
sample = "s2d1"

# %%
output_path

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

set.seed(1234)


# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-multiome-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(RNA_Dir, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

set.seed(1234)

# %%
RNA_counts

# %%
rna_subset_meta <- metadata[metadata$batch==sample, ]
table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)

rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"],
                       dataset["task_type"],
                       sample,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ADT
ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-ADT-counts.mtx")
ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

# %%
adt_subset_meta <- metadata[rownames(metadata) %in% rownames(rna_subset_meta), ]
table(adt_subset_meta$batch)
table(adt_subset_meta$cell_type)

adt_subset_meta <- as.data.frame(adt_subset_meta)
rownames(adt_subset_meta) <- adt_subset_meta$barcode

# %%
# subset adt counts matrix
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# save raw adt to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       sample,
                       "ADT", "peaks.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(adt_subset_meta, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "ADT", "peaks.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = adt_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample, 
                                   "ADT", "peaks.h5Seurat", sep = "-")))

# save raw adt to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample, "ADT", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## subset N cells

# %%
input_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/RawData"
output_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/c50k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/dataset.json"))

# 设置提取数量
Ncell = 50000
process = "c50k"

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)


set.seed(1234)

# %%
dim(RNA_counts)

# %%
# random sample Ncells of data
rna_subset_meta <- metadata %>% slice_sample(n=Ncell)

table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)


# %%
rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ADT
ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-ADT-counts.mtx")
ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)


# %%
adt_subset_meta <- metadata[rownames(metadata) %in% rownames(rna_subset_meta), ]
table(adt_subset_meta$batch)
table(adt_subset_meta$cell_type)

adt_subset_meta <- as.data.frame(adt_subset_meta)
rownames(adt_subset_meta) <- adt_subset_meta$barcode

# %%
# subset adt counts matrix
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# save raw adt to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ADT", "peaks.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(adt_subset_meta, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "peaks.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = adt_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "peaks.h5Seurat", sep = "-")))

# save raw adt to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ADT", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %% [markdown]
# ## sample N + N cells

# %%
input_path <- "/home/wsg/BM/data/BMMC/RNA+ADT/RawData"
output_path <- "/home/wsg/BM/data/BMMC/RNA+ADT/c20k_c20k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/BMMC/RNA+ADT/RawData/dataset.json"))

# %%
# 设置提取样本
#  s1d1  s1d2  s1d3  s2d1  s2d4  s2d5  s3d1  s3d6  s3d7  s4d1  s4d8  s4d9 
#  5227  4978  6106 10465  5584  9122  9521 11035 11473  5456  3929  7365 

process = "c20k_c20k"

# %%
# # Load RNA
# RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-RNA-counts.mtx")
# RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# %%
# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

set.seed(1234)

# %%
set.seed(1234)
# random sample Ncells of data
# rna_subset_meta <- metadata %>% slice_sample(n=Ncell)
meta_qry <- metadata %>% slice_sample(n=20000)

# %%
metadata_others <- metadata[!(metadata$barcode %in% meta_qry$barcode), ]
meta_ref <- metadata_others %>% slice_sample(n=20000)

meta_qry$data_size <- "c20k_1"
meta_ref$data_size <- "c20k_2"

rna_subset_meta <- rbind(meta_qry, meta_ref)

table(rna_subset_meta$data_size)
dim(rna_subset_meta)

# %%
max(table(rna_subset_meta$barcode))

# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# # Load ADT
# ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-ADT-counts.mtx")
# ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

# %%
adt_subset_meta <- rna_subset_meta

# subset adt counts matrix
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# save raw adt to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ADT", "peaks.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(adt_subset_meta, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "peaks.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = adt_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "peaks.h5Seurat", sep = "-")))

# save raw adt to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ADT", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## seubset N cells (N > rowdata)

# %%
input_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/RawData"
output_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/c500k"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/dataset.json"))

# 设置提取数量
Ncell = 500000
process = "c500k"

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)


# %%
# random sample Ncells more than raw data
set.seed(1)
rna_subset_meta_1 <- metadata %>% slice_sample(n=50000)
set.seed(2)
rna_subset_meta_2 <- metadata %>% slice_sample(n=50000)
set.seed(3)
rna_subset_meta_3 <- metadata %>% slice_sample(n=50000)
set.seed(4)
rna_subset_meta_4 <- metadata %>% slice_sample(n=50000)
set.seed(5)
rna_subset_meta_5 <- metadata %>% slice_sample(n=50000)
set.seed(6)
rna_subset_meta_6 <- metadata %>% slice_sample(n=50000)
set.seed(7)
rna_subset_meta_7 <- metadata %>% slice_sample(n=50000)
set.seed(8)
rna_subset_meta_8 <- metadata %>% slice_sample(n=50000)
set.seed(9)
rna_subset_meta_9 <- metadata %>% slice_sample(n=50000)
set.seed(10)
rna_subset_meta_10 <- metadata %>% slice_sample(n=50000)

# rna_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2)
rna_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2, rna_subset_meta_3, 
                         rna_subset_meta_4, rna_subset_meta_5, rna_subset_meta_6, 
                         rna_subset_meta_7, rna_subset_meta_8, rna_subset_meta_9, 
                         rna_subset_meta_10)
dim(rna_subset_meta)

# %%
table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)

# %%
rna_subset_meta <- as.data.frame(rna_subset_meta)

# subset counts matrix
RNA_subset_counts <- RNA_counts[, rna_subset_meta$barcode]
colnames(RNA_subset_counts) <- make.unique(colnames(RNA_subset_counts))

rna_subset_meta$barcode <- make.unique(rna_subset_meta$barcode)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# %%
dim(RNA_subset_counts)

# %%
rna_subset_meta <- as.data.frame(rna_subset_meta)
rownames(rna_subset_meta) <- rna_subset_meta$barcode

# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ADT
ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-ADT-counts.mtx")
ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
# 添加barcode到metadata
metadata['barcode'] <- rownames(metadata)

# %%
# adt_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2)
adt_subset_meta <- rbind(rna_subset_meta_1, rna_subset_meta_2, rna_subset_meta_3, 
                         rna_subset_meta_4, rna_subset_meta_5, rna_subset_meta_6, 
                         rna_subset_meta_7, rna_subset_meta_8, rna_subset_meta_9, 
                         rna_subset_meta_10)
dim(adt_subset_meta)

table(adt_subset_meta$batch)
table(adt_subset_meta$cell_type)

# %%
adt_subset_meta <- as.data.frame(adt_subset_meta)

# subset counts matrix
ADT_subset_counts <- ADT_counts[, adt_subset_meta$barcode]
colnames(ADT_subset_counts) <- make.unique(colnames(ADT_subset_counts))

adt_subset_meta$barcode <- make.unique(adt_subset_meta$barcode)
rownames(adt_subset_meta) <- adt_subset_meta$barcode

# %%
dim(ADT_subset_counts)

# %%
# subset adt counts matrix
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# save raw adt to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ADT", "peaks.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(adt_subset_meta, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "peaks.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = adt_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "peaks.h5Seurat", sep = "-")))

# save raw adt to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ADT", "peaks.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %% [markdown]
# ## split sites and donors

# %%
input_path <- "/home/wsg/BM/data/BMMC/RNA+ADT/RawData"
output_path <- "/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/BMMC/RNA+ADT/RawData/dataset.json"))

# %%
# 设置提取样本
#  s1d1  s1d2  s1d3  s2d1  s2d4  s2d5  s3d1  s3d6  s3d7  s4d1  s4d8  s4d9 
#  5227  4978  6106 10465  5584  9122  9521 11035 11473  5456  3929  7365 

# sample = "s4d9"
sample = "s2d1_s3d6"

# %%
# # Load RNA
# RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-RNA-counts.mtx")
# RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# %%
# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = 1)
metadata['barcode'] <- rownames(metadata)

set.seed(1234)

# %%
# rna_subset_meta <- metadata[metadata$batch==sample, ]

samples = strsplit(sample, split = "_")[[1]]
rna_subset_meta <- metadata[metadata$batch==samples[1] | metadata$batch==samples[2], ]

table(rna_subset_meta$batch)
table(rna_subset_meta$cell_type)


# %%
# subset counts matrix
RNA_subset_counts <- RNA_counts[, colnames(RNA_counts) %in% rna_subset_meta$barcode]

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"],
                       dataset["task_type"],
                       sample,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(rna_subset_meta, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = rna_subset_meta)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample,
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample,"RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# Load ADT
ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-raw-ADT-counts.mtx")
ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

# %%
ADT_counts

# %%
adt_subset_meta = rna_subset_meta
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# %%
ADT_subset_counts

# %%
adt_subset_meta = rna_subset_meta

# subset adt counts matrix
ADT_subset_counts <- ADT_counts[ , colnames(ADT_counts) %in% adt_subset_meta$barcode]

# save raw adt to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       sample,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(adt_subset_meta, here(output_path, "metadata.csv"))

# save raw adt to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          sample,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = adt_subset_meta)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   sample, 
                                   "ADT", "counts.h5Seurat", sep = "-")))

# save raw adt to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], sample, "ADT", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%

# %% [markdown]
# # Data Manipulation: Load and Downsample SHARE

# %%
library(scuttle)

# %% [markdown]
# ## Downsample n percents counts of SHARE data

# %%
library(scuttle)
input_path <- "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/SHARE/RNA+ATAC/RawData/RawData.json"))

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/SHARE-multiome-raw-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

# %%
# 设置提取比例
proportion = 0.10
process = "R10_A10"

output_path <- paste0("/home/wsg/BM/data/SHARE/RNA+ATAC/", process)

# %%
# RNA_counts
dim(RNA_counts)
nnzero(RNA_counts)
sum(RNA_counts)

# %%
RNA_subset_counts <- downsampleMatrix(RNA_counts, prop = proportion, bycol = T)

# %%
# RNA_subset_counts
dim(RNA_subset_counts)
nnzero(RNA_subset_counts)
sum(RNA_subset_counts)

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%
# # Load binarized ATAC
# ATAC_Dir <- paste0(input_path, "/SHARE-multiome-raw-ATAC-peaks.mtx")
# ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

# metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
# # 添加barcode到metadata
# metadata['barcode'] <- rownames(metadata)

# %%
# ATAC_counts
dim(ATAC_counts)
nnzero(ATAC_counts)
sum(ATAC_counts)

# %%
ATAC_subset_counts <- downsampleMatrix(ATAC_counts, prop = proportion, bycol = T)

# %%
# ATAC_subset_counts
dim(ATAC_subset_counts)
nnzero(ATAC_subset_counts)
sum(ATAC_subset_counts)

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = metadata)

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

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# ## Downsample n percents counts of BMMC multiome single sample

# %%
library(scuttle)
input_path <- "/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10/s1d1_s3d10.json"))

# %%
# 设置提取比例
proportion = 0.75
process = "s1d1_R75_A75_s3d10_R75_A75"

output_path <- paste0("/home/wsg/BM/data/BMMC/RNA+ATAC/s1d1_s3d10_robust/", process)

# %%
# Load RNA
RNA_Dir <- list.files(input_path, pattern = "\\RNA-counts.mtx$", full.names = TRUE)
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# add metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

# %%
# RNA_counts
dim(RNA_counts)
nnzero(RNA_counts)
sum(RNA_counts)

# %%
RNA_subset_counts <- downsampleMatrix(RNA_counts, prop = proportion, bycol = T)

# %%
# RNA_subset_counts
dim(RNA_subset_counts)
nnzero(RNA_subset_counts)
sum(RNA_subset_counts)

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%
# Load ATAC
ATAC_Dir <- list.files(input_path, pattern = "\\ATAC-peaks.mtx$", full.names = TRUE)
ATAC_counts <- Read10X(data.dir = ATAC_Dir, gene.column = 1)

# add metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

# %%
# ATAC_counts
dim(ATAC_counts)
nnzero(ATAC_counts)
sum(ATAC_counts)

# %%
ATAC_subset_counts <- downsampleMatrix(ATAC_counts, prop = proportion, bycol = T)

# %%
# ATAC_subset_counts
dim(ATAC_subset_counts)
nnzero(ATAC_subset_counts)
sum(ATAC_subset_counts)

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ATAC", "peaks.mtx", sep = "-"))
write10xCounts(x = ATAC_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ATAC_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ATAC", "peaks.rds", sep = "-")))

ATAC_subset <- CreateSeuratObject(counts = ATAC_subset_counts, assay = "ATAC", meta.data = metadata)

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

# %%

# %%

# %%

# %%

# %% [markdown]
# ## Downsample n percents counts of BMMC CITE-seq p10 data

# %%
input_path <- "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/p10"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/p10.json"))

# %%
# Load RNA
RNA_Dir <- paste0(input_path, "/BMMC-CITE_seq-p10-RNA-counts.mtx")
RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

# %%
# 设置提取比例
proportion = 0.5
process = "ds50"

output_path <- paste0("/home/wsg/BM/pipeline/data/BMMC/RNA+ADT/CITE-seq/p10/downsample/", process)

# %%
# RNA_counts
dim(RNA_counts)
nnzero(RNA_counts)
sum(RNA_counts)

# %%
RNA_subset_counts <- downsampleMatrix(RNA_counts, prop = proportion, bycol = T)

# %%
# RNA_subset_counts
dim(RNA_subset_counts)
nnzero(RNA_subset_counts)
sum(RNA_subset_counts)

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%
# Load ADT
ADT_Dir <- paste0(input_path, "/BMMC-CITE_seq-p10-ADT-counts.mtx")
ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

# 添加barcode到metadata
metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
metadata['barcode'] <- rownames(metadata)

# %%
ADT_counts

# %%
downsampleMatrix(ADT_counts, prop = proportion, bycol = T)

# %%
# ADT_counts
dim(ADT_counts)
nnzero(ADT_counts)
sum(ADT_counts)

# %%
ADT_subset_counts <- downsampleMatrix(ADT_counts, prop = proportion, bycol = T)

# %%
# ADT_subset_counts
dim(ADT_subset_counts)
nnzero(ADT_subset_counts)
sum(ADT_subset_counts)

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = metadata)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "counts.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ADT", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %% [markdown]
#  ## Downsample n percents counts of BMMC CITE-seq single sample

# %%
library(scuttle)
input_path <- "/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6"
dataset <- unlist(fromJSON(file = "/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6/s2d1_s3d6.json"))

# %%
# 设置提取比例
proportion = 0.10
process = "s2d1_R10_A10_s3d6_R10_A10"

output_path <- paste0("/home/wsg/BM/data/BMMC/RNA+ADT/s2d1_s3d6_robust/", process)

# %%
# # Load RNA
# RNA_Dir <- list.files(input_path, pattern = "\\RNA-counts.mtx$", full.names = TRUE)
# RNA_counts <- Read10X(data.dir = RNA_Dir, gene.column = 1)

# # add metadata
# metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
# metadata['barcode'] <- rownames(metadata)

# %%
# RNA_counts
dim(RNA_counts)
nnzero(RNA_counts)
sum(RNA_counts)

# %%
RNA_subset_counts <- downsampleMatrix(RNA_counts, prop = proportion, bycol = T)

# %%
# RNA_subset_counts
dim(RNA_subset_counts)
nnzero(RNA_subset_counts)
sum(RNA_subset_counts)

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# save raw rna to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "RNA", "counts.mtx", sep = "-"))
write10xCounts(x = RNA_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw rna to rds
saveRDS(RNA_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "RNA", "counts.rds", sep = "-")))
# Create Seurat Object
RNA_subset <- CreateSeuratObject(counts = RNA_subset_counts, meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(RNA_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "RNA", "counts.h5Seurat", sep = "-")))

# Convert h5Seurat to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "RNA", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%
# # Load ADT
# ADT_Dir <- list.files(input_path, pattern = "\\ADT-counts.mtx$", full.names = TRUE)
# ADT_counts <- Read10X(data.dir = ADT_Dir, gene.column = 1)

# # add metadata
# metadata <- read.csv(paste0(input_path, "/metadata.csv"), row.names = "barcode")
# metadata['barcode'] <- rownames(metadata)

# %%
# ADT_counts
dim(ADT_counts)
nnzero(ADT_counts)
sum(ADT_counts)

# %%
ADT_subset_counts <- downsampleMatrix(ADT_counts, prop = proportion, bycol = T)

# %%
# ADT_subset_counts
dim(ADT_subset_counts)
nnzero(ADT_subset_counts)
sum(ADT_subset_counts)

# %%
# save raw atac to mtx
data_path <- here(output_path,
                 paste(dataset["data_name"], 
                       dataset["task_type"], 
                       process,
                       "ADT", "counts.mtx", sep = "-"))
write10xCounts(x = ADT_subset_counts, path = data_path, version = "3")
write_csv(metadata, here(output_path, "metadata.csv"))

# save raw atac to rds
saveRDS(ADT_subset_counts, 
        file = here(output_path, 
                    paste(dataset["data_name"], 
                          dataset["task_type"], 
                          process,
                          "ADT", "counts.rds", sep = "-")))

ADT_subset <- CreateSeuratObject(counts = ADT_subset_counts, assay = "ADT", meta.data = metadata)

# save raw rna to h5Seurat
SaveH5Seurat(ADT_subset, overwrite = TRUE, 
             filename = here(output_path, 
                             paste(dataset["data_name"], 
                                   dataset["task_type"], 
                                   process, 
                                   "ADT", "counts.h5Seurat", sep = "-")))

# save raw atac to h5ad
setwd(output_path)
Convert(paste(dataset["data_name"], dataset["task_type"], process, "ADT", "counts.h5Seurat", sep = "-"), 
        dest = "h5ad")

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# # Function: convert mtx to rds and h5
#

# %%
mtx_to_rds_h5 <- function(input_path,
                          output_path,
                          dataset){
    # check the species of data
    if (dataset['species'] == "human") { 
        genome = "GRCh38"
    } else if (dataset['species'] == "mouse") {
        genome = "mm10"
    } else {
        stop(paste0("species should be human or mouse, not ", dataset['species']))
    }

    # load the rna data
    rna <- readMM(here(input_path, dataset["gene_expression"]))
    rna_gene <- read.table(here(input_path, dataset["gene_names"]))
    cells_label <- read.table(here(input_path, dataset["gene_barcodes"]))

    rownames(rna) <- rna_gene[, 1]
    colnames(rna) <- cells_label[, 1]
    
    # load the atac data
    atac <- readMM(here(input_path, dataset["atac_expression"]))
    atac_peak <- read.table(here(input_path, dataset["atac_names"]))
    cells_label <- read.table(here(input_path, dataset["atac_barcodes"]))

    rownames(atac) <- atac_peak[, 1]
    colnames(atac) <- cells_label[, 1]
    
    # save raw rna to rds
    saveRDS(rna, 
            file = here(output_path, 
                        paste(dataset["data_name"], 
                              "raw",
                              dataset["task_type"], 
                              "RNA", "count.rds", sep = "-"))
           )
    # save raw atac to rds
    saveRDS(atac, 
            file = here(output_path, 
                        paste(dataset["data_name"], 
                              "raw",
                              dataset["task_type"], 
                              "ATAC", "peaks.rds", sep = "-"))
           )
    
    # save raw rna to h5Seurat
    rna_seurat <- CreateSeuratObject(counts = rna, project = "snare_p0_rna")
    SaveH5Seurat(rna_seurat, overwrite = TRUE, 
                 filename = here(output_path, 

                                 paste(dataset["data_name"], 
                                       "raw", 
                                       dataset["task_type"], 
                                       "RNA", "counts.h5Seurat", sep = "-")))
    # save raw atac to h5Seurat
    chrom_assay <- CreateChromatinAssay(
       counts = atac,
       sep = c("-", "-")
    )
    atac_seurat <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", project = "snare_p0_atac")
    SaveH5Seurat(atac_seurat, overwrite = TRUE, 
                 filename = here(output_path, 
                                 paste(dataset["data_name"], 
                                       "raw", 
                                       dataset["task_type"], 
                                       "ATAC", "peaks.h5Seurat", sep = "-")))
    
    # save raw rna to h5ad
    setwd(output_path)
    Convert(paste(dataset["data_name"],"raw", dataset["task_type"], "RNA", "counts.h5Seurat", sep = "-"), 
            dest = "h5ad")
    # save raw atac to h5ad
    Convert(paste(dataset["data_name"],"raw", dataset["task_type"], "ATAC", "peaks.h5Seurat", sep = "-"), 
            dest = "h5ad")

}

# %%
input_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/RawData" 
output_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/example/snare.json"))

# %%
mtx_to_rds_h5(input_path, output_path, dataset)

# %% [markdown]
# # Function: create rds data of MOFA2
#

# %%
MOFA2_rds_data <- function(input_path,
                           output_path,
                           dataset){
    library(MOFA2)
    library(Seurat)
    library(Signac)
    library(tidyverse)
    
    # check the species of data
    if (dataset['species'] == "human") { 
        genome = "hg38"
        library(EnsDb.Hsapiens.v86)
        lib <- EnsDb.Hsapiens.v86
    } else if (dataset['species'] == "mouse") {
        genome = "mm10"
        library(EnsDb.Mmusculus.v79)
        lib <- EnsDb.Mmusculus.v79
    } else {
        stop(paste0("species should be human or mouse, not ", dataset['species']))
    }
    
    # Load Data
    rna <- readRDS(here(output_path, 
                    paste(dataset["data_name"], 
                          "raw", 
                          dataset["task_type"],
                          "RNA", "counts.rds", sep = "-")
                   )
              )
    
    atac <- readRDS(here(output_path, 
                    paste(dataset["data_name"], 
                          "raw", 
                          dataset["task_type"],
                          "ATAC", "peaks.rds", sep = "-")
                   )
              )
    
    # SCTransform on RNA
    mofa_data <- CreateSeuratObject(counts = rna)
    DefaultAssay(mofa_data) <- "RNA"
    mofa_data <- SCTransform(mofa_data, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    mofa_data <- FindVariableFeatures(mofa_data, selection.method = "vst", nfeatures = 3000) 
    
    # Annotate ATAC
    chrom_assay <- CreateChromatinAssay(counts = atac, sep = c("-", "-"), )
    mofa_data[["ATAC"]] <- chrom_assay    
    
    DefaultAssay(mofa_data) <- "ATAC"
    annotations <- GetGRangesFromEnsDb(ensdb = lib)
    seqlevelsStyle(annotations) <- 'Ensembl'
    genome(annotations) <- genome
    Annotation(mofa_data) <- annotations
    
    # Filter Data
    mofa_data <- subset(x = mofa_data, 
                        subset = nCount_ATAC < 70000 & nCount_ATAC > 10 & 
                        nCount_RNA < 25000 & nCount_RNA > 10)
    
    # RunTFIDF on ATAC
    DefaultAssay(mofa_data) <- "ATAC"
    mofa_data <- RunTFIDF(mofa_data)
    mofa_data <- FindTopFeatures(mofa_data, min.cutoff = 'q98')
    
    # Merge Data
    mofa <- create_mofa(mofa_data, assays = c("SCT","ATAC"))
    print(mofa)
    plot_data_overview(mofa)
    
    # Save Data
    saveRDS(mofa, 
            file = here(output_path, 
                        paste(dataset["data_name"], 
                              "MOFA2",
                              dataset["task_type"], 
                              "multi", "filtered.rds", sep = "-"))
           )   
    
}

# %%
input_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess" 
output_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/BMMC.json"))

# %%
MOFA2_rds_data(input_path, output_path, dataset)

# %% [markdown]
# # Function: create mtx data of scDEC
#

# %%
scDEC_mtx_data <- function(input_path,
                           output_path,
                           dataset){
    # Make Dir
    system(paste0("mkdir -p ", output_path))   
    data_dir=paste(dataset["data_name"], "scDEC", 
                   dataset["task_type"], "multi", "raw.mtx", sep = "-")
    data_path=here(output_path, data_dir)
    system(paste0("mkdir -p ", data_path))
    
    # Load Data
    rna <- as.matrix(readMM(here(input_path, dataset["gene_expression"])))
    rna_gene <- read.table(here(input_path, dataset["gene_names"]))

    atac <- as.matrix(readMM(here(input_path, dataset["atac_expression"])))
    atac_peak <- read.table(here(input_path, dataset["atac_names"]))

    # Merge Data
    merge_mat <- rbind(rna, atac)
    out_tab <- as(as.matrix(merge_mat), "dgCMatrix")

    feat_tab <- data.frame(a = 0, 
                           name = c(rna_gene[,1], 
                                    atac_peak[,1]),
                           group = c(rep("Gene Expression", nrow(rna_gene)),
                                     rep("Peaks", nrow(atac_peak))))

    # Save Data
    writeMM(out_tab, here(data_path, "matrix.mtx"))
    write.table(feat_tab, file = paste0(data_path, "/features.tsv"), 
                row.names = F, col.names = F, sep='\t', quote=F)
    system(paste("cp", 
                 here(input_path, dataset["gene_barcodes"]), 
                 here(data_path, "barcodes.tsv.gz")))
    
}

# %%
input_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess" 
output_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/BMMC.json"))

# %%
scDEC_mtx_data(input_path, output_path, dataset)

# %% [markdown]
# # Function: create rds data of SeuratV4

# %%
SeuratV4_rds_data <- function(input_path, 
                              output_path,
                              dataset){
    library(Seurat)
    library(Signac)
    library(tidyverse)
    
    # check the species of data
    if (dataset['species'] == "human") { 
        genome = "hg38"
        library(EnsDb.Hsapiens.v86)
        lib <- EnsDb.Hsapiens.v86
    } else if (dataset['species'] == "mouse") {
        genome = "mm10"
        library(EnsDb.Mmusculus.v79)
        lib <- EnsDb.Mmusculus.v79
    } else {
        stop(paste0("species should be human or mouse, not ", dataset['species']))
    }
    
    # Load Data
    rna <- readRDS(here(output_path, 
                    paste(dataset["data_name"], 
                          "raw", 
                          dataset["task_type"],
                          "RNA", "counts.rds", sep = "-")
                   )
              )
    
    atac <- readRDS(here(output_path, 
                    paste(dataset["data_name"], 
                          "raw", 
                          dataset["task_type"],
                          "ATAC", "peaks.rds", sep = "-")
                   )
              )
    
    # Merge Data
    wnn_data <- CreateSeuratObject(counts = rna)
    
    # Annotate ATAC
    chrom_assay <- CreateChromatinAssay(counts = atac, sep = c("-", "-"), )
    wnn_data[["ATAC"]] <- chrom_assay
    DefaultAssay(wnn_data) <- "ATAC"
    
    annotations <- GetGRangesFromEnsDb(ensdb = lib)
    seqlevelsStyle(annotations) <- 'Ensembl'
    genome(annotations) <- genome
    Annotation(mofa_data) <- annotations
    
    # Filter Data
    wnn_data <- subset(x = wnn_data,
        subset = nCount_ATAC < 7e4 & nCount_ATAC > 10 &
        nCount_RNA < 25000 & nCount_RNA > 10 
    )
    
    # Save Data
    saveRDS(wnn_data, 
            file = here(output_path, 
                        paste(dataset["data_name"], 
                              "SeuratV4",
                              dataset["task_type"], 
                              "multi", "filtered.rds", sep = "-"))
           )   
    
}

# %%
input_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess" 
output_path <- "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess"
dataset <- unlist(fromJSON(file = "/Data/wangsg/BM/pipeline/results/BMMC/pair/data_preprocess/BMMC.json"))

# %%
SeuratV4_rds_data(input_path, output_path, dataset)

# %%
