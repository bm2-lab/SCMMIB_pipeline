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

library(Matrix)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(data.table)

# %%
# input_path <- "/home/wsg/BM/data/test"
# output_path <- "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_MaxFuse"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/test/c1k.json"))

# %%
input_path <- args[1]
output_path <- args[2]
config_path <- args[3]
config <- unlist(fromJSON(file=config_path))

# %%
# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# %%
# Load Data
rna <- readRDS(here(input_path, config["rna_rds_filename"]))
atac <- readRDS(here(input_path, config["atac_rds_filename"]))

# %%
# Load RNA
rna <- CreateSeuratObject(counts = rna)
rna <- PercentageFeatureSet(rna, pattern = "^MT-", col.name = "percent.mt")
rna <- SCTransform(rna, vars.to.regress = "percent.mt", verbose = FALSE)

# Save SCTransform Matrix
rna_SCTransform = as.data.frame(t(rna@assays$SCT@data))

# Convert sparse matrix (dgCMatrix) to matrix
rna_SCTransform <- as(rna_SCTransform, "matrix")

# Save RNA SCTransform Matrix
rna_file <- gzfile(paste(output_path, "scRNA_SCTransform.csv.gz", sep='/'), "w")
write.csv(rna_SCTransform, rna_file, quote = F)
close(rna_file)

# %%
# Run LSI on ATAC 
atac <- CreateSeuratObject(counts = atac)
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac)

# Save LSI Matrix
atac_LSI = atac@reductions$lsi@cell.embeddings[,c(2:50)]
# Save RNA SCTransform Matrix
atac_file <- gzfile(paste(output_path, "scATAC_LSI.csv.gz", sep='/'), "w")
write.csv(atac_LSI, atac_file, quote = F)
close(atac_file)

# %%
