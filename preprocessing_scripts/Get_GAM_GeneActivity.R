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

library(Seurat)
library(Signac)
library(SeuratDisk)

library(patchwork)
library(pbapply)
library(future)

# %%
input_path <- args[1]
config <- unlist(fromJSON(file=args[2]))
data_name = args[3]

# input_path <- "/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/10x_mouse_brain/RNA+ATAC/2p5/2p5.json"))
# data_name = "brain-multiome-2p5-ATAC-gam"

# input_path <- "/home/wsg/BM/data/BMMC/RNA+ATAC/c1k"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/BMMC/RNA+ATAC/c1k/c1k.json"))
# data_name = "BMMC-multiome-c1k-ATAC-gam"

# %%
peaks <- readRDS(here(input_path, config["atac_rds_filename"]))

# %%
# check the species of data
if (config['specie'] == "human") { 
    seq_levels = c(paste0("chr", 1:22), "chrX", "chrY")
} else if (config['specie'] == "mouse") {
    seq_levels = c(paste0("chr", 1:19), "chrX", "chrY")
} else {
    stop(paste0("species should be human or mouse, not ", dataset['species']))
}

# check the species of data
if (config['specie'] == "human") { 
    genome = "hg38"
    library(EnsDb.Hsapiens.v86)
    lib <- EnsDb.Hsapiens.v86
} else if (config['specie'] == "mouse") {
    genome = "mm10"
    library(EnsDb.Mmusculus.v79)
    lib <- EnsDb.Mmusculus.v79
} else {
    stop(paste0("species should be human or mouse, not ", config['specie']))
}

# %%
metadata <- read.csv(here(input_path, config["metadata"]))
rownames(metadata) <- metadata$barcode

# %%
chrom_assay <- CreateChromatinAssay(
  counts = peaks,
  sep = c("-", "-"),
  fragments = here(input_path, config["fragments_filename"]),
  min.cells = 0,
  min.features = 0
)

seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# seurat[['peaks']]

# granges(seurat)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = lib)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- genome

# add the gene information to the object
Annotation(seurat) <- annotations

# %%
gene.activities <- GeneActivity(seurat)

# %%
saveRDS(gene.activities, 
        file = here(input_path, 
                    paste(data_name, "rds", sep = "."))
       )

# Create Seurat Object
GAM <- CreateSeuratObject(counts = gene.activities, assay = "ATAC", meta.data = metadata)

# save Seurat to h5Seurat
SaveH5Seurat(GAM, overwrite = TRUE, 
             filename = here(input_path, 
                             paste(data_name, "h5Seurat", sep = "."))
            )

# Convert h5Seurat to h5ad
setwd(input_path)
Convert(paste(data_name, "h5Seurat", sep = "."), 
        dest = "h5ad")

# %%
