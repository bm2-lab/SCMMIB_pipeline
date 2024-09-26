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
library(dplyr)
library(ggplot2)

library(Seurat)
library(Signac)
library(tidyverse)

# library(EnsDb.Hsapiens.v86) # hg38
# library(EnsDb.Mmusculus.v79) # mm10
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1

# %%
seurat_module <- function(input_path,
                          output_path, 
                          config){
    # Make Dir
    if (!dir.exists(output_path)){
        dir.create(output_path, 
                   recursive = TRUE)
    }

    # Load Data
    rna_counts <- readRDS(here(input_path, config["rna_rds_filename"]))
    atac_counts <- readRDS(here(input_path, config["atac_rds_filename"]))

    metadata <- read.csv(here(input_path, config["metadata"]), row.names = config["barcode_key"])

    # Create Seurat object
    seurat <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

    if (config['fragments_filename'] != '') {
        # check the species of data
        if (config['specie'] == "human") { 
            genome = "hg38"
            library(EnsDb.Hsapiens.v86)
            EnsDb <- EnsDb.Hsapiens.v86
        } else if (config['specie'] == "mouse") {
            genome = "mm10"
            library(EnsDb.Mmusculus.v79)
            EnsDb <- EnsDb.Mmusculus.v79
        } else {
            stop(paste0("specie should be human or mouse, not ", config['specie']))
        }
        # Now add in the ATAC-seq data
        # we'll only use peaks in standard chromosomes
        grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
        grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
        atac_counts <- atac_counts[as.vector(grange.use), ]
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb)
        seqlevelsStyle(annotations) <- 'UCSC'
        genome(annotations) <- genome

        frag.file <- here(input_path, config["fragments_filename"])
        chrom_assay <- CreateChromatinAssay(
            counts = atac_counts,
            sep = c("-", "-"),
            genome = genome,
            fragments = frag.file,
            min.cells = 10,
            annotation = annotations
         )
    } else {
        chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c("-", "-"))
    }

    seurat[["ATAC"]] <- chrom_assay

    # Filter Data
    # seurat <- subset(x = seurat,
    #                subset = nCount_ATAC < 7e4 &
    #                nCount_ATAC > 5e3 &
    #                nCount_RNA < 25000 &
    #                nCount_RNA > 1000 &
    #                percent.mt < 20)

    # Data preprocess
    ## RNA analysis
    DefaultAssay(seurat) <- "RNA"
    seurat <- SCTransform(seurat, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

    ## ATAC analysis
    DefaultAssay(seurat) <- "ATAC"
    seurat <- RunTFIDF(seurat)
    seurat <- FindTopFeatures(seurat, min.cutoff = 'q0')
    seurat <- RunSVD(seurat)
    seurat <- RunUMAP(seurat, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

    # WNN
    seurat <- FindMultiModalNeighbors(seurat, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    seurat <- RunUMAP(seurat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    seurat <- FindClusters(seurat, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

    ## save UMAP
    umap <- as.data.frame(seurat@reductions$wnn.umap@cell.embeddings)
    colnames(umap) <- c("UMAP1", "UMAP2")
    umap[, "cluster"] <- seurat@meta.data$seurat_clusters
    write.table(umap, 
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4-multi-umap.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)

    ## save graphs
    write.table(seurat@graphs$wsnn,
                file = here(output_path, paste(config["output_prefix"], 
                                               "SeuratV4-multi-graph.csv", sep = "-")), 
                row.names =T, col.names = T, sep=',', quote=F)
}

# %%
seurat_module(input_path = args[1],
              output_path = args[2],
              config = unlist(fromJSON(file=args[3]))
             )

# %%

# %%

# %%

# %%
# input_path <- "/home/wsg/BM/data/test"
# output_path <- "/home/wsg/BM/results/task/scRNA+scATAC/accuracy/SHARE/RawData/rep_1/run_SeuratV4"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/test/c1k.json"))

# %%
# seurat_module(input_path,
#               output_path,
#               config)

# %%
