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
# SeruatV3 Function for GAM
CreateGeneActivityMatrix <- function(peak.matrix, 
                                     annotation.file, 
                                     seq.levels = c(paste0("chr", 1:22), "chrX", "chrY"),
                                     include.body = TRUE,
                                     upstream = 2000, 
                                     downstream = 0, 
                                     keep.sparse = TRUE, 
                                     verbose = TRUE) {
    .Deprecated(new = "Signac::GeneActivity", 
                msg = paste("CreateGeneActivityMatrix functionality is being moved to Signac.",
                            "Equivalent functionality can be", 
                            "achieved via the Signac::GeneActivity function; for",
                            "more information on Signac, please see https://github.com/timoast/Signac"))
    if (!PackageCheck("GenomicRanges", error = FALSE)) {
        stop("Please install GenomicRanges from Bioconductor.")
    }
    if (!PackageCheck("rtracklayer", error = FALSE)) {
        stop("Please install rtracklayer from Bioconductor.")
    }
    peak.df <- rownames(x = peak.matrix)
    peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":",
        replacement = "-"), split = "-"))
    peak.df <- as.data.frame(x = peak.df)
    colnames(x = peak.df) <- c("chromosome", "start", "end")
    peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
    BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
    gtf <- rtracklayer::import(con = annotation.file)
    gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = "coarse")
    if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
        GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
    }
    gtf.genes <- gtf[gtf$type == "gene"]
    if (include.body) {
        gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
    } else {
        gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream,
            downstream = downstream)
    }
    gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
    keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance ==
        0]
    peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
    gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
    gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
    peak.ids$gene.name <- gene.ids$gene_name
    peak.ids <- as.data.frame(x = peak.ids)
    peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
    annotations <- peak.ids[, c("peak", "gene.name")]
    colnames(x = annotations) <- c("feature", "new_feature")
    if (!keep.sparse) {
        peak.matrix <- as(object = peak.matrix, Class = "matrix")
    }
    all.features <- unique(x = annotations$new_feature)
    if (nbrOfWorkers() > 1) {
        mysapply <- future_sapply
    } else {
        mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
    }
    newmat.list <- mysapply(X = 1:length(x = all.features), FUN = function(x) {
        features.use <- annotations[annotations$new_feature == all.features[[x]],
            ]$feature
        submat <- peak.matrix[features.use, ]
        if (length(x = features.use) > 1) {
            submat <- Matrix::colSums(submat)
        }
        if (keep.sparse) {
            return(as(object = as.matrix(x = submat), Class = "dgCMatrix"))
        } else {
            return(as.matrix(x = submat))
        }
    }, simplify = FALSE)
    newmat = do.call(what = cbind, args = newmat.list)
    newmat <- t(x = newmat)
    rownames(x = newmat) <- all.features
    return(as(object = newmat, Class = "dgCMatrix"))
}

# %%
input_path <- args[1]
config <- unlist(fromJSON(file=args[2]))
data_name = args[3]

# input_path <- "/home/wsg/BM/data/task/scRNA+scATAC/accuracy/SHARE/RawData"
# config <- unlist(fromJSON(file="/home/wsg/BM/data/task/scRNA+scATAC/accuracy/SHARE/RawData/RawData.json"))

# data_name = "SHARE-multiome-raw-ATAC-gam"

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

# %%
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, 
                                            annotation.file = config['gtf_file'], 
                                            seq.levels = seq_levels, 
                                            upstream = 2000, verbose = TRUE)

# %%
metadata <- read.csv(here(input_path, config['metadata']))
rownames(metadata) <- metadata[,'barcode']

# %%
saveRDS(activity.matrix, 
        file = here(input_path, 
                    paste(data_name, "rds", sep = "."))
       )

# Create Seurat Object
GAM <- CreateSeuratObject(counts = activity.matrix, assay = "ATAC", meta.data = metadata)

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
