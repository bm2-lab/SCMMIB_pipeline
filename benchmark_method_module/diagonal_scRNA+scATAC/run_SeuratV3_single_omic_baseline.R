library(here)
library(rjson)

library(Seurat)
library(Matrix)
library(ggplot2)

task = "ATAC"
input_path = paste0("/home/wsg/BM/pipeline/data/BMMC/RNA+ATAC/multiome/p10/downsample/R10_", task)
output_path = paste0("/home/wsg/BM/pipeline/results/BMMC/RNA+ATAC/multiome/downsample/single_omic_baseline/ATAC/", task)

# Make Dir
if (!dir.exists(output_path)){
    dir.create(output_path)
}

# Load Data
atac_rds <- list.files(input_path, "*ATAC-peaks.rds", full.names = T)
atac <- readRDS(atac_rds)

gtf_file <- "/home/wsg/BM/pipeline/config/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz"
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = atac, 
                                            annotation.file = gtf_file,
                                            seq.levels = c(paste0("chr", 1:22), "chrX", "chrY"), 
                                            upstream = 2000, 
                                            verbose = TRUE)

SeuratV3.atac <- CreateSeuratObject(counts = atac, assay = "ATAC", project = "SeuratV3_ATAC")
SeuratV3.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
SeuratV3.atac$tech <- "atac"

DefaultAssay(SeuratV3.atac) <- "ACTIVITY" 
SeuratV3.atac <- FindVariableFeatures(SeuratV3.atac)
SeuratV3.atac <- NormalizeData(SeuratV3.atac)
SeuratV3.atac <- ScaleData(SeuratV3.atac)

DefaultAssay(SeuratV3.atac) <- "ATAC"
VariableFeatures(SeuratV3.atac) <- names(which(Matrix::rowSums(SeuratV3.atac) > 100))
SeuratV3.atac <- RunLSI(SeuratV3.atac, n = 50, scale.max = NULL)
SeuratV3.atac <- RunUMAP(SeuratV3.atac, reduction = "lsi", dims = 1:15)

# Save Latent
latent <- as.data.frame(SeuratV3.atac@reductions$lsi@cell.embeddings)

write.table(latent, 
            file = here(output_path, 
                        paste("SeuratV3_single_omic_baseline",
                              "ATAC",
                              task, "latent.csv", sep = "-")), 
            row.names = T, 
            col.names = T, 
            sep = ",", quote = F)

# Save UMAP
umap <- as.data.frame(SeuratV3.atac@reductions$umap@cell.embeddings)
colnames(umap) <- c("UMAP1", "UMAP2")

umap[, "cluster"] <- SeuratV3.atac@meta.data$seurat_clusters
write.table(umap, 
            file = here(output_path, 
                        paste("SeuratV3_single_omic_baseline",
                              "ATAC",
                              task, "UMAP.csv", sep = "-")), 
            row.names = T, 
            col.names = T, 
            sep = ",", quote = F)