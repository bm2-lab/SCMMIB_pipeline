# %%
.libPaths("/home/wsg/R/conda-library/4.2")
library(here)
library(rjson)

library(Seurat)
library(SeuratDisk)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

# https://satijalab.org/seurat/articles/seurat5_integration_bridge
# https://github.com/satijalab/seurat/issues/3354
# https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
library(glmGamPoi)
library(BPCells)
library(presto)
library(sctransform)
sessionInfo()

args <- commandArgs(T) 
input_path = args[1]
output_path = args[2]
config = unlist(fromJSON(file=args[3]))


###############
RNA_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["rna_rds_filename"]))
ATAC_h5ad_path = paste0(input_path, "/", gsub(".rds",".h5Seurat",config["atac_rds_filename"]))
RNA_data = LoadH5Seurat(RNA_h5ad_path)
ATAC_data = LoadH5Seurat(ATAC_h5ad_path)

if("batch" %in% colnames(RNA_data@meta.data)){
    data_names = unique(RNA_data@meta.data[, "batch"])
    if('s1d1' %in% data_names){
        paired_dataset= "s1d1"
        unpaired_dataset= "s3d10"
    }
    else {
        paired_dataset= "s3d6"
        unpaired_dataset= "s2d1"
    }


}

        


preprocess_RNA=TRUE
preprocess_ADT=TRUE
preprocess=TRUE
imputation=TRUE


# 修改使用的batch！！！
if("data_size" %in% colnames(RNA_data@meta.data)){
    data_names = unique(RNA_data@meta.data[, "data_size"])
    if('c5k' %in% data_names){
        unpaired_dataset = "c5k"
        paired_dataset = data_names[data_names != "c5k"][1]
    }
    else if('c5k_1' %in% data_names){
        unpaired_dataset = "c5k_1"
        paired_dataset = "c5k_2"
    }
    else{
        data_size = strsplit(data_names[1], split = "_")[[1]][1]
        unpaired_dataset = paste0(data_size, "_1")
        paired_dataset = paste0(data_size, "_2")
    }
    RNA_data@meta.data$batch = RNA_data@meta.data$data_size
    ATAC_data@meta.data$batch = ATAC_data@meta.data$data_size


}
print(unpaired_dataset)
print(paired_dataset)

###### 
rna_counts = RNA_data
atac_counts = ATAC_data





ref_batch = paired_dataset
gap_batch = unpaired_dataset
rna_counts_multi = subset(rna_counts,batch==ref_batch) 
atac_counts_multi = subset(atac_counts,batch==ref_batch) 
rna_counts_ref = subset(rna_counts,batch==gap_batch) 
atac_counts_ref = subset(atac_counts,batch==gap_batch) 
new_RNA_names = paste0(colnames(rna_counts_ref@assays$RNA),'_rna')
new_ATAC_names = paste0(colnames(atac_counts_ref@assays$ATAC),'_atac')
rna_counts_ref <- RenameCells(
    rna_counts_ref,
    new.names = new_RNA_names 
)
atac_counts_ref <- RenameCells(
    atac_counts_ref,
    new.names = new_ATAC_names
)
rna_counts = rna_counts_multi@assays$RNA@counts
atac_counts = atac_counts_multi@assays$ATAC@counts
obj.multi <- CreateSeuratObject(counts = rna_counts)
obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

# add the ATAC-seq assay
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
# 人类的注释
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Change style to UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c("-", "-"),
  genome = 'hg38',
  min.cells = 10,
  annotation = annotations
)
# Add the ATAC assay to the multiome object
obj.multi[["ATAC"]] <- chrom_assay

# obj.multi <- subset(
#   x = obj.multi,
#   subset = nCount_ATAC < 7e4 &
#     nCount_ATAC > 5e3 &
#     nCount_RNA < 25000 &
#     nCount_RNA > 1000 &
#     percent.mt < 20
# )
atac_counts = atac_counts_ref@assays$ATAC@counts
atac_counts = atac_counts[rownames(obj.multi@assays$ATAC$counts),]
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Change to UCSC style 
seqlevelsStyle(annotation) <- 'UCSC'
ATAC_assay <- CreateChromatinAssay(
  counts = atac_counts,
#   fragments = fragpath,
  annotation = annotation
)
# Create Seurat sbject
obj.atac  <- CreateSeuratObject(counts = ATAC_assay,assay = 'ATAC')
# obj.atac[['peak.orig']] <- atac_pbmc_assay
# obj.atac <- subset(obj.atac, subset = nCount_ATAC < 7e4 & nCount_ATAC > 2000)
obj.rna = rna_counts_ref
obj.rna = SCTransform(object = obj.rna,min_cells=0)%>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',return.model = TRUE)
DefaultAssay(obj.multi) <- "RNA"
obj.multi <- SCTransform(obj.multi, verbose = TRUE,min_cells=0)
# normalize multiome ATAC
DefaultAssay(obj.multi) <- "ATAC"
obj.multi <- RunTFIDF(obj.multi)
obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
# normalize query
obj.atac <- RunTFIDF(obj.atac)
dims.atac <- 2:50# atac的第一维往往是测序深度
dims.rna <- 1:50
DefaultAssay(obj.multi) <-  "RNA"
DefaultAssay(obj.rna) <- "SCT"
obj.rna.ext <- PrepareBridgeReference(
  reference = obj.rna, bridge = obj.multi,
  reference.reduction = "pca", reference.dims = dims.rna,
  normalization.method = "SCT")

bridge.anchor <- FindBridgeIntegrationAnchors(
  extended.reference = obj.rna.ext, query = obj.atac,
  reduction = "lsiproject", dims = dims.atac)

obj.multi.integrate <- IntegrateEmbeddings(anchorset = bridge.anchor, reference = obj.rna.ext,
  query = obj.atac)

latent_embed_path = paste0(output_path, "/", paste(config["output_prefix"], 
                                               "SeuratV5-multi-latent.csv", sep = "-")) 
latent_embed_path_multi = sub('-multi-latent.csv','_multi_lap_latent.csv',latent_embed_path)
latent_embed_path_rna = sub('-multi-latent.csv','_rna_reduc_latent.csv',latent_embed_path)
latent_embed_path_atac = sub('-multi-latent.csv','_atac_reduc_latent.csv',latent_embed_path)
write.csv(data.frame(obj.rna.ext@reductions$lap@cell.embeddings),latent_embed_path_multi)
write.csv(data.frame(obj.multi.integrate@reductions$integrated_dr@cell.embeddings)[new_RNA_names,],latent_embed_path_rna)
write.csv(data.frame(obj.multi.integrate@reductions$integrated_dr@cell.embeddings)[new_ATAC_names,],latent_embed_path_atac)
# write.csv(data.frame(obj.multi.integrate@reductions$integrated_dr@cell.embeddings),latent_embed_path)
