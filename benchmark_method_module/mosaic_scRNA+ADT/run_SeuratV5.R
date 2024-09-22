# %%
.libPaths("/home/wsg/R/conda-library/4.2")
# library(here)
library(rjson)

args <- commandArgs(T) 
input_path = args[1]
output_path = args[2]
print(args[3])
config = unlist(fromJSON(file=args[3]))

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




###############
RNA_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["rna_rds_filename"]))
ADT_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["adt_rds_filename"]))
RNA_data = LoadH5Seurat(RNA_h5ad_path)
ADT_data = LoadH5Seurat(ADT_h5ad_path)


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
    ADT_data@meta.data$batch = ADT_data@meta.data$data_size


}
print(unpaired_dataset)
print(paired_dataset)

##################

rna_counts = RNA_data
adt_counts = ADT_data




rna_counts@assays$RNA # Assay data with 13953 features for 21500 cells
adt_counts@assays$ADT# Assay data with 134 features for 21500 cells
ref_batch = paired_dataset
gap_batch = unpaired_dataset
rna_counts_multi = subset(rna_counts,batch==ref_batch) 
adt_counts_multi = subset(adt_counts,batch==ref_batch) 
rna_counts_ref = subset(rna_counts,batch==gap_batch) 
adt_counts_ref = subset(adt_counts,batch==gap_batch) 
new_RNA_names = paste0(colnames(rna_counts_ref@assays$RNA),'_rna')
new_ADT_names = paste0(colnames(adt_counts_ref@assays$ADT),'_adt')
rna_counts_ref <- RenameCells(
    rna_counts_ref,
    new.names = new_RNA_names 
)
adt_counts_ref <- RenameCells(
    adt_counts_ref,
    new.names = new_ADT_names
)
rna_counts = rna_counts_multi@assays$RNA@counts
adt_counts = adt_counts_multi@assays$ADT@counts
obj.multi <- CreateSeuratObject(counts = rna_counts)
obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

adt_assay <- CreateAssayObject(counts = adt_counts )
obj.multi[["ADT"]] <- adt_assay
# obj.multi <- subset(
#   x = obj.multi,
#   subset = nCount_RNA < 25000 &
#     nCount_RNA > 1000 &
#     percent.mt < 20
# )
adt_counts = adt_counts_ref@assays$ADT@counts
# Create Seurat sbject
obj.adt  <- CreateSeuratObject(counts = adt_counts ,assay = 'ADT')
# obj.atac[['peak.orig']] <- atac_pbmc_assay
# obj.adt <- subset(obj.adt)
obj.rna = rna_counts_ref
obj.rna = SCTransform(object = obj.rna,min_cells=0)%>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',return.model = TRUE)
DefaultAssay(obj.multi) <- "RNA"
obj.multi <- SCTransform(obj.multi, verbose = TRUE,min_cells=0)
# normalize multiome ATAC
DefaultAssay(obj.multi) <- "ADT"
obj.multi <- NormalizeData(obj.multi, normalization.method = 'CLR', margin = 2) %>% ScaleData() 
obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
# normalize query
obj.adt <- NormalizeData(obj.adt, normalization.method = 'CLR', margin = 2) %>% ScaleData() 
# dims.atac <- 2:50# atac的第一维往往是测序深度
dims.rna <- 1:50
dims.adt <- 1:20
DefaultAssay(obj.multi) <-  "RNA"
DefaultAssay(obj.rna) <- "SCT"
obj.rna.ext <- PrepareBridgeReference(
  reference = obj.rna, bridge = obj.multi,
  reference.reduction = "pca", reference.dims = dims.rna,
  supervised.reduction = 'spca',
  bridge.ref.assay = "RNA",
  bridge.query.assay = "ADT",
  normalization.method = "SCT")

bridge.anchor <- FindBridgeIntegrationAnchors(
  extended.reference = obj.rna.ext, query = obj.adt,
  reduction = "pcaproject", dims = dims.adt)


obj.multi.integrate <- IntegrateEmbeddings(anchorset = bridge.anchor, reference = obj.rna.ext,
  query = obj.adt)

latent_embed_path = paste0(output_path,"/", paste(config["output_prefix"], 
                                               "SeuratV5-multi-latent.csv", sep = "-"))  
latent_embed_path_multi = sub('-multi-latent.csv','_multi_lap_latent.csv',latent_embed_path)
latent_embed_path_rna = sub('-multi-latent.csv','_rna_reduc_latent.csv',latent_embed_path)
latent_embed_path_adt = sub('-multi-latent.csv','_adt_reduc_latent.csv',latent_embed_path)
write.csv(data.frame(obj.rna.ext@reductions$lap@cell.embeddings),latent_embed_path_multi)
write.csv(data.frame(obj.multi.integrate@reductions$integrated_dr@cell.embeddings)[new_RNA_names,],latent_embed_path_rna)
write.csv(data.frame(obj.multi.integrate@reductions$integrated_dr@cell.embeddings)[new_ADT_names,],latent_embed_path_adt)

# write.csv(data.frame(obj.multi.integrate@reductions$integrated_dr@cell.embeddings),latent_embed_path)

