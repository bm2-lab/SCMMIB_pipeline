.libPaths("/home/wsg/R/conda-library/4.2")
library(here)
library(rjson)

args <- commandArgs(T) 
input_path = args[1]
output_path = args[2]
# print(args[3])
config = unlist(fromJSON(file=args[3]))
imputation = as.logical(args[4])
# print(imputation)

preprocess_RNA=TRUE

library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(StabMap)
library(SingleCellMultiModal)
library(scran)
library(ExperimentHub)

RNA_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["rna_rds_filename"]))
ATAC_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["atac_rds_filename"]))
RNA_data = LoadH5Seurat(RNA_h5ad_path)
ATAC_data = LoadH5Seurat(ATAC_h5ad_path)

print("hello")

if ("batch" %in% colnames(RNA_data@meta.data)){
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
if ("data_size" %in% colnames(RNA_data@meta.data)){
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


rna_counts = RNA_data
atac_counts = ATAC_data
rna_counts = as.SingleCellExperiment(rna_counts)
atac_counts = as.SingleCellExperiment(atac_counts)
if (preprocess_RNA){
    rna_counts <- logNormCounts(rna_counts)
    # Feature selection
    decomp <- modelGeneVar(rna_counts)
    decomp<-subset(decomp, mean>0.01)
    decomp = decomp[order(decomp$p.value),]
    hvgs <-  intersect(rownames(decomp[1:min(600,length(rownames(decomp))),]), rownames(rna_counts))
    # hvgs <- rownames(decomp)[decomp$mean>0.01 & decomp$p.value <= 0.05]
    rna_counts <- rna_counts[hvgs,]
    print(length(hvgs))
}
atac_counts <- logNormCounts(atac_counts)

# Feature selection using highly variable peaks
# And adding matching peaks to genes
decomp <- modelGeneVar(atac_counts)
decomp<-subset(decomp, mean>0.01)
decomp = decomp[order(decomp$p.value),]
# hvgs <- rownames(decomp[1:min(600,length(rownames(decomp))),])
hvgs <-  intersect(rownames(decomp[1:min(600,length(rownames(decomp))),]), rownames(atac_counts))
# hvgs <- rownames(decomp)[decomp$mean>0.05
                        #  & decomp$p.value <= 0.05]
print(length(hvgs))

atac_counts <- atac_counts[hvgs,]
rna_counts = as.Seurat(rna_counts)
atac_counts = as.Seurat(atac_counts)
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
rna_counts_multi = as.SingleCellExperiment(rna_counts_multi)
atac_counts_multi = as.SingleCellExperiment(atac_counts_multi)
rna_counts_ref = as.SingleCellExperiment(rna_counts_ref)
atac_counts_ref = as.SingleCellExperiment(atac_counts_ref)
logcounts_multi = rbind(logcounts(rna_counts_multi), logcounts(atac_counts_multi))
logcounts_rna = logcounts(rna_counts_ref)
logcounts_atac = logcounts(atac_counts_ref)
assay_list = list(
  RNA = logcounts_rna,
  Multiome = logcounts_multi,
  ATAC = logcounts_atac
)
stab = stabMap(assay_list,
               reference_list = c("Multiome"),
               plot = FALSE)


latent_embed_path = paste0(output_path,"/", paste("StabMap",config["output_prefix"], 
                                               "latent.csv", sep = "_")) 
write.csv(data.frame(stab),latent_embed_path)



if (imputation){
    imp_atac = imputeEmbedding(
        assay_list,
        stab,
        reference = colnames(assay_list[["Multiome"]]),
        query = colnames(assay_list[["RNA"]]))

    imp_rna = imputeEmbedding(
        assay_list,
        stab,
        reference = colnames(assay_list[["Multiome"]]),
        query = colnames(assay_list[["ATAC"]]))

    # latent_embed_path_multi = sub('-multi-latent.csv','_multi_lap_latent.csv',latent_embed_path)
    atac_imputation_path = sub('_latent.csv','_imputation_atac.csv',latent_embed_path)
    rna_imputation_path = sub('_latent.csv','_imputation_rna.csv',latent_embed_path)

    write.csv(t(round(as.data.frame(imp_atac$Multiome)[rownames(logcounts_atac),],3)), atac_imputation_path,row.names=T,col.names=T,quote=F)
    write.csv(t(round(as.data.frame(imp_rna$Multiome)[rownames(logcounts_rna),],3)), rna_imputation_path,row.names=T,col.names=T,quote=F)
}


