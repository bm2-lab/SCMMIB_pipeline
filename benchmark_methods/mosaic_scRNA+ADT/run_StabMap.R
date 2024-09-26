.libPaths("/home/wsg/R/conda-library/4.2")
library(here)
library(rjson)

args <- commandArgs(T) 
# input_path = args[1]
# output_path = args[2]
# # print(args[3])
# config = unlist(fromJSON(file=args[3]))
# imputation = as.logical(args[4])



# preprocess_ADT=TRUE
# imputation=TRUE


library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(StabMap)
library(SingleCellMultiModal)
library(scran)
library(ExperimentHub)

## Related links
# https://marionilab.github.io/StabMap/articles/stabMap_PBMC_Multiome.html
# https://github.com/MarioniLab/StabMap
# https://github.com/MarioniLab/StabMap2021
# https://www.nature.com/articles/s41587-023-01766-z

StabMap_module <- function(input_path,
                            output_path, 
                            config,
                            impute=TRUE){

    RNA_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["rna_rds_filename"]))
    ADT_h5ad_path = paste0(input_path,"/", gsub(".rds",".h5Seurat",config["adt_rds_filename"]))
    RNA_data = LoadH5Seurat(RNA_h5ad_path)
    ADT_data = LoadH5Seurat(ADT_h5ad_path)

    preprocess_RNA=TRUE
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
        ADT_data@meta.data$batch = ADT_data@meta.data$data_size
    }
    print(unpaired_dataset)
    print(paired_dataset)

    rna_counts = RNA_data
    adt_counts = ADT_data
    rna_counts = as.SingleCellExperiment(rna_counts)
    adt_counts = as.SingleCellExperiment(adt_counts)
    if (preprocess_RNA){
        rna_counts <- logNormCounts(rna_counts)
        # Feature selection
        decomp <- modelGeneVar(rna_counts)
        decomp<-subset(decomp, mean>0.01)
        decomp = decomp[order(decomp$p.value),]
        hvgs <-  intersect(rownames(decomp[1:min(600,length(rownames(decomp))),]), rownames(rna_counts))
        rna_counts <- rna_counts[hvgs,]
        print(length(hvgs))
    }

    adt_counts <- logNormCounts(adt_counts)

    rna_counts = as.Seurat(rna_counts)
    adt_counts = as.Seurat(adt_counts)
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
    rna_counts_multi = as.SingleCellExperiment(rna_counts_multi)
    adt_counts_multi = as.SingleCellExperiment(adt_counts_multi)
    rna_counts_ref = as.SingleCellExperiment(rna_counts_ref)
    adt_counts_ref = as.SingleCellExperiment(adt_counts_ref)
    logcounts_multi = rbind(logcounts(rna_counts_multi), logcounts(adt_counts_multi))
    logcounts_rna = logcounts(rna_counts_ref)
    logcounts_adt = logcounts(adt_counts_ref)
    assay_list = list(
    RNA = logcounts_rna,
    CITE = logcounts_multi,
    ADT = logcounts_adt
    )
    stab = stabMap(assay_list,
                reference_list = c("CITE"),
                plot = FALSE)

    latent_embed_path = paste0(output_path,"/", paste("StabMap",config["output_prefix"], 
                                                "latent.csv", sep = "_")) 
    write.csv(data.frame(stab),latent_embed_path)

    if (impute){
        imp_adt = imputeEmbedding(
            assay_list,
            stab,
            reference = colnames(assay_list[["CITE"]]),
            query = colnames(assay_list[["RNA"]]))

        imp_rna = imputeEmbedding(
            assay_list,
            stab,
            reference = colnames(assay_list[["CITE"]]),
            query = colnames(assay_list[["ADT"]]))

        # latent_embed_path_multi = sub('-multi-latent.csv','_multi_lap_latent.csv',latent_embed_path)
        adt_imputation_path = sub('_latent.csv','_imputation_adt.csv',latent_embed_path)
        rna_imputation_path = sub('_latent.csv','_imputation_rna.csv',latent_embed_path)

        # readr::write_csv(round(as.data.frame(imp_adt$CITE)[rownames(logcounts_adt),],3), adt_imputation_path)
        write.csv(t(round(as.data.frame(imp_adt$CITE)[rownames(logcounts_adt),],3)), adt_imputation_path,row.names=T,col.names=T,quote=F)
        write.csv(t(round(as.data.frame(imp_rna$CITE)[hvgs,],3)), rna_imputation_path,row.names=T,col.names=T,quote=F)
    }
}

StabMap_module(input_path = args[1],
                     output_path = args[2],
                     config = unlist(fromJSON(file=args[3])),
                     impute = as.logical(args[4])
                    )