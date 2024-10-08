


# Must import python in mefisto env, as mofa/mefisto requires correction mofapy2 version!
Sys.setenv(RETICULATE_PYTHON = "/home/shaliu_fu/miniconda3/envs/mefisto_env/bin/python")
# .libPaths("/home/wsg/R/conda-library/4.2")
library(MOFA2)
library(tidyverse)
library(cowplot)
library(magrittr)

# library(EnsDb.Hsapiens.v86) # hg38
# library(EnsDb.Mmusculus.v79) # mm10
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(stringr)

library(Seurat)
library(SeuratDisk)
# library(anndata)
library(here)
.libPaths()
packageVersion('MOFA2') # >1.1.4

# ADT_h5ad_path = '/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node/lymph_node-CITE_seq-raw-ADT-counts.rds'
# RNA_h5ad_path = '/home/wsg/BM/data/SPATIAL/RNA+ADT/lymph_node/lymph_node-CITE_seq-raw-RNA-counts.rds'
# output_path = "/home/sirm/project/SCMMIB/task/spatial_scRNA+ADT/accuracy/lymph_node/mefisto/RUN_1" # embedding

mefisto_module <- function(input_path, 
                        output_path, 
                        config
                       ){
  # Make Dir
  if (!dir.exists(output_path)){
      dir.create(output_path, 
                  recursive = TRUE)
  }

  # Load Data
  rna_counts <- readRDS(here(input_path, config["rna_rds_filename"]))
  adt_counts <- readRDS(here(input_path, config["adt_rds_filename"]))

  # rna_counts = readRDS(RNA_h5ad_path)
  # adt_counts = readRDS(ADT_h5ad_path)

  spatial_loc = rna_counts@meta.data[,c("X","Y")]
  head(spatial_loc)
  adt_counts
  obj.multi <- CreateSeuratObject(counts = rna_counts@assays[['Spatial_RNA']]@counts,meta.data = rna_counts@meta.data)
  obj.multi[["ADT"]] <- CreateAssayObject(counts = adt_counts@assays[['Spatial_ADT']]@counts ,assay = 'ADT')
  # obj.multi <- CreateSeuratObject(counts = rna_counts@assays[['RNA']]@counts)

  # obj.multi[["percent.mt"]] <- PercentageFeatureSet(obj.multi, pattern = "^MT-")

  # obj.multi[["ADT"]] <- CreateSeuratObject(counts = adt_counts@assays[['ADT']]@counts)

  # rna_counts = readRDS(RNA_rds_path)
  # adt_counts = readRDS(ADT_rds_path)
  obj.multi
  DefaultAssay(obj.multi) <- "RNA"
  obj.multi  <- SCTransform(obj.multi , verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  # rna_scale=obj.multi [["SCT"]]@scale.data
  # rna_scale_gene<-rownames(rna_scale)
  obj.multi <- FindVariableFeatures(obj.multi, selection.method = "vst", nfeatures = 3000) 

  DefaultAssay(obj.multi) <- "ADT"
  obj.multi <- NormalizeData(obj.multi, normalization.method = 'CLR', margin = 2) %>% ScaleData() 
  obj.multi<- FindVariableFeatures(obj.multi)
  mofa <- create_mofa(obj.multi, assays = c("SCT","ADT"))
  spatial_loc2 = data.frame(spatial_loc)
  rownames(spatial_loc2)=colnames(rna_counts)
  colnames(spatial_loc2)=c("coord1","coord2")
  mofa <- set_covariates(mofa, t(spatial_loc2))
  data_opts <- get_default_data_options(mofa)

  model_opts <- get_default_model_options(mofa)
  model_opts$num_factors <- 4

  train_opts <- get_default_training_options(mofa)
  train_opts$maxiter <- 1000

  mefisto_opts <- get_default_mefisto_options(mofa)

  mofa <- prepare_mofa(mofa, model_options = model_opts,
                    mefisto_options = mefisto_opts,
                    training_options = train_opts,
                    data_options = data_opts)


  mofa <- run_mofa(mofa,outfile =paste0(output_path,"/mefisto_model.hdf5"),use_basilisk = FALSE)

  mofa <-load_model(paste0(output_path,"/mefisto_model.hdf5"),remove_inactive_factors = FALSE)
  factors <- 1:get_dimensions(mofa)[["K"]]

  mofa <- run_umap(mofa, 
    factors = factors, 
    n_neighbors = 15,  
    min_dist = 0.30
  )

  out_tab=mofa@expectations$Z$group1
  write.table(out_tab,file = paste0(output_path,"/mefisto_latent.csv"),sep=",",row.names=T,col.names=T)
}

SeuratV5_module(input_path = args[1],
                     output_path = args[2],
                     config = unlist(fromJSON(file=args[3]))
                    )
