library(harmony)
library(Seurat)

library(data.table)
library(Matrix)

ADT_h5ad_path = '/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen/spleen-CITE_seq-raw-ADT-counts.rds'
RNA_h5ad_path = '/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen/spleen-CITE_seq-raw-RNA-counts.rds'

cell_meta = read.table("/home/wsg/BM/data/SPATIAL/RNA+ADT/spleen/metadata.csv",sep=",",header=T,row.names=2)

rna_counts = readRDS(RNA_h5ad_path)
adt_counts = readRDS(ADT_h5ad_path) # ADT用不到

options(Seurat.object.assay.version = "v3")
wnn_data<-CreateSeuratObject(counts = rna_counts, meta.data=cell_meta)
DefaultAssay(wnn_data) <- "RNA"
wnn_data <- SCTransform(wnn_data, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

wnn_data <- wnn_data %>% RunHarmony("batch")

harmony_latent =  Embeddings(wnn_data,"harmony")

write.table(harmony_latent,file="../output/harmony/spleen_latent.csv",sep=',',row.names=T,col.names=T,quote=F)