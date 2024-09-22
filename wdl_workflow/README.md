### WDL usage

cromwell.jar file can be downloaded from https://github.com/broadinstitute/cromwell. <br>
A breif introduction of WDL and cromwell is available at https://cromwell.readthedocs.io/en/stable/CommandLine/. 


### inputs.json usage demostration

```python
{
  "main.input_path": "/home/wsg/BM/data/BMMC/RNA+ATAC/test", # input path of h5ad and rds file
  "main.config": {
    # file names of h5ad files and rds file. rds files are prepared for method in R.
    "rna_h5ad_filename": "BMMC-multiome-test-RNA-counts.h5ad", 
    "atac_h5ad_filename": "BMMC-multiome-test-ATAC-peaks.h5ad", 
    "gam_h5ad_filename": "BMMC-multiome-test-ATAC-gam.h5ad",
    "rna_rds_filename": "BMMC-multiome-test-RNA-counts.rds",
    "atac_rds_filename": "BMMC-multiome-test-ATAC-peaks.rds",
    "gam_rds_filename": "BMMC-multiome-test-ATAC-gam.rds", 
    "fragments_filename": "", # fragment file name
    "metadata": "metadata.csv",
    "barcode_key": "barcode", # column of barcode
    "celltype_key": "cell_type", # column of cell types
    "batch_key": "batch", # column of batch
    "output_prefix": "BMMC-multiome-test-diagonal_scRNA+scATAC", # output prefix for all out files
    "specie": "human", # species, mouse or human 
    "gene_id": "symbol", # column name of gene id used in gtf file below.
    "gtf_file": "/home/wsg/BM/pipeline/config/gencode.v41.chr_patch_hapl_scaff.annotation.gtf.gz",
    "utils_path": "/home/wsg/BM/pipeline/scripts/utils" # path for get_GAM function.
  },
  "main.output_path": "/home/wsg/BM/results/task/diagonal_scRNA+scATAC/scalability/BMMC_GAM/test", # output folder
  "main.method_path": "/home/wsg/BM/pipeline/scripts/task/diagonal_scRNA+scATAC", # path of the algorithm module
  "main.method": {  # setting all algorithms for execution as true. all true method will be run in parallel. At least 1 method should be true.
    "cobolt": false, 
    "DCCA": false, 
    "MOFA2": false,
    "multigrate": false,
    "multigrate_batch": false,
    "multivi": false,
    "multivi_batch": false,
    "scAI": false,
    "scDEC": false,
    "scMDC": false,
    "scMVAE": false,
    "scMVP": false,
    "SeuratV4": false,
    "uniPort": false,
    
    "GLUE": false,
    "LIGER": false,
    "LIGER_INMF": false,
    "LIGER_UINMF": false,
    "pamona": false,
    "SeuratV3": false,
    "unionCom": false,
    "Pseudo": false,

    "bindSC": true,
    "CiteFuse": false,
    "SCALEX": false,
    "sciPENN": false,
    "sciPENN_batch": false,
    "totalVI": false,
    "totalVI_batch": false,

    "DeepMAPS": false,
    "MaxFuse": true, # In this example, only maxfuse will be executed 
    "scJoint": false,
    "scMM": false,
    "SIMBA": false,
    "LIGER_OINMF": false,
    "SeuratV4_RPCA": false,
    "SeuratV3_CCA": false,
    "MultiMAP": false,
    
    "SeuratV5": false,
    "MIDAS": false,
    "StabMap": false,
    "scVAEIT":false,
    "SpatialGlue":false,
    "scmomat":false,
    "mefisto":false
  },

  "main.software": "/home/wsg/BM/pipeline/software" # where cromwell-50.jar and time deposited. 
}

```

