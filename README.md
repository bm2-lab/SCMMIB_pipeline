# SCMMIB pipeline for 65 benchmark methods

## Introduction
This SCMMIB_pipeline folder is a sub-project of SCMMIB project. The main project github folder is at https://github.com/bm2-lab/SCMMI_Benchmark. 

**SCMMIB_pipeline** contains:
1. [preprocessing_scripts](preprocessing_scripts/): including scripts to generate the gene activity matrix and downsample and scalability simulation datasets from all benchmark datasets, and a introduction to these scripts.

2. [envs](envs/): including conda environment yaml files for 40 benchmark algorithms, and an example to use the env files.

3. [wdl_workflow](wdl_workflow/): including the wdl worklfow file and input json file for all benchmark method, and an example for input json configuration.

4. [benchmark_methods](benchmark_methods/): including task specific module scripts for 6 single-cell multimodal integration types. These scripts can be executed with uniform pipeline of [wdl_workflow](wdl_workflow/) and algorithm specific [envs](envs/).


The output of all benchmark methods can be evaluated with `scmmib` python package in our project folder https://github.com/bm2-lab/SCMMI_Benchmark. 

## Data availablity
The pre-processed project data (h5ad, rds and rSeurat format) is available at figshare folder (https://doi.org/10.6084/m9.figshare.27161451.v1). Simluation datasets in SCMMIB study can be generated from R scripts in the [data_simulation](preprocessing_scripts/data_simulation/) folder.

## Tutorials 
1. Tutorial for reproducing all methods and all tasks in SCMMIB project [tutorial 1](docs/tutorial_scmmib.md).

2. Tutorial for applying scmmib pipeline to new integration methods [tutorial 2](docs/tutorial_new_methods.md).

3. Tutorial for applying scmmib pipeline to new benchmark datasets [tutorial 3](docs/tutorial_new_datasets.md).


## Related manuscript
Our stage 1 proposal manuscript is available at https://springernature.figshare.com/articles/journal_contribution/Benchmarking_single-cell_multi-modal_data_integrations/26789572. 

Our stage 2 manuscript is submitted. 