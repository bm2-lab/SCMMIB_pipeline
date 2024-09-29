### 1. WDL usage

cromwell.jar file can be downloaded from https://github.com/broadinstitute/cromwell. <br>
A breif introduction of WDL and cromwell is available at https://cromwell.readthedocs.io/en/stable/CommandLine/. 

## 2.Demo of wdl workflow execution
A example for running single/multiple tasks with wdl script
```Bash
java -jar cromwell-50.jar run -i inputs.json main.wdl
```
Only inputs.json file is required to be configured for different algorithms and different tasks. 


### 3.Demostration of inputs.json configuration 

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

  "main.software": "/home/wsg/BM/pipeline/software" # where cromwell-50.jar deposited. 
}

```

### 4. Demo output of the workflow
1. output latent files, imputation files (only part of algorithms).
2. Log file ended with time and memory comsumpation, located in the  `./monitor` folder, an example as following:
```Bash
Elapsed time: 29:52.18
Memory usage: 6624512 KB
CPU usage: 104%
```
Algorithms with GPU accelartion will generate a log file in the same path as latent. The content of this log file is as following:
```Bash
,gpu_memory
device 2,3009MiB
```


### 5. Demostration of adding new algorithms to the wdl pipeline

Here we add a new python method call "new_method" to the `taskit.wdl` as follows:
```python
task run_new_method {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_new_method.txt ]
        then
            echo "run_new_method has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_new_method
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate new_method_env ## run with you own environment 
            ## then execute the new method, and summarize the resouce comsumption at the same time.

            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \ 
            python ~{method_path}/run_new_method.py \
            ~{input_path} \
            ~{output_path}/run_new_method \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_new_method.txt

            conda deactivate

        fi
    }
    output {
    }
}
```
Then add the "new_method" to main.wdl in `workflow main` section as follows:
```
if (method["new_method"]) { 
        call taskit.run_uniPort {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }
```
Finally, add the new_method to `inputs.json` in `main.methods` section for execution as `"new_method":true` and set other algorithms as `false`.  