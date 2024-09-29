## Tutorial 2: Step by step tutorial to use SCMMIB pipeline to evaluate new algorithms

### Step 1. Install your new method

### Step2. Prepare the input datasets

- For accuracy tasks, users can downloaded the processed datasets in public figshare folder: https://figshare.com/account/home#/projects/221476.
- For robustness and scalability simulation, users can generated the simulation datasets with downloaded process datasets and `data_preprocess_R.R` for each dataset. The conda environment for preprocess code is same as scmmib package: [scmmib_env file](https://github.com/bm2-lab/SCMMI_Benchmark/blob/main/scmmib_env.yml). 


### Step3. Run you own methods for all tasks with WDL workflow

1. Add your new method to WDL workflow by steps in section 5 of [documents](../wdl_workflow/README.md):  
2. And then run with command line as:
```Bash
java -jar cromwell-50.jar run -i inputs.json main.wdl
```

### Step4. Calculate the metrics for selected tasks by `scmmib` package

The scmmib package is available at  [github](https://github.com/bm2-lab/SCMMI_Benchmark/). The installation manual is available at [documents](https://github.com/bm2-lab/SCMMI_Benchmark/blob/main/README.md).

We provided three examples for differenet tasks at: [demos](https://github.com/bm2-lab/SCMMI_Benchmark/blob/main/docs/scmmib_demo.py). <br> 

The reference manual for all functions and parameters is available at [documents](https://github.com/bm2-lab/SCMMI_Benchmark/blob/main/docs/scmmib_py_manual.md). <br>

### Step5. Plot summary table for new method along with SCMMIB benchmark methods
The metrics summary tables for all SCMMIB datasets are available at [scmmib github folder](https://github.com/bm2-lab/SCMMI_Benchmark/tree/main/manuscript_figure_script_and_data/stage2_res/SCMMIB_metrics_final).<br>

And the analysis code for summary table in each integration task at [reproducible folder](https://github.com/bm2-lab/SCMMI_Benchmark/tree/main/manuscript_figure_script_and_data/stage2_script). Users can compared the metrics output of new_method with our benchmark methods with the reproducible scripts. <br>

A simple example of output figure in our manuscript is here:
![rank_plot](./pair_RNA_ATAC_robustness.png).

Then compare the new method with other SCMMIB benchmark method following [a quick demo](https://github.com/bm2-lab/SCMMI_Benchmark/blob/main/docs/scmmib_summary_table_demo.r) and [function manuals](https://github.com/bm2-lab/SCMMI_Benchmark/blob/main/docs/scmmib_tab_r_manual.md) in `scmmib` package. 