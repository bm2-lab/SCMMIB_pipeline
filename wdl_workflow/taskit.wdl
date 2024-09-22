version 1.0

## This WDL pipeline contains all the tasks needed in the benchmark workflow.  
##
## Software version requirements 
##
## Cromwell version support 
## - Successfully tested on v78
##
# task DEFINITION

task run_bindSC {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_bindSC.txt ]
        then
            echo "run_bindSC has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_bindSC
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_bindSC.R \
            ~{input_path} \
            ~{output_path}/run_bindSC \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_bindSC.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_cobolt {
    input {
        String input_path
        String output_path
        String method_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_cobolt.txt ]
        then
            echo "run_cobolt has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_cobolt

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_cobolt.py \
            ~{input_path} \
            ~{output_path}/run_cobolt \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_cobolt.txt

            conda deactivate

        fi
    }

    output {

    }
}

task run_CiteFuse {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_CiteFuse.txt ]
        then
            echo "run_CiteFuse has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_CiteFuse
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_CiteFuse.R \
            ~{input_path} \
            ~{output_path}/run_CiteFuse \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_CiteFuse.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_DCCA {
    input {
        String input_path
        String output_path
        String method_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_DCCA.txt ]
        then
            echo "run_DCCA has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_DCCA

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate DCCA 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_DCCA.py \
            ~{input_path} \
            ~{output_path}/run_DCCA \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_DCCA.txt

            conda deactivate

        fi
    }

    output {

    }
}

task run_GLUE {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_GLUE.txt ]
        then
            echo "run_GLUE has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_GLUE

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scglue 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_GLUE.py \
            ~{input_path} \
            ~{output_path}/run_GLUE \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_GLUE.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_MOFA2 {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_MOFA2.txt ]
        then
            echo "run_MOFA2 has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_MOFA2

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_MOFA2.R \
            ~{input_path} \
            ~{output_path}/run_MOFA2 \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_MOFA2.txt      

            conda deactivate

        fi

    }

    output {

    }

}

task run_multigrate {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_multigrate.txt ]
        then
            echo "run_multigrate has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_multigrate

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate multigrate 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_multigrate.py \
            ~{input_path} \
            ~{output_path}/run_multigrate \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_multigrate.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_multigrate_batch {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_multigrate_batch.txt ]
        then
            echo "run_multigrate_batch has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_multigrate_batch

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate multigrate 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_multigrate_batch.py \
            ~{input_path} \
            ~{output_path}/run_multigrate_batch \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_multigrate_batch.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_multivi {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_multivi.txt ]
        then
            echo "run_multivi has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_multivi

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scvi-env
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_multivi.py \
            ~{input_path} \
            ~{output_path}/run_multivi \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_multivi.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_multivi_batch {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_multivi_batch.txt ]
        then
            echo "run_multivi_batch has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_multivi_batch

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scvi-env 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_multivi_batch.py \
            ~{input_path} \
            ~{output_path}/run_multivi_batch \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_multivi_batch.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_LIGER_INMF {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_LIGER_INMF.txt ]
        then
            echo "run_LIGER_INMF has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_LIGER_INMF
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate SeuratV3 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_LIGER_INMF.R \
            ~{input_path} \
            ~{output_path}/run_LIGER_INMF \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_LIGER_INMF.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_LIGER_OINMF {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_LIGER_OINMF.txt ]
        then
            echo "run_LIGER_OINMF has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_LIGER_OINMF
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate SeuratV3 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_LIGER_OINMF.R \
            ~{input_path} \
            ~{output_path}/run_LIGER_OINMF \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_LIGER_OINMF.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_LIGER_UINMF {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_LIGER_UINMF.txt ]
        then
            echo "run_LIGER_UINMF has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_LIGER_UINMF
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate SeuratV3 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_LIGER_UINMF.R \
            ~{input_path} \
            ~{output_path}/run_LIGER_UINMF \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_LIGER_UINMF.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_MaxFuse {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_MaxFuse.txt ]
        then
            echo "run_MaxFuse has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_MaxFuse

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate maxfuse
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_MaxFuse.py \
            ~{input_path} \
            ~{output_path}/run_MaxFuse \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_MaxFuse.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_MultiMAP {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_MultiMAP.txt ]
        then
            echo "run_MultiMAP has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_MultiMAP

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate MultiMAP
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_MultiMAP.py \
            ~{input_path} \
            ~{output_path}/run_MultiMAP \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_MultiMAP.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_pamona {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_pamona.txt ]
        then
            echo "run_pamona has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_pamona

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_pamona.py \
            ~{input_path} \
            ~{output_path}/run_pamona \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_pamona.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_scAI {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scAI.txt ]
        then
            echo "run_scAI has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scAI

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_scAI.R \
            ~{input_path} \
            ~{output_path}/run_scAI \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scAI.txt      

            conda deactivate

        fi

    }

    output {

    }

}

task run_SCALEX {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command <<<
        if [ -f ~{output_path}/monitor/run_SCALEX.txt ]
        then
            echo "run_SCALEX has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SCALEX

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate SCALEX
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_SCALEX.py \
            ~{input_path} \
            ~{output_path}/run_SCALEX \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SCALEX.txt

            conda deactivate

        fi

    >>>

    output {

    }
}

task run_scDEC {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scDEC.txt ]
        then
            echo "run_scDEC has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scDEC
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_scDEC.py \
            ~{input_path} \
            ~{output_path}/run_scDEC \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scDEC.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_sciPENN {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_sciPENN.txt ]
        then
            echo "run_sciPENN has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_sciPENN

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_sciPENN.py \
            ~{input_path} \
            ~{output_path}/run_sciPENN \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_sciPENN.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_sciPENN_batch {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_sciPENN_batch.txt ]
        then
            echo "run_sciPENN_batch has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_sciPENN_batch

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_sciPENN_batch.py \
            ~{input_path} \
            ~{output_path}/run_sciPENN_batch \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_sciPENN_batch.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_scMDC {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scMDC.txt ]
        then
            echo "run_scMDC has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scMDC
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_scMDC.py \
            ~{input_path} \
            ~{output_path}/run_scMDC \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scMDC.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_scMVAE {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scMVAE.txt ]
        then
            echo "run_scMVAE has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scMVAE

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_scMVAE.py \
            ~{input_path} \
            ~{output_path}/run_scMVAE \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scMVAE.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_scMVP {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scMVP.txt ]
        then
            echo "run_scMVP has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scMVP

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_scMVP.py \
            ~{input_path} \
            ~{output_path}/run_scMVP \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scMVP.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_SeuratV3_CCA {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_SeuratV3_CCA.txt ]
        then
            echo "run_SeuratV3_CCA has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SeuratV3_CCA
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate SeuratV3 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_SeuratV3_CCA.R \
            ~{input_path} \
            ~{output_path}/run_SeuratV3_CCA \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SeuratV3_CCA.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_SeuratV4_RPCA {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_SeuratV4_RPCA.txt ]
        then
            echo "run_SeuratV4_RPCA has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SeuratV4_RPCA
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_SeuratV4_RPCA.R \
            ~{input_path} \
            ~{output_path}/run_SeuratV4_RPCA \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SeuratV4_RPCA.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_SeuratV4 {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_SeuratV4.txt ]
        then
            echo "run_SeuratV4 has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SeuratV4
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_SeuratV4.R \
            ~{input_path} \
            ~{output_path}/run_SeuratV4 \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SeuratV4.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_SIMBA {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_SIMBA.txt ]
        then
            echo "run_SIMBA has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SIMBA

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate simba 
            
            export TMPDIR=/home/wsg/BM/pipeline/script/method/RNA+ATAC/unpair/simba

            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_SIMBA.py \
            ~{input_path} \
            ~{output_path}/run_SIMBA \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SIMBA.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_totalVI {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_totalVI.txt ]
        then
            echo "run_totalVI has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_totalVI

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scvi-env 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_totalVI.py \
            ~{input_path} \
            ~{output_path}/run_totalVI \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_totalVI.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_totalVI_batch {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_totalVI_batch.txt ]
        then
            echo "run_totalVI_batch has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_totalVI_batch

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scvi-env 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_totalVI_batch.py \
            ~{input_path} \
            ~{output_path}/run_totalVI_batch \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_totalVI_batch.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_unionCom {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_unionCom.txt ]
        then
            echo "run_unionCom has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_unionCom

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_unionCom.py \
            ~{input_path} \
            ~{output_path}/run_unionCom \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_unionCom.txt

            conda deactivate

        fi

    }

    output {

    }
}

task run_uniPort {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_uniPort.txt ]
        then
            echo "run_uniPort has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_uniPort

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate benchmark
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            python ~{method_path}/run_uniPort.py \
            ~{input_path} \
            ~{output_path}/run_uniPort \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_uniPort.txt

            conda deactivate

        fi

    }

    output {

    }
}


task run_mefisto {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_mefisto.txt ]
        then
            echo "run_mefisto has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_mefisto

            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate r4 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_mefisto.R \
            ~{input_path} \
            ~{output_path}/run_mefisto \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_mefisto.txt      

            conda deactivate

        fi

    }

    output {

    }

}



task run_SeuratV5 {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_SeuratV5.txt ]
        then
            echo "run_SeuratV5 has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SeuratV5
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate SeuratV5 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_SeuratV5.R \
            ~{input_path} \
            ~{output_path}/run_SeuratV5 \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SeuratV5.txt

            conda deactivate

        fi

    }

    output {

    }
}


task run_StabMap {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_StabMap.txt ]
        then
            echo "run_StabMap has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_StabMap
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate StabMap 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_StabMap.R \
            ~{input_path} \
            ~{output_path}/run_StabMap \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_StabMap.txt

            conda deactivate

        fi

    }

    output {

    }
}


task run_MIDAS {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_MIDAS.txt ]
        then
            echo "run_MIDAS has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_MIDAS
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scmidas 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_MIDAS.R \
            ~{input_path} \
            ~{output_path}/run_MIDAS \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_MIDAS.txt

            conda deactivate

        fi

    }

    output {

    }
}


task run_scVAEIT {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scVAEIT.txt ]
        then
            echo "run_scVAEIT has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scVAEIT
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scVAEIT 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_scVAEIT.R \
            ~{input_path} \
            ~{output_path}/run_scVAEIT \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scVAEIT.txt

            conda deactivate

        fi

    }

    output {

    }
}


task run_scmomat {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_scmomat.txt ]
        then
            echo "run_scmomat has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_scmomat
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate scmomat 
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_scmomat.R \
            ~{input_path} \
            ~{output_path}/run_scmomat \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_scmomat.txt

            conda deactivate

        fi

    }

    output {

    }
}


task run_SpatialGlue {
    input {
        String method_path
        String input_path
        String output_path
        String software

        Map[String, String] config
    }
    
    command {
        if [ -f ~{output_path}/monitor/run_SpatialGlue.txt ]
        then
            echo "run_SpatialGlue has already been successfully executed and therefore skipped"
        else
            mkdir -p ~{output_path}/run_SpatialGlue
            
            source ~/software/miniconda3/etc/profile.d/conda.sh
            conda activate spglue_env
            
            ~{software}/time -f 'Elapsed time: %E\nMemory usage: %M KB\nCPU usage: %P' \
            Rscript ~{method_path}/run_SpatialGlue.R \
            ~{input_path} \
            ~{output_path}/run_SpatialGlue \
            ~{write_json(config)} \
            &> ~{output_path}/monitor/run_SpatialGlue.txt

            conda deactivate

        fi

    }

    output {

    }
}


