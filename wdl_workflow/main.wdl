version 1.0

## This WDL pipeline contains execute the benchmark workflow based on your input.json file.  
##
## Software version requirements 
##
## Cromwell version support 
## - Successfully tested on v78
##
# WORKFLOW DEFINITION 

import "/home/wsg/BM/pipeline/workflow/taskit.wdl" as taskit

workflow main {
    input {
        String method_path
        String input_path
        String output_path
        String software
         
        Map[String, String] config
        Map[String, Boolean] method
    }  

    call monitor {
        input:
            output_path = output_path
    } 

    if (method["Pseudo"]) { 
        call Pseudo {
            input: 
                input_path = input_path,
                output_path = output_path,
                method_path = method_path,
                software = software,
                config = config
                
        }
    }

    if (method["bindSC"]) {
        call taskit.run_bindSC {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["cobolt"]) { 
        call taskit.run_cobolt {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["CiteFuse"]) { 
        call taskit.run_CiteFuse {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["DCCA"]) { 
        call taskit.run_DCCA {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    } 

    if (method["GLUE"]) {
        call taskit.run_GLUE {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["MOFA2"]) {
        call taskit.run_MOFA2 {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["multigrate"]) {
        call taskit.run_multigrate {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["multigrate_batch"]) {
        call taskit.run_multigrate_batch {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["multivi"]) {
        call taskit.run_multivi {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["multivi_batch"]) {
        call taskit.run_multivi_batch {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["LIGER_INMF"]) {
        call taskit.run_LIGER_INMF {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["LIGER_OINMF"]) {
        call taskit.run_LIGER_OINMF {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["LIGER_UINMF"]) {
        call taskit.run_LIGER_UINMF {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["MaxFuse"]) { 
        call taskit.run_MaxFuse {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["MultiMAP"]) { 
        call taskit.run_MultiMAP {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["pamona"]) { 
        call taskit.run_pamona {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["scAI"]) {
        call taskit.run_scAI {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["SCALEX"]) { 
        call taskit.run_SCALEX {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["scDEC"]) {
        call taskit.run_scDEC {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["sciPENN"]) { 
        call taskit.run_sciPENN {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["sciPENN_batch"]) { 
        call taskit.run_sciPENN_batch {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["scMDC"]) {
        call taskit.run_scMDC {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["scMVAE"]) { 
        call taskit.run_scMVAE {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["scMVP"]) { 
        call taskit.run_scMVP {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["SeuratV3_CCA"]) {
        call taskit.run_SeuratV3_CCA {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["SeuratV4_RPCA"]) {
        call taskit.run_SeuratV4_RPCA {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["SeuratV4"]) { 
        call taskit.run_SeuratV4 {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["SIMBA"]) { 
        call taskit.run_SIMBA {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["totalVI"]) { 
        call taskit.run_totalVI {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["totalVI_batch"]) { 
        call taskit.run_totalVI_batch {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["unionCom"]) { 
        call taskit.run_unionCom {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

    if (method["uniPort"]) { 
        call taskit.run_uniPort {
            input: 
                input_path = input_path,
                output_path = output_path,
                config = config,
                method_path = method_path,
                software = software
        }
    }

}

# task monitor
task monitor{

    input{
        String output_path
    }

	command{
		mkdir -p ${output_path}/monitor
	}	

}

# task Pseudo
task Pseudo {
    input {
        String method_path
        String input_path
        String output_path
        String software
        Map[String, String] config
    }
    
    command {
        
    }
    
    output {

    }
}