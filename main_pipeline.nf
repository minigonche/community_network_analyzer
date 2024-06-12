#!/usr/bin/env nextflow

process split_files {
    input:
    val cycle
    path input_folder
    path metadata_file
    val split_by_lc

    output:
    path "*.csv"

    script:
    """            
    file_splitter.R --cycle=$cycle --input_folder=$input_folder --metadata_file=$metadata_file --split_by_lc=$split_by_lc   
    """    

}

process build_network {

    publishDir "${params.outDir}", mode: 'copy'

    input:
    path input_file
    path metadata_file
    val metadata_cols
    val quantitative

    
    output:
    path "*.graphml", emit: graph
    path "*.RData", optional: true

    script:
    """            
    network_builder.R --input_file=$input_file --metadata_file=$metadata_file --metadata_cols=$metadata_cols --quantitative=$quantitative    
    """    

}


process PRINT_PATH {
  debug true
  output:
    stdout
  script:
  """
  which Rscript
  """
}

workflow {

    // Loads the parameters
    split_by_lc = params.split_by_lc
    quantitative = params.quantitative
    input_folder = params.input_folder
    metadata_file = params.metadata_file
    metadata_cols = params.metadata_cols
    

    // Splits the files into cycles
    cycles = channel.from( params.cycles)
    
    divided_file_paths = split_files(cycles, input_folder, metadata_file, split_by_lc)

    graph_paths = build_network(divided_file_paths.flatten(), metadata_file, metadata_cols, quantitative)

    //graph_paths.view { it }


}


// workflow {
//   PRINT_PATH()
// }



// On completion
workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Duration              : ${workflow.duration}"
    println "Execution status      : ${workflow.success ? 'All done!' : 'Failed' }"
}

// On error
workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

