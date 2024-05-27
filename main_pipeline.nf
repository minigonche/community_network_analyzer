#!/usr/bin/env nextflow

process pca_impute {
    input:
    path metadata_file


    output:
    path "*.csv"

    script:
    """            
    pca_impute.R --metadata_file=$metadata_file    
    """    

}

process split_files {
    input:
    val cycle
    path input_folder
    path metadata_file


    output:
    path "*.csv"

    script:
    """            
    file_splitter.R --cycle=$cycle --input_folder=$input_folder --metadata_file=$metadata_file    
    """    

}

process build_network {
    input:
    path input_file
    path metadata_file


    output:
    path "*.graphml"

    script:
    """            
    network_builder.R --input_file=$input_file --metadata_file=$metadata_file    
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
    input_folder = params.input_folder
    metadata_file = params.metadata_file

    metadata_file_imputed = pca_impute(metadata_file)

    // Splits the files into cycles
    cycles = channel.from( params.cycles)
    divided_file_paths = split_files(cycles, input_folder, metadata_file_imputed)

    graph_paths = build_network(divided_file_paths.flatten(), metadata_file_imputed)

    graph_paths.view { it }


}


// workflow {
//   PRINT_PATH()
// }
