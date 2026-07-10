process COMPILE_RESULTS {
    tag "Compiling Master Table"
    label 'process_single'
    
    // We reuse the exact same container from ANNOTATE_DELLY because it already has Pandas!
    conda "conda-forge::pandas=2.1.0 bioconda::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bc25c1566cf2c2bc016a6d6168ceb26090a2a5e:275217823f95f4ad6111fdbd783d085df7ffb2a0-0' :
        'quay.io/biocontainers/mulled-v2-0bc25c1566cf2c2bc016a6d6168ceb26090a2a5e:275217823f95f4ad6111fdbd783d085df7ffb2a0-0' }"

    // IMPORTANT: 'path' here accepts a LIST of hundreds of files, not a single file!
    input:
    path csv_files

    output:
    path "is_sentinel_results.csv.gz"         , emit: raw_table
    path "is_sentinel_results_formatted.csv"  , emit: formatted_table

    script:
    """
    # 1. Trick your python script's glob pattern!
    # We create a dummy folder that matches the "results/*/*.csv" expected pattern
    mkdir -p dummy_results/sample_folder/
    mv *_is6110_candidates.csv dummy_results/sample_folder/

    # 2. Run your unchanged Python script
    compile_results.py \\
        --results dummy_results \\
        --raw-out is_sentinel_results.csv \\
        --formatted-out is_sentinel_results_formatted.csv \\
        --max-dist ${params.max_dist}
    """
}
