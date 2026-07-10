process FILTER_REGION {
    tag "Filtering mmpL5-Rv0678"
    label 'process_single'
    
    // We reuse the exact same container from COMPILE_RESULTS since it already has Pandas
    conda "conda-forge::pandas=2.1.0 bioconda::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bc25c1566cf2c2bc016a6d6168ceb26090a2a5e:275217823f95f4ad6111fdbd783d085df7ffb2a0-0' :
        'quay.io/biocontainers/mulled-v2-0bc25c1566cf2c2bc016a6d6168ceb26090a2a5e:275217823f95f4ad6111fdbd783d085df7ffb2a0-0' }"

    input:
    path compiled_csv // This will accept the 'raw_table' emitted by COMPILE_RESULTS

    output:
    path "is_sentinel_roi.csv"           , emit: raw_roi
    path "is_sentinel_roi_formatted.csv" , emit: formatted_roi

    script:
    """
    filter_region.py \\
        --input ${compiled_csv} \\
        --raw-out is_sentinel_roi.csv \\
        --formatted-out is_sentinel_roi_formatted.csv \\
        --max-dist ${params.max_dist}
    """
}
