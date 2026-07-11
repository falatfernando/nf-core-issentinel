process FILTER_REGION {
    tag "Filtering ${params.filter_name}"
    label 'process_single'
    
    conda "conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    path compiled_csv // This will accept the 'raw_table' emitted by COMPILE_RESULTS

    output:
    path "is_sentinel_${params.filter_name}.csv"           , emit: raw_roi
    path "is_sentinel_${params.filter_name}_formatted.csv" , emit: formatted_roi

    script:
    """
    filter_region.py \\
        --input ${compiled_csv} \\
        --raw-out is_sentinel_${params.filter_name}.csv \\
        --formatted-out is_sentinel_${params.filter_name}_formatted.csv \\
        --chrom ${params.filter_chrom} \\
        --start ${params.filter_start} \\
        --end ${params.filter_end} \\
        --region-name ${params.filter_name} \\
        --max-dist ${params.max_dist}
    """
}
