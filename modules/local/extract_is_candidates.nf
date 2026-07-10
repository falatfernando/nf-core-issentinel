process EXTRACT_IS_CANDIDATES {
    tag "Annotating $sample_id"
    label 'process_single'
    
    // 1. Strict, version-pinned dependencies. 
    // nf-core requires exact versions to guarantee reproducibility.
    conda "conda-forge::pandas=2.1.0 bioconda::biopython=1.81"
    
    // 2. The standard nf-core container directive block.
    // If you enable Wave (-with-wave), Nextflow ignores this and builds a custom image from the Conda block above.
    // Otherwise, it falls back to a pre-built 'mulled' container that has both pandas and biopython.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0bc25c1566cf2c2bc016a6d6168ceb26090a2a5e:275217823f95f4ad6111fdbd783d085df7ffb2a0-0' :
        'quay.io/biocontainers/mulled-v2-0bc25c1566cf2c2bc016a6d6168ceb26090a2a5e:275217823f95f4ad6111fdbd783d085df7ffb2a0-0' }"

    input:
    tuple val(sample_id), path(vcf)
    path gbk_file

    output:
    tuple val(sample_id), path("${sample_id}_is_candidates.csv"), emit: candidates

    script:
    """
    extract_is_candidates.py \\
        --vcf ${vcf} \\
        --gbk ${gbk_file} \\
        --out ${sample_id}_is_candidates.csv \\
        --max-dist ${params.max_dist} \\
        --min-pe-iso ${params.min_pe_iso} \\
        --min-sr-iso ${params.min_sr_iso} \\
        --min-total-clust ${params.min_total_clust} \\
        --decoy ${params.decoy_name}
    """
}
