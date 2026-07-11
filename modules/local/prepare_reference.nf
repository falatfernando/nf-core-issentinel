process PREPARE_REFERENCE {
    tag "Engineering Reference"
    label 'process_single' // Uses minimal CPU/RAM from base.config

    conda "conda-forge::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path gbk
    path fasta
    path is6110_fasta

    output:
    // We name the output channel 'masked_fasta' so we can grab it easily
    path "masked_reference.fasta", emit: masked_fasta

    script:
    """
    # Calls your python script from the bin/ folder
    prepare_masked_reference.py \\
        --gbk ${gbk} \\
        --fasta ${fasta} \\
        --is6110 ${is6110_fasta} \\
        --out masked_reference.fasta
    """
}
