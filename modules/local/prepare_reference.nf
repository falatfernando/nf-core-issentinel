process PREPARE_REFERENCE {
    tag "Engineering Reference"
    label 'process_single' // Uses minimal CPU/RAM from base.config

    // Biocontainer that comes with Biopython pre-installed
    container "quay.io/biocontainers/biopython:1.79--py38h7f98852_2"

    input:
    path gbk
    path fasta
    path is6110_fasta

    output:
    // Output channel named 'masked_fasta' so we can grab it easily
    path "masked_reference.fasta", emit: masked_fasta

    script:
    """
    # Calls the python script from the bin/ folder
    prepare_masked_reference.py \\
        --gbk ${gbk} \\
        --fasta ${fasta} \\
        --is6110 ${is6110_fasta} \\
        --out masked_reference.fasta
    """
}
