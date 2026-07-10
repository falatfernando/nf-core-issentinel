process DELLY_CALL {
    tag "SV Calling $sample_id"
    label 'process_medium'
    
    container "quay.io/biocontainers/delly:1.1.6--h270b39a_0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_indexes 

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    """
    # Delly natively outputs VCF if you specify the .vcf extension, saving us the bcftools step!
    delly call -q 0 -g ${ref_fasta} -o ${sample_id}.vcf ${bam}
    """
}
