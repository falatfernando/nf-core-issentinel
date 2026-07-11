process DELLY_CALL {
    tag "SV Calling $sample_id"
    label 'process_medium'
    
    conda "bioconda::delly=1.1.6 bioconda::bcftools=1.15.1 bioconda::samtools=1.15.1"
    container "quay.io/biocontainers/mulled-v2-f1f590052697e4ef32c9e2fa8376efc72f366d05:9b23d11a789c7a78a0cc4a6bb7f41103097dd04f-0"

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref_fasta
    path ref_indexes 

    output:
    tuple val(sample_id), path("${sample_id.id}.vcf"), emit: vcf

    script:
    """
    delly call -q 0 -g ${ref_fasta} -o ${sample_id.id}.bcf ${bam}
    bcftools view ${sample_id.id}.bcf > ${sample_id.id}.vcf
    """
}
