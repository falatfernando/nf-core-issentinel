process BWA_ALIGN {
    tag "Aligning $sample_id"
    label 'process_high'
    
    // 'mulled' container that has BOTH bwa and samtools installed
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"

    input:
    tuple val(sample_id), path(reads)
    path ref_fasta
    path ref_indexes // We pass the .bwt, .fai, etc. down the channel so it doesn't re-index!

    output:
    tuple val(sample_id), path("${sample_id.id}.sorted.bam"), path("${sample_id.id}.sorted.bam.bai"), emit: bam

    script:
    """
    # \$task.cpus to dynamically scale threads
    bwa mem -t ${task.cpus} ${ref_fasta} ${reads[0]} ${reads[1]} | \\
    samtools view -Sb - | \\
    samtools sort -@ ${task.cpus} -o ${sample_id.id}.sorted.bam -
    
    samtools index ${sample_id.id}.sorted.bam
    """
}
