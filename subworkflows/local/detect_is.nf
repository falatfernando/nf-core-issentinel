include { BWA_ALIGN           } from '../../modules/local/bwa_align'
include { DELLY_CALL          } from '../../modules/local/delly_call'
include { EXTRACT_IS_CANDIDATES } from '../../modules/local/extract_is_candidates'

workflow DETECT_IS {
    take:
    reads_ch      // Channel: tuple val(sample_id), path(reads)
    ref_fasta     // Value channel: path to masked FASTA
    ref_indexes   // Value channel: path to reference indexes (BWA + FAI)
    gbk_file      // Value channel: path to reference GenBank file

    main:
    // 1. Align the raw reads (outputs BAM)
    BWA_ALIGN(reads_ch, ref_fasta, ref_indexes)

    // 2. Call Structural Variants using the BAM from BWA (outputs VCF)
    DELLY_CALL(BWA_ALIGN.out.bam, ref_fasta, ref_indexes)

    // 3. Extract and annotate candidate IS insertions from the VCF (outputs CSV)
    EXTRACT_IS_CANDIDATES(DELLY_CALL.out.vcf, gbk_file)

    emit:
    // Expose the final annotated CSVs so the main workflow can collect them
    candidates_csv = EXTRACT_IS_CANDIDATES.out.candidates
}
