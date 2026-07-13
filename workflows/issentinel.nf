/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// nf-core utilities
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_issentinel_pipeline'

// IS Sentinel custom components
include { PREPARE_REFERENCE      } from '../modules/local/prepare_reference'
include { INDEX_REFERENCE        } from '../modules/local/index_reference'
include { DETECT_IS              } from '../subworkflows/local/detect_is'
include { COMPILE_RESULTS        } from '../modules/local/compile_results'
include { FILTER_REGION          } from '../modules/local/filter_region'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISSENTINEL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input (acts as our ch_reads)
    outdir         // the output directory path

    main:

    def ch_versions = channel.empty()

    // ==========================================================
    // 1. Initialize static reference channels from the params
    // ==========================================================
    ch_gbk       = Channel.fromPath(params.gbk_file).first()
    ch_is_fasta  = Channel.fromPath(params.is_fasta).first()
    ch_raw_ref   = Channel.fromPath(params.ref_genome).first()

    // ==========================================================
    // STEP 0: PREPARE REFERENCE
    // ==========================================================
    PREPARE_REFERENCE(ch_gbk, ch_raw_ref, ch_is_fasta)
    // The active reference is now the newly masked output containing the decoy contig
    ch_active_ref = PREPARE_REFERENCE.out.masked_fasta


    // Always index the active reference before alignment
    INDEX_REFERENCE(ch_active_ref)
    ch_ref_indexes = INDEX_REFERENCE.out.indexes

    // ==========================================================
    // STEP 1: DETECT IS (PARALLEL SUBWORKFLOW)
    // ==========================================================
    ch_candidates = Channel.empty()
    
    if (params.run_pipeline) {
        // This spins up N containers in parallel (one for each sample)
        DETECT_IS(ch_samplesheet, ch_active_ref, ch_ref_indexes, ch_gbk)
        
        // Isolate just the CSV files from the output tuple
        ch_candidates = DETECT_IS.out.candidates_csv.map { meta, csv -> csv }
    }

    // ==========================================================
    // STEP 2: COMPILE RESULTS
    // ==========================================================
    ch_compiled_raw = Channel.empty()

    if (params.run_compile) {
        // .collect() waits for all parallel tasks to finish and bundles the CSVs into one list!
        ch_all_csvs = ch_candidates.collect()
        
        COMPILE_RESULTS(ch_all_csvs)
        ch_compiled_raw = COMPILE_RESULTS.out.raw_table
    }

    // ==========================================================
    // STEP 3: FILTER REGION (mmpL5-Rv0678)
    // ==========================================================
    if (params.run_filter) {
        // Pipes the compiled master table directly into the filter
        FILTER_REGION(ch_compiled_raw)
    }

    // ==========================================================
    // nf-core: Collate and save software versions
    // ==========================================================
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'issentinel_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
        
    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
