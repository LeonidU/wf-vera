//
// Read-based variant calling against reference genomes.
//
//   1. per sample, choose references (MMseqs2 top hits of the viral contigs, or a fixed
//      FASTA panel) and collapse near-identical ones
//   2. map the cleaned reads competitively to all of the sample's references (bowtie2)
//   3. keep primary reads with MAPQ >= variant_min_mapq, measure coverage/callable bases
//   4. call variants with LoFreq or iVar (params.variant_caller)
//   5. summarise into caller-independent tables + a compact JSON for the LLM report
//

include { PREPARE_VARIANT_REFERENCES                } from '../../modules/local/variants/prepare_references/main'
include { BOWTIE2_BUILD as VARIANT_BOWTIE2_BUILD    } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as VARIANT_BOWTIE2_ALIGN    } from '../../modules/nf-core/bowtie2/align/main'
include { VARIANT_BAM_QC                            } from '../../modules/local/variants/bam_qc/main'
include { LOFREQ_CALL                               } from '../../modules/local/variants/lofreq_call/main'
include { IVAR_VARIANTS                             } from '../../modules/local/variants/ivar_variants/main'
include { SUMMARIZE_VARIANTS                        } from '../../modules/local/variants/summarize/main'

workflow VARIANT_CALLING {

    take:
    ch_reads   // [ meta, reads ]          cleaned reads (CAT_FASTQ output)
    ch_tophit  // [ meta, tophit_aln ]     MMSEQS_NUCLEOTIDES top hits (unused in fasta_all)
    ch_taxa    // [ meta, filtered_taxa ]  PARSE_MMSEQS_NUCLEOTIDE_RESULTS table

    main:
    ch_versions = Channel.empty()

    def caller = params.variant_caller.toString().toLowerCase()
    def mode   = params.variant_ref_source.toString().toLowerCase()
    if ( !(caller in ['lofreq', 'ivar']) ) {
        error "variant_caller must be 'lofreq' or 'ivar', got '${params.variant_caller}'"
    }
    if ( !(mode in ['mmseqs', 'fasta', 'fasta_all']) ) {
        error "variant_ref_source must be 'mmseqs', 'fasta' or 'fasta_all', got '${params.variant_ref_source}'"
    }
    if ( !(params.sequencing_platform in ['illumina', 'iontorrent']) ) {
        error "sequencing_platform must be 'illumina' or 'iontorrent', got '${params.sequencing_platform}'"
    }
    if ( mode != 'mmseqs' && !params.variant_ref_fasta ) {
        error "variant_ref_source='${mode}' needs variant_ref_fasta"
    }

    ref_source = mode == 'mmseqs' ?
        file(params.mmseqs_easy_taxonomy_refseq_viral_db, checkIfExists: true) :
        file(params.variant_ref_fasta, checkIfExists: true)

    //
    // 1. references per sample
    //
    if ( mode == 'fasta_all' ) {
        ch_prepare_in = ch_reads.map { meta, reads -> [ meta, [], [] ] }
    } else {
        ch_prepare_in = ch_tophit.join(ch_taxa)
    }

    PREPARE_VARIANT_REFERENCES ( ch_prepare_in, ref_source, mode )
    ch_versions = ch_versions.mix(PREPARE_VARIANT_REFERENCES.out.versions.first())

    ch_refs = PREPARE_VARIANT_REFERENCES.out.fasta

    //
    // 2. competitive mapping to all references of the sample
    //
    VARIANT_BOWTIE2_BUILD ( ch_refs )
    ch_versions = ch_versions.mix(VARIANT_BOWTIE2_BUILD.out.versions.first())

    ch_map_in = ch_reads
        .join(VARIANT_BOWTIE2_BUILD.out.index)
        .join(ch_refs)
        .multiMap { meta, reads, index, fasta ->
            reads: [ meta, reads ]
            index: [ meta, index ]
            fasta: [ meta, fasta ]
        }

    VARIANT_BOWTIE2_ALIGN (
        ch_map_in.reads,
        ch_map_in.index,
        ch_map_in.fasta,
        false,   // save_unaligned
        true     // sort_bam
    )
    ch_versions = ch_versions.mix(VARIANT_BOWTIE2_ALIGN.out.versions.first())

    //
    // 3. read filter + coverage / callable bases
    //
    VARIANT_BAM_QC ( VARIANT_BOWTIE2_ALIGN.out.bam.join(ch_refs) )
    ch_versions = ch_versions.mix(VARIANT_BAM_QC.out.versions.first())

    ch_call_in = VARIANT_BAM_QC.out.bam.join(ch_refs)   // [ meta, bam, bai, fasta ]

    //
    // 4. variant calling
    //
    if ( caller == 'lofreq' ) {
        LOFREQ_CALL ( ch_call_in )
        ch_calls    = LOFREQ_CALL.out.vcf
        ch_versions = ch_versions.mix(LOFREQ_CALL.out.versions.first())
    } else {
        IVAR_VARIANTS ( ch_call_in )
        ch_calls    = IVAR_VARIANTS.out.tsv
        ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions.first())
    }

    //
    // 5. tables for the report
    //
    ch_summary_in = PREPARE_VARIANT_REFERENCES.out.table
        .join(ch_refs)
        .join(VARIANT_BAM_QC.out.coverage)
        .join(VARIANT_BAM_QC.out.callable)
        .join(ch_calls)

    SUMMARIZE_VARIANTS ( ch_summary_in, caller )
    ch_versions = ch_versions.mix(SUMMARIZE_VARIANTS.out.versions.first())

    emit:
    variants   = SUMMARIZE_VARIANTS.out.variants   // [ meta, *_variants.tsv ]
    summary    = SUMMARIZE_VARIANTS.out.summary    // [ meta, *_variant_summary.tsv ]
    json       = SUMMARIZE_VARIANTS.out.json       // [ meta, *_variants.json ]
    references = PREPARE_VARIANT_REFERENCES.out.table
    calls      = ch_calls
    bowtie2_log = VARIANT_BOWTIE2_ALIGN.out.log
    versions   = ch_versions
}
