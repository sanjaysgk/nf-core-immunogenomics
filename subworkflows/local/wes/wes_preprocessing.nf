/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WES PREPROCESSING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                              } from '../../modules/nf-core/fastqc/main'
include { BWA_MEM                             } from '../../modules/local/bwa/mem'
include { GATK4_SORTSAM                       } from '../../modules/local/gatk4/sortsam'
include { GATK4_MARKDUPLICATES                } from '../../modules/local/gatk4/markduplicates'
include { GATK4_COLLECTALIGNMENTSUMMARYMETRICS } from '../../modules/local/gatk4/collectalignmentsummarymetrics'
include { GATK4_BASERECALIBRATOR              } from '../../modules/local/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR                     } from '../../modules/local/gatk4/applybqsr'
include { GATK4_ANALYZECOVARIATES             } from '../../modules/local/gatk4/analyzecovariates'

workflow WES_PREPROCESSING {

    take:
    ch_samples           // channel: [ val(meta), [ path(reads) ] ]
    ch_fasta            // channel: [ path(fasta) ]
    ch_fai              // channel: [ path(fai) ]
    ch_dict             // channel: [ path(dict) ]
    ch_bwa_index        // channel: [ path(bwa_index) ]
    ch_dbsnp            // channel: [ path(dbsnp) ]
    ch_dbsnp_tbi        // channel: [ path(dbsnp_tbi) ]
    ch_mills            // channel: [ path(mills) ]
    ch_mills_tbi        // channel: [ path(mills_tbi) ]
    ch_indels           // channel: [ path(indels) ]
    ch_indels_tbi       // channel: [ path(indels_tbi) ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: FastQC for raw reads
    //
    FASTQC (
        ch_samples
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: BWA-MEM2 alignment
    //
    BWA_MEM (
        ch_samples,
        ch_bwa_index,
        ch_fasta,
        true  // sort output
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Sort SAM files
    //
    GATK4_SORTSAM (
        BWA_MEM.out.bam,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_SORTSAM.out.versions.first())

    //
    // MODULE: Mark duplicates
    //
    GATK4_MARKDUPLICATES (
        GATK4_SORTSAM.out.bam,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    //
    // MODULE: Collect alignment summary metrics
    //
    GATK4_COLLECTALIGNMENTSUMMARYMETRICS (
        GATK4_MARKDUPLICATES.out.bam,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_COLLECTALIGNMENTSUMMARYMETRICS.out.versions.first())

    //
    // MODULE: Base Quality Score Recalibration (BQSR) - First pass
    //
    ch_known_sites = ch_dbsnp.combine(ch_dbsnp_tbi)
                           .combine(ch_mills.combine(ch_mills_tbi))
                           .combine(ch_indels.combine(ch_indels_tbi))

    GATK4_BASERECALIBRATOR (
        GATK4_MARKDUPLICATES.out.bam,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    //
    // MODULE: Apply BQSR
    //
    GATK4_APPLYBQSR (
        GATK4_MARKDUPLICATES.out.bam.join(GATK4_BASERECALIBRATOR.out.table),
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

    //
    // MODULE: Base Quality Score Recalibration (BQSR) - Second pass for QC
    //
    GATK4_BASERECALIBRATOR (
        GATK4_APPLYBQSR.out.bam,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites
    )

    //
    // MODULE: Analyze covariates for BQSR QC plots
    //
    GATK4_ANALYZECOVARIATES (
        GATK4_BASERECALIBRATOR.out.table.join(GATK4_BASERECALIBRATOR.out.table)  // before and after tables
    )
    ch_versions = ch_versions.mix(GATK4_ANALYZECOVARIATES.out.versions.first())

    emit:
    recalibrated_bam     = GATK4_APPLYBQSR.out.bam        // channel: [ val(meta), path(bam), path(bai) ]
    fastqc_html          = FASTQC.out.html                // channel: [ val(meta), path(html) ]
    fastqc_zip           = FASTQC.out.zip                 // channel: [ val(meta), path(zip) ]
    duplicate_metrics    = GATK4_MARKDUPLICATES.out.metrics // channel: [ val(meta), path(metrics) ]
    alignment_metrics    = GATK4_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics // channel: [ val(meta), path(metrics) ]
    bqsr_table          = GATK4_BASERECALIBRATOR.out.table // channel: [ val(meta), path(table) ]
    bqsr_plots          = GATK4_ANALYZECOVARIATES.out.pdf  // channel: [ val(meta), path(pdf) ]
    versions            = ch_versions                      // channel: [ path(versions.yml) ]
}