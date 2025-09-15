/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WES VARIANT CALLING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_MUTECT2                  } from '../../modules/local/gatk4/mutect2'
include { GATK4_LEARNREADORIENTATIONMODEL } from '../../modules/local/gatk4/learnreadorientationmodel'
include { GATK4_GETPILEUPSUMMARIES       } from '../../modules/local/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION   } from '../../modules/local/gatk4/calculatecontamination'
include { GATK4_FILTERMUTECTCALLS        } from '../../modules/local/gatk4/filtermutectcalls'
include { GATK4_SELECTVARIANTS           } from '../../modules/local/gatk4/selectvariants'

workflow WES_VARIANT_CALLING {

    take:
    ch_tumor_normal_pairs    // channel: [ val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai) ]
    ch_fasta                // channel: [ path(fasta) ]
    ch_fai                  // channel: [ path(fai) ]
    ch_dict                 // channel: [ path(dict) ]
    ch_pon                  // channel: [ path(pon) ]
    ch_pon_tbi              // channel: [ path(pon_tbi) ]
    ch_germline_resource    // channel: [ path(germline_resource) ]
    ch_germline_resource_tbi // channel: [ path(germline_resource_tbi) ]
    ch_small_exac           // channel: [ path(small_exac) ]
    ch_small_exac_tbi       // channel: [ path(small_exac_tbi) ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Mutect2 somatic variant calling
    //
    GATK4_MUTECT2 (
        ch_tumor_normal_pairs,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_pon,
        ch_pon_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    //
    // MODULE: Learn read orientation model
    //
    GATK4_LEARNREADORIENTATIONMODEL (
        GATK4_MUTECT2.out.f1r2
    )
    ch_versions = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions.first())

    //
    // MODULE: Get pileup summaries for tumor samples
    //
    ch_tumor_bams = ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
        def tumor_meta = [
            sample: meta.tumor_sample,
            patient_id: meta.patient_id,
            sample_type: 'tumor'
        ]
        return [tumor_meta, tumor_bam, tumor_bai]
    }

    GATK4_GETPILEUPSUMMARIES (
        ch_tumor_bams,
        ch_small_exac,
        ch_small_exac_tbi
    )
    ch_tumor_pileups = GATK4_GETPILEUPSUMMARIES.out.table

    //
    // MODULE: Get pileup summaries for normal samples
    //
    ch_normal_bams = ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
        def normal_meta = [
            sample: meta.normal_sample,
            patient_id: meta.patient_id,
            sample_type: 'normal'
        ]
        return [normal_meta, normal_bam, normal_bai]
    }

    GATK4_GETPILEUPSUMMARIES (
        ch_normal_bams,
        ch_small_exac,
        ch_small_exac_tbi
    )
    ch_normal_pileups = GATK4_GETPILEUPSUMMARIES.out.table

    //
    // MODULE: Calculate contamination
    //
    ch_contamination_input = ch_tumor_pileups
        .join(ch_normal_pileups, by: [0])
        .map { meta, tumor_table, normal_table ->
            def pair_meta = [
                patient_id: meta.patient_id,
                tumor_sample: meta.sample_type == 'tumor' ? meta.sample : null,
                normal_sample: meta.sample_type == 'normal' ? meta.sample : null
            ]
            return [pair_meta, tumor_table, normal_table]
        }

    GATK4_CALCULATECONTAMINATION (
        ch_contamination_input
    )
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions.first())

    //
    // MODULE: Filter Mutect2 calls
    //
    ch_filter_input = GATK4_MUTECT2.out.vcf
        .join(GATK4_MUTECT2.out.stats)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior)
        .join(GATK4_CALCULATECONTAMINATION.out.contamination)
        .join(GATK4_CALCULATECONTAMINATION.out.segmentation)

    GATK4_FILTERMUTECTCALLS (
        ch_filter_input,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())

    //
    // MODULE: Select high-confidence variants
    //
    GATK4_SELECTVARIANTS (
        GATK4_FILTERMUTECTCALLS.out.vcf,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())

    emit:
    vcf             = GATK4_MUTECT2.out.vcf           // channel: [ val(meta), path(vcf) ]
    filtered_vcf    = GATK4_FILTERMUTECTCALLS.out.vcf // channel: [ val(meta), path(vcf) ]
    pass_vcf        = GATK4_SELECTVARIANTS.out.vcf    // channel: [ val(meta), path(vcf) ]
    contamination   = GATK4_CALCULATECONTAMINATION.out.contamination // channel: [ val(meta), path(contamination) ]
    versions        = ch_versions                     // channel: [ path(versions.yml) ]
}