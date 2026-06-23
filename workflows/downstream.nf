/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_TAPS_CONVERSION as METHYLSEQ } from '../subworkflows/nf-core/bam_taps_conversion/main'

// TODO: include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'
// TODO: include { TVC                   } from '../modules/local/tvc/main'
// TODO: include { MULTIQC               } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DOWNSTREAM ANALYSIS WORKFLOW
    Receives the final consensus-filtered BAM from FASTQUORUM and runs:
      - METHYLSEQ  : rastair methylation calling (TAPS)
      - GATK HaplotypeCaller : variant calling            (TODO)
      - TVC (Watchmaker)     : variant calling            (TODO)
      - MultiQC              : aggregate QC report        (TODO - run last)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNSTREAM {

    take:
    ch_bam      // channel: [ val(meta), path(bam) ]
    ch_bai      // channel: [ val(meta), path(bai) ]
    ch_fasta    // channel: [ val(meta), path(fasta) ]  (value channel)
    ch_fasta_fai // channel: [ val(meta), path(fai) ]   (value channel)

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Build broadcast inputs: combine per-sample bam/bai with reference channels
    //
    ch_taps_inputs = ch_bam
        .join(ch_bai)
        .combine(ch_fasta)
        .combine(ch_fasta_fai)
        .multiMap { meta, bam, bai, _meta_fasta, fasta, _meta_fai, fai ->
            bam:         [ meta, bam ]
            bai:         [ meta, bai ]
            fasta:       [ meta, fasta ]
            fasta_index: [ meta, fai ]
        }

    //
    // SUBWORKFLOW: TAPS methylation calling with Rastair
    //   rastair mbias  -> per-position C->T bias
    //   rastair call   -> genome-wide methylation calls
    //   rastair methylkit -> methylKit-format output
    //
    METHYLSEQ(
        ch_taps_inputs.bam,
        ch_taps_inputs.bai,
        ch_taps_inputs.fasta,
        ch_taps_inputs.fasta_index
    )
    ch_versions = ch_versions.mix(METHYLSEQ.out.versions)

    // TODO: GATK4_HAPLOTYPECALLER — variant calling in parallel with rastair
    // TODO: TVC (Watchmaker)      — variant calling in parallel with rastair

    //
    // Collate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'downstream_software_mqc_versions.yml',
            sort: true,
            newLine: true,
        )

    // TODO: MultiQC — collect metrics from all steps above and run at the end
    // MULTIQC(
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     [],
    //     [],
    //     [],
    //     [],
    // )

    emit:
    rastair_mbias     = METHYLSEQ.out.mbias     // channel: [ val(meta), path(*.txt)    ]
    rastair_call      = METHYLSEQ.out.call      // channel: [ val(meta), path(*.txt)    ]
    rastair_methylkit = METHYLSEQ.out.methylkit // channel: [ val(meta), path(*.txt.gz) ]
    versions          = ch_versions
}
