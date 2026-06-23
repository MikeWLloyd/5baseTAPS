/*
 * TAPS methylation conversion subworkflow (rastair 2.1.1)
 *
 * rastair call     -> per-position BED (+ VCF) from BAM     ─> methylKit
 * rastair per-read -> per-read BED from BAM (parallel)      ─> tabix ─> mbias HTML
 */

include { RASTAIR_CALL         } from '../../../modules/nf-core/rastair/call/main'
include { METH_SUMMARY } from '../../../modules/local/rastair_call_summary/main'
include { PERREAD      } from '../../../modules/nf-core/rastair/perread/main'
include { TABIX_TABIX          } from '../../../modules/local/tabix_tabix/main'
include { PERCHROMOSOME_MBIAS as MBIAS } from '../../local/perChr_mbias/main'
include { METHYLKIT    } from '../../../modules/nf-core/rastair/methylkit/main'

workflow BAM_TAPS_CONVERSION {

    take:
    ch_bam         // channel: [ val(meta), path(bam) ]
    ch_bai         // channel: [ val(meta), path(bai) ]
    ch_fasta       // channel: [ val(meta), path(fa)  ]
    ch_fasta_index // channel: [ val(meta), path(fai) ]

    main:
    ch_versions = Channel.empty()

    log.info "Running TAPS conversion with Rastair 2.1.1: call + per-read -> mbias + methylKit"

    //
    // STEP 1a: call methylated positions; emit per-position BED and VCF
    //
    RASTAIR_CALL(
        ch_bam,
        ch_bai,
        ch_fasta,
        ch_fasta_index,
    )
    ch_versions = ch_versions.mix(RASTAIR_CALL.out.versions)

    //
    // STEP 1b: per-read methylation calls for M-bias analysis (runs in parallel with CALL)
    //
    PERREAD(
        ch_bam,
        ch_bai,
        ch_fasta,
        ch_fasta_index,
    )
    ch_versions = ch_versions.mix(PERREAD.out.versions)

    //
    // STEP 2a: tabix-index the per-read BED (rastair container lacks tabix)
    //
    TABIX_TABIX(
        PERREAD.out.bed
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    //
    // STEP 2b: M-bias HTML report — per-chromosome two-stage strategy
    //          (avoids R 2^31 vector-length crash on large samples)
    //
    MBIAS(
        PERREAD.out.bed.join(TABIX_TABIX.out.tbi),
        RASTAIR_CALL.out.vcf,
        ch_fasta
    )

    //
    // STEP 2c: methylKit-format table from per-position BED (pure shell, no container)
    //
    METHYLKIT(
        RASTAIR_CALL.out.bed
    )

    //
    // STEP 2d: global CpG methylation summary from per-position BED
    //
    METH_SUMMARY(
        RASTAIR_CALL.out.bed
    )
    ch_versions = ch_versions.mix(METH_SUMMARY.out.versions)

    emit:
    mbias            = MBIAS.out.html                   // channel: [ val(meta), path("*.html")    ]
    bed              = RASTAIR_CALL.out.bed                     // channel: [ val(meta), path("*.bed.gz")  ]
    vcf              = RASTAIR_CALL.out.vcf                     // channel: [ val(meta), path("*.vcf.gz")  ]
    perread          = PERREAD.out.bed                  // channel: [ val(meta), path("*.bed.gz")  ]
    methylkit        = METHYLKIT.out.methylkit          // channel: [ val(meta), path("*.txt.gz")  ]
    meth_summary     = METH_SUMMARY.out.tsv             // channel: [ val(meta), path("*.tsv")     ]
    versions         = ch_versions
}
