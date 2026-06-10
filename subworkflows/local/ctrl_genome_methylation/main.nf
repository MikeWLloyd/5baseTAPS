/*
 * Control genome methylation subworkflow (rd mode, lambda + pUC19)
 *
 * Runs the full fgbio UMI-collapse pipeline against a spike-in control genome,
 * then calls methylation with rastair and emits a global methylation summary.
 * Mirrors the FASTQUORUM rd-mode path exactly (CALL* → ALIGN → FILTER → RASTAIR).
 *
 * Two identically-structured workflows with distinct module aliases to avoid
 * Nextflow DSL2 process-name collisions when both are called from the same scope:
 *   LAMBDA_GENOME_METHYLATION  — lambda phage (negative control)
 *   PUC19_GENOME_METHYLATION   — pUC19 plasmid (positive control)
 *
 * Only --mode rd is supported; non-rd (CALLANDFILTER*) paths are not wired.
 */

// ── Lambda-specific module imports ─────────────────────────────────────────
include { ALIGN_BAM as LAMBDA_ALIGN_RAW_BAM                               } from '../../../modules/local/align_bam/main'
include { ALIGN_BAM as LAMBDA_ALIGN_CONSENSUS_BAM                         } from '../../../modules/local/align_bam/main'
include { SAMTOOLS_MERGE as LAMBDA_MERGE_BAM                               } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FLAGSTAT as LAMBDA_SAMTOOLS_FLAGSTAT                   } from '../../../modules/local/samtools/flagstat/main'
include { FGBIO_GROUPREADSBYUMI as LAMBDA_GROUPREADSBYUMI                 } from '../../../modules/local/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLDDUPLEXCONSENSUSREADS as LAMBDA_CALLDDUPLEXCONSENSUSREADS   } from '../../../modules/local/fgbio/callduplexconsensusreads/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as LAMBDA_CALLMOLECULARCONSENSUSREADS } from '../../../modules/local/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS as LAMBDA_FILTERCONSENSUSREADS       } from '../../../modules/local/fgbio/filterconsensusreads/main'
include { RASTAIR_CALL as LAMBDA_RASTAIR_CALL                             } from '../../../modules/nf-core/rastair/call/main'
include { RASTAIR_CALL_SUMMARY as LAMBDA_RASTAIR_CALL_SUMMARY             } from '../../../modules/local/rastair_call_summary/main'

// ── pUC19-specific module imports ──────────────────────────────────────────
include { ALIGN_BAM as PUC19_ALIGN_RAW_BAM                               } from '../../../modules/local/align_bam/main'
include { ALIGN_BAM as PUC19_ALIGN_CONSENSUS_BAM                         } from '../../../modules/local/align_bam/main'
include { SAMTOOLS_MERGE as PUC19_MERGE_BAM                               } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FLAGSTAT as PUC19_SAMTOOLS_FLAGSTAT                    } from '../../../modules/local/samtools/flagstat/main'
include { FGBIO_GROUPREADSBYUMI as PUC19_GROUPREADSBYUMI                 } from '../../../modules/local/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLDDUPLEXCONSENSUSREADS as PUC19_CALLDDUPLEXCONSENSUSREADS   } from '../../../modules/local/fgbio/callduplexconsensusreads/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as PUC19_CALLMOLECULARCONSENSUSREADS } from '../../../modules/local/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS as PUC19_FILTERCONSENSUSREADS       } from '../../../modules/local/fgbio/filterconsensusreads/main'
include { RASTAIR_CALL as PUC19_RASTAIR_CALL                             } from '../../../modules/nf-core/rastair/call/main'
include { RASTAIR_CALL_SUMMARY as PUC19_RASTAIR_CALL_SUMMARY             } from '../../../modules/local/rastair_call_summary/main'

// ──────────────────────────────────────────────────────────────────────────

workflow LAMBDA_GENOME_METHYLATION {

    take:
    ch_unmapped_bam    // channel: [ val(meta), path(unmapped.bam) ] — one per samplesheet row
    ch_fasta           // value channel: [ val(meta), path(fasta) ]
    ch_fasta_fai       // value channel: [ val(meta), path(fai) ]
    ch_dict            // value channel: [ val(meta), path(dict) ]
    ch_bwamem2         // value channel: [ val(meta), path(bwamem2_dir) ]

    main:
    ch_versions = Channel.empty()

    LAMBDA_ALIGN_RAW_BAM(
        ch_unmapped_bam,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_bwamem2,
        "template-coordinate"
    )
    ch_versions = ch_versions.mix(LAMBDA_ALIGN_RAW_BAM.out.versions.first())

    bam_to_merge = LAMBDA_ALIGN_RAW_BAM.out.bam
        .map { meta, bam ->
            def meta_no_lane = meta.findAll { k, v -> k != 'lane' }
            [groupKey(meta_no_lane, meta_no_lane.n_samples), bam]
        }
        .groupTuple()
        .branch { meta, bam ->
            single:   meta.n_samples <= 1
                return [meta, bam[0]]
            multiple: meta.n_samples > 1
        }

    LAMBDA_MERGE_BAM(bam_to_merge.multiple, [[], []], [[], []])
    ch_versions = ch_versions.mix(LAMBDA_MERGE_BAM.out.versions.first())

    bam_all = LAMBDA_MERGE_BAM.out.bam.mix(bam_to_merge.single)

    // Run flagstat on merged BAM (lane-stripped meta) so ch_lambda_flagstat key
    // matches ch_prededup_flagstat for the TAPS_MAPPING_METRICS join.
    LAMBDA_SAMTOOLS_FLAGSTAT(bam_all)
    ch_versions = ch_versions.mix(LAMBDA_SAMTOOLS_FLAGSTAT.out.versions.first())

    def umi_strategy = params.groupreadsbyumi_strategy ?: (params.duplex_seq ? 'Paired' : 'Adjacency')
    log.info("[fgbio GroupReadsByUmi] strategy='${umi_strategy}' (duplex_seq=${params.duplex_seq})")
    LAMBDA_GROUPREADSBYUMI(bam_all, umi_strategy, params.groupreadsbyumi_edits)
    ch_versions = ch_versions.mix(LAMBDA_GROUPREADSBYUMI.out.versions.first())

    if (params.duplex_seq) {
        LAMBDA_CALLDDUPLEXCONSENSUSREADS(
            LAMBDA_GROUPREADSBYUMI.out.bam,
            params.call_min_reads,
            params.call_min_baseq
        )
        ch_versions = ch_versions.mix(LAMBDA_CALLDDUPLEXCONSENSUSREADS.out.versions.first())
        ch_consensus_bam = LAMBDA_CALLDDUPLEXCONSENSUSREADS.out.bam
    }
    else {
        LAMBDA_CALLMOLECULARCONSENSUSREADS(
            LAMBDA_GROUPREADSBYUMI.out.bam,
            params.call_min_reads,
            params.call_min_baseq
        )
        ch_versions = ch_versions.mix(LAMBDA_CALLMOLECULARCONSENSUSREADS.out.versions.first())
        ch_consensus_bam = LAMBDA_CALLMOLECULARCONSENSUSREADS.out.bam
    }

    LAMBDA_ALIGN_CONSENSUS_BAM(
        ch_consensus_bam,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_bwamem2,
        "none"
    )
    ch_versions = ch_versions.mix(LAMBDA_ALIGN_CONSENSUS_BAM.out.versions.first())

    LAMBDA_FILTERCONSENSUSREADS(
        LAMBDA_ALIGN_CONSENSUS_BAM.out.bam,
        ch_fasta,
        params.filter_min_reads,
        params.filter_min_baseq,
        params.filter_max_base_error_rate
    )
    ch_versions = ch_versions.mix(LAMBDA_FILTERCONSENSUSREADS.out.versions.first())

    LAMBDA_RASTAIR_CALL(
        LAMBDA_FILTERCONSENSUSREADS.out.bam,
        LAMBDA_FILTERCONSENSUSREADS.out.bai,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(LAMBDA_RASTAIR_CALL.out.versions)

    LAMBDA_RASTAIR_CALL_SUMMARY(
        LAMBDA_RASTAIR_CALL.out.bed
    )
    ch_versions = ch_versions.mix(LAMBDA_RASTAIR_CALL_SUMMARY.out.versions)

    emit:
    meth_summary  = LAMBDA_RASTAIR_CALL_SUMMARY.out.tsv       // channel: [ val(meta), path("*.methylation_summary.tsv") ]
    flagstat      = LAMBDA_SAMTOOLS_FLAGSTAT.out.flagstat      // channel: [ val(meta), path("*.flagstat") ]
    versions      = ch_versions
}

// ──────────────────────────────────────────────────────────────────────────

workflow PUC19_GENOME_METHYLATION {

    take:
    ch_unmapped_bam
    ch_fasta
    ch_fasta_fai
    ch_dict
    ch_bwamem2

    main:
    ch_versions = Channel.empty()

    PUC19_ALIGN_RAW_BAM(
        ch_unmapped_bam,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_bwamem2,
        "template-coordinate"
    )
    ch_versions = ch_versions.mix(PUC19_ALIGN_RAW_BAM.out.versions.first())

    bam_to_merge = PUC19_ALIGN_RAW_BAM.out.bam
        .map { meta, bam ->
            def meta_no_lane = meta.findAll { k, v -> k != 'lane' }
            [groupKey(meta_no_lane, meta_no_lane.n_samples), bam]
        }
        .groupTuple()
        .branch { meta, bam ->
            single:   meta.n_samples <= 1
                return [meta, bam[0]]
            multiple: meta.n_samples > 1
        }

    PUC19_MERGE_BAM(bam_to_merge.multiple, [[], []], [[], []])
    ch_versions = ch_versions.mix(PUC19_MERGE_BAM.out.versions.first())

    bam_all = PUC19_MERGE_BAM.out.bam.mix(bam_to_merge.single)

    // Run flagstat on merged BAM (lane-stripped meta) — matches ch_prededup_flagstat join key.
    PUC19_SAMTOOLS_FLAGSTAT(bam_all)
    ch_versions = ch_versions.mix(PUC19_SAMTOOLS_FLAGSTAT.out.versions.first())

    def puc19_umi_strategy = params.groupreadsbyumi_strategy ?: (params.duplex_seq ? 'Paired' : 'Adjacency')
    log.info("[fgbio GroupReadsByUmi] pUC19 strategy='${puc19_umi_strategy}' (duplex_seq=${params.duplex_seq})")
    PUC19_GROUPREADSBYUMI(bam_all, puc19_umi_strategy, params.groupreadsbyumi_edits)
    ch_versions = ch_versions.mix(PUC19_GROUPREADSBYUMI.out.versions.first())

    if (params.duplex_seq) {
        PUC19_CALLDDUPLEXCONSENSUSREADS(
            PUC19_GROUPREADSBYUMI.out.bam,
            params.call_min_reads,
            params.call_min_baseq
        )
        ch_versions = ch_versions.mix(PUC19_CALLDDUPLEXCONSENSUSREADS.out.versions.first())
        ch_consensus_bam = PUC19_CALLDDUPLEXCONSENSUSREADS.out.bam
    }
    else {
        PUC19_CALLMOLECULARCONSENSUSREADS(
            PUC19_GROUPREADSBYUMI.out.bam,
            params.call_min_reads,
            params.call_min_baseq
        )
        ch_versions = ch_versions.mix(PUC19_CALLMOLECULARCONSENSUSREADS.out.versions.first())
        ch_consensus_bam = PUC19_CALLMOLECULARCONSENSUSREADS.out.bam
    }

    PUC19_ALIGN_CONSENSUS_BAM(
        ch_consensus_bam,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_bwamem2,
        "none"
    )
    ch_versions = ch_versions.mix(PUC19_ALIGN_CONSENSUS_BAM.out.versions.first())

    PUC19_FILTERCONSENSUSREADS(
        PUC19_ALIGN_CONSENSUS_BAM.out.bam,
        ch_fasta,
        params.filter_min_reads,
        params.filter_min_baseq,
        params.filter_max_base_error_rate
    )
    ch_versions = ch_versions.mix(PUC19_FILTERCONSENSUSREADS.out.versions.first())

    PUC19_RASTAIR_CALL(
        PUC19_FILTERCONSENSUSREADS.out.bam,
        PUC19_FILTERCONSENSUSREADS.out.bai,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(PUC19_RASTAIR_CALL.out.versions)

    PUC19_RASTAIR_CALL_SUMMARY(
        PUC19_RASTAIR_CALL.out.bed
    )
    ch_versions = ch_versions.mix(PUC19_RASTAIR_CALL_SUMMARY.out.versions)

    emit:
    meth_summary  = PUC19_RASTAIR_CALL_SUMMARY.out.tsv        // channel: [ val(meta), path("*.methylation_summary.tsv") ]
    flagstat      = PUC19_SAMTOOLS_FLAGSTAT.out.flagstat       // channel: [ val(meta), path("*.flagstat") ]
    versions      = ch_versions
}
