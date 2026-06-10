#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TheJacksonLaboratory/5baseTAPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/TheJacksonLaboratory/5baseTAPS
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQUORUM          } from './workflows/fastquorum'
include { RASTAIR_METHYLSEQ   } from './workflows/rastair'
include { GATK_VARIANTCALL as GATK } from './workflows/gatk_variantcall'
include { TAPS_QC_REPORT      } from './modules/local/taps_qc_report/main'
include { TAPS_MULTIQC        } from './modules/local/taps_multiqc/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_fastquorum_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_fastquorum_pipeline'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_fastquorum_pipeline'

if (!params.fasta     && getGenomeAttribute('fasta'))     params.replace('fasta',     getGenomeAttribute('fasta'))
if (!params.fasta_fai && getGenomeAttribute('fasta_fai')) params.replace('fasta_fai', getGenomeAttribute('fasta_fai'))
if (!params.dict      && getGenomeAttribute('dict'))      params.replace('dict',      getGenomeAttribute('dict'))
if (!params.bwamem2   && getGenomeAttribute('bwamem2'))   params.replace('bwamem2',   getGenomeAttribute('bwamem2'))

// If local index paths no longer exist,
// reset to null so SAMTOOLS_FAIDX / SAMTOOLS_DICT / BWAMEM2_INDEX rebuild from the FASTA.
// fasta is intentionally excluded — it resolves to iGenomes S3 and is always accessible.
if (params.fasta_fai && !file(params.fasta_fai).exists()) {
    log.warn "fasta_fai not found at ${params.fasta_fai} — will rebuild from FASTA"
    params.replace('fasta_fai', null)
}
if (params.dict && !file(params.dict).exists()) {
    log.warn "dict not found at ${params.dict} — will rebuild from FASTA"
    params.replace('dict', null)
}
if (params.bwamem2) {
    def bwa_dir = file(params.bwamem2)
    if (!bwa_dir.exists() || bwa_dir.list()?.size() == 0) {
        log.warn "bwamem2 index not found or empty at ${params.bwamem2} — will rebuild from FASTA"
        params.replace('bwamem2', null)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    //
    // WORKFLOW: Run main workflow
    //
    JAXGT_5BASE_TAPS(
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        JAXGT_5BASE_TAPS.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow JAXGT_5BASE_TAPS {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    // Initialize fasta file with meta map:
    fasta = params.fasta ? Channel.fromPath(params.fasta).map { it -> [[id: it.baseName], it] }.collect() : Channel.empty()

    // Set various consensus calling and filtering parameters if not given
    if (params.duplex_seq) {
        if (params.groupreadsbyumi_strategy == '') {
            params.replace("groupreadsbyumi_strategy", 'Paired')
        }
        else if (params.groupreadsbyumi_strategy != 'Paired') {
            log.error("config groupreadsbyumi_strategy must be 'Paired' for duplex-sequencing data")
            exit(1)
        }
        if (params.call_min_reads == '') {
            params.replace("call_min_reads", '1 0 0')
        }
        if (!params.filter_min_reads) {
            params.replace("filter_min_reads", '3 1 1')
        }
    }
    else {
        if (params.groupreadsbyumi_strategy == '') {
            params.replace("groupreadsbyumi_strategy", 'Adjacency')
        }
        else if (params.groupreadsbyumi_strategy == 'Paired') {
            log.error("config groupreadsbyumi_strategy cannot be 'Paired' for non-duplex-sequencing data")
            exit(1)
        }
        if (params.call_min_reads == '') {
            params.replace("call_min_reads", '1')
        }
        if (params.filter_min_reads == '') {
            params.replace("filter_min_reads", '3')
        }
    }

    // WORKFLOW: build indexes if needed
    PREPARE_GENOME(fasta)

    // Gather built indices or get them from the params
    // Built from the fasta file:
    dict = params.dict
        ? Channel.fromPath(params.dict).map { it -> [[id: 'dict'], it] }.collect()
        : PREPARE_GENOME.out.dict
    fasta_fai = params.fasta_fai
        ? Channel.fromPath(params.fasta_fai).map { it -> [[id: 'fai'], it] }.collect()
        : PREPARE_GENOME.out.fasta_fai
    bwamem2 = params.bwamem2
        ? Channel.fromPath(params.bwamem2, type: 'dir').map { it -> [[id: 'bwamem2'], it] }.collect()
        : PREPARE_GENOME.out.bwamem2
    //
    // WORKFLOW: Run pipeline
    //
    FASTQUORUM(
        params,
        samplesheet,
        bwamem2,
        dict,
        fasta,
        fasta_fai,
    )

    //
    // WORKFLOW: Downstream analysis — methylation calling, variant calling, QC
    //
    RASTAIR_METHYLSEQ(
        FASTQUORUM.out.bam,
        FASTQUORUM.out.bai,
        fasta,
        fasta_fai,
        FASTQUORUM.out.unmapped_bam,
        FASTQUORUM.out.prededup_flagstat,
        FASTQUORUM.out.mosdepth_summary,
        FASTQUORUM.out.mosdepth_dist,
    )

    //
    // WORKFLOW: GATK variant calling — non-CpG variants, CpG sites excluded via rastair mask
    //   Skip with --run_gatk false for methylation-only runs.
    //
    ch_vc_metrics_csv  = Channel.empty()
    ch_gatk_mqc_files  = Channel.empty()

    if (params.run_gatk) {
        ch_cpg_mask = RASTAIR_METHYLSEQ.out.cpg_mask
            .join(RASTAIR_METHYLSEQ.out.cpg_mask_tbi)

        GATK(
            FASTQUORUM.out.bam,
            FASTQUORUM.out.bai,
            fasta,
            fasta_fai,
            dict,
            ch_cpg_mask,
            RASTAIR_METHYLSEQ.out.meth_summary,
        )

        ch_vc_metrics_csv  = GATK.out.vc_metrics_csv
        ch_gatk_mqc_files  = GATK.out.mqc_files
    }

    //
    // MODULE: TAPS_MULTIQC — run run_multiQC.py after all workflows complete.
    // Collects lightweight trigger signals to enforce ordering; the script then
    // searches params.outdir for all published metric files automatically.
    // Adding new metrics only requires editing bin/run_multiQC.py — no channel
    // wiring changes needed.
    //
    ch_mqc_trigger = FASTQUORUM.out.multiqc_files
        .mix(RASTAIR_METHYLSEQ.out.multiqc_files)
        .mix(ch_gatk_mqc_files)
        .collect()

    TAPS_MULTIQC(
        ch_mqc_trigger,
        Channel.value(params.genome ?: "Sample"),
    )

    emit:
    multiqc_report = TAPS_MULTIQC.out.report.toList()  // channel: /path/to/multiqc_report.html
}
