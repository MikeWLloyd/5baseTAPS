//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWAMEM2_INDEX } from '../../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_DICT } from '../../../modules/nf-core/samtools/dict/main'

workflow PREPARE_GENOME {
    take:
    fasta // channel: [mandatory] fasta

    main:
    versions = Channel.empty()

    if (!params.bwamem2)   BWAMEM2_INDEX(fasta)
    if (!params.fasta_fai) SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []])
    if (!params.dict)      SAMTOOLS_DICT(fasta)

    // Gather versions of all tools used
    if (!params.bwamem2)   versions = versions.mix(BWAMEM2_INDEX.out.versions)
    if (!params.fasta_fai) versions = versions.mix(SAMTOOLS_FAIDX.out.versions)

    emit:
    bwamem2   = params.bwamem2   ? Channel.empty() : BWAMEM2_INDEX.out.index.collect()
    dict      = params.dict      ? Channel.empty() : SAMTOOLS_DICT.out.dict.collect()
    fasta_fai = params.fasta_fai ? Channel.empty() : SAMTOOLS_FAIDX.out.fai.collect()
    versions // channel: [ versions.yml ]
}
