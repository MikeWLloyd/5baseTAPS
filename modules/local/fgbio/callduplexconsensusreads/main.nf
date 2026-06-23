process FGBIO_CALL_DUPLEX {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::fgbio=2.4.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'
        : 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'}"

    input:
    tuple val(meta), path(grouped_bam)
    val min_reads
    val min_baseq

    output:
    tuple val(meta), path("*.cons.unmapped.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info('[fgbio CallDuplexConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.')
    }
    else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CallDuplexConsensusReads \\
        --input ${grouped_bam} \\
        --output ${prefix}.cons.unmapped.bam \\
        --min-reads ${min_reads} \\
        --min-input-base-quality ${min_baseq} \\
        --threads ${task.cpus} \\
        ${args};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cons.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
