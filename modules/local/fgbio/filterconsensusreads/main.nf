process FGBIO_FILTERCONSENSUSREADS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::fgbio=2.4.0 bioconda::samtools=1.21"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'
        : 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'}"

    input:
    tuple val(meta), path(consensus_bam)
    tuple val(genome), path(fasta)
    val min_reads
    val min_baseq
    val max_base_error_rate

    output:
    tuple val(meta), path("*.cons.filtered.bam"), emit: bam
    tuple val(meta), path("*.cons.filtered.bam.bai"), emit: bai
    path "versions.yml", emit: versions

    script:
    def fgbio_args        = task.ext.fgbio_args        ?: ''
    def samtools_args     = task.ext.samtools_args     ?: ''
    def samtools_threads  = task.ext.samtools_threads  ?: 4
    def prefix            = task.ext.prefix            ?: "${meta.id}"
    def mem_gb            = task.ext.fgbio_mem_gb      ?: 8

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --async-io=true \\
        --tmp-dir=. \\
        --compression=0 \\
        FilterConsensusReads \\
        --input ${consensus_bam} \\
        --output /dev/stdout \\
        --ref ${fasta} \\
        --min-reads ${min_reads} \\
        --min-base-quality ${min_baseq} \\
        --max-base-error-rate ${max_base_error_rate} \\
        ${fgbio_args} \\
        | samtools sort \\
        --threads ${samtools_threads} \\
        -o ${prefix}.cons.filtered.bam##idx##${prefix}.cons.filtered.bam.bai \\
        --write-index \\
        ${samtools_args};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cons.filtered.bam
    touch ${prefix}.cons.filtered.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
