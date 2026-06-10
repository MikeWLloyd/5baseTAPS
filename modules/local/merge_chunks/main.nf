process MERGE_CHUNKS {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::samtools=1.21"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'
        : 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.merged.bam"), emit: bam
    tuple val(meta), path("*.merged.bam.bai"), emit: bai, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_mem = task.ext.sort_mem ?: '8G'
    // bams is a list; pass each explicitly to samtools cat
    def bam_list = bams instanceof List ? bams.join(' ') : bams
    """
    samtools cat ${bam_list} \\
        | samtools sort \\
            --template-coordinate \\
            -@ ${task.cpus} \\
            -m ${sort_mem} \\
            -o ${prefix}.merged.bam##idx##${prefix}.merged.bam.bai \\
            -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.bam
    touch ${prefix}.merged.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
