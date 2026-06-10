process GATK4_MERGEVCFS {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0'
        : 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0' }"

    input:
    tuple val(meta),  path(vcfs)
    tuple val(meta2), path(dict)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}.haplotypecaller"
    def avail_mem = task.memory ? (task.memory.toGiga() - 1) : 6
    def input_args = vcfs.collect { "--INPUT ${it}" }.join(' \\\n        ')
    """
    gatk --java-options "-Xmx${avail_mem}g" \\
        MergeVcfs \\
        ${input_args} \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --OUTPUT ${prefix}.vcf.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.haplotypecaller"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//' || echo "stub")
    END_VERSIONS
    """
}
