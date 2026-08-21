process GATK4_COMPOSESTRTABLEFILE {
    tag "genome"
    label 'process_low'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0'
        : 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(dict)

    output:
    tuple val(meta), path("*.str.table"), emit: str_table
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.toGiga() - 2) : 14
    """
    gatk --java-options "-Xmx${avail_mem}g" \\
        ComposeSTRTableFile \\
        -R ${fasta} \\
        -O ${prefix}.str.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.str.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//' || echo "stub")
    END_VERSIONS
    """
}
