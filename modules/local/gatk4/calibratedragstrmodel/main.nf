process GATK4_CALIBRATEDRAGSTRMODEL {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0'
        : 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(dict)
    tuple val(meta6), path(str_table)

    output:
    tuple val(meta), path("*.dragstr_model.txt"), emit: dragstr_model
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.toGiga() - 2) : 14
    """
    gatk --java-options "-Xmx${avail_mem}g" \\
        CalibrateDragstrModel \\
        --threads ${task.cpus} \\
        -R ${fasta} \\
        --str-table-path ${str_table} \\
        -I ${bam} \\
        -O ${prefix}.dragstr_model.txt \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dragstr_model.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//' || echo "stub")
    END_VERSIONS
    """
}
