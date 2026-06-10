process MOSDEPTH {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::mosdepth=0.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/biocontainers/mosdepth:0.3.4--hd299d5a_0'
        : 'quay.io/biocontainers/mosdepth:0.3.4--hd299d5a_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),     emit: summary
    tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --quantize 0:1:5:10:20: \\
        ${args} \\
        ${prefix} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mosdepth.summary.txt
    touch ${prefix}.mosdepth.global.dist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //' || echo "stub")
    END_VERSIONS
    """
}
