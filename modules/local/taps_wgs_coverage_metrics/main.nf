process WGS_QC {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'biocontainers/multiqc:1.35--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(mosdepth_summary), path(mosdepth_dist), path(methyl_summary)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wgs_coverage_metrics.py \\
        ${prefix} \\
        ${mosdepth_summary} \\
        ${mosdepth_dist} \\
        ${methyl_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_coverage_metrics_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
