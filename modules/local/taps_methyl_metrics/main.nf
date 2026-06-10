process TAPS_METHYL_METRICS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'biocontainers/multiqc:1.35--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(main_tsv), path(lambda_tsv), path(puc19_tsv)

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome = task.ext.genome ?: "Sample"
    """
    methyl_metrics.py \\
        ${prefix} \\
        ${main_tsv} \\
        ${lambda_tsv} \\
        ${puc19_tsv} \\
        "${genome}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_methyl_controls_mqc.tsv
    touch ${prefix}_methyl_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
