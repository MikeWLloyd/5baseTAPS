process TAPS_MAPPING_METRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'biocontainers/multiqc:1.35--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(main_flagstat,  stageAs: "prededup.flagstat"), \
                     path(lambda_flagstat, stageAs: "lambda_negCtrl.flagstat"), \
                     path(puc19_flagstat,  stageAs: "puc19_posCtrl.flagstat")

    output:
    tuple val(meta), path("*_mqc.tsv"), emit: mqc
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mapping_metrics.py \\
        ${prefix} \\
        ${main_flagstat} \\
        ${lambda_flagstat} \\
        ${puc19_flagstat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_read_partitioning_mqc.tsv
    touch ${prefix}_mapping_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
