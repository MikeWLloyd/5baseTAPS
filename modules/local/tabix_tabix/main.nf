process TABIX_TABIX {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bed_gz)

    output:
    tuple val(meta), path("*.tbi"), emit: tbi
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tabix -p bed ${bed_gz}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //')
    END_VERSIONS
    """

    stub:
    """
    touch ${bed_gz}.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: stub
    END_VERSIONS
    """
}
