process GATK4_VARIANTFILTRATION {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0'
        : 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0' }"

    input:
    tuple val(meta),  path(vcf)
    tuple val(meta2), path(tbi)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)
    tuple val(meta5), path(dict)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? (task.memory.toGiga() - 1) : 6
    """
    gatk --java-options "-Xmx${avail_mem}g" \\
        VariantFiltration \\
        -R ${fasta} \\
        -V ${vcf} \\
        --filter-expression "QUAL < 10.4139" \\
        --filter-name "DRAGENHardQUAL" \\
        -O ${prefix}.filtered.vcf.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//' || echo "stub")
    END_VERSIONS
    """
}
