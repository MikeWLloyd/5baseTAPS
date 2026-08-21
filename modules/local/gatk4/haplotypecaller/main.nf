process GATK4_HAPLOTYPECALLER {
    tag "${meta.id}:${interval}"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0'
        : 'oras://ghcr.io/tjongh/dragmap-gatk:4.6.1.0' }"

    input:
    tuple val(meta),  path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fasta_fai)
    tuple val(meta5), path(dict)
    tuple val(meta6), path(cpg_mask), path(cpg_mask_tbi)
    tuple val(meta7), path(dragstr_model)
    val   interval

    output:
    tuple val(meta), path("*.haplotypecaller.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.haplotypecaller.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}.${interval}"
    def avail_mem   = task.memory ? (task.memory.toGiga() - 1) : 14
    def hmm_threads = Math.min(4, Math.max(1, task.cpus - 2))
    """
    gatk --java-options "-Xmx${avail_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads=2" \\
        HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${prefix}.haplotypecaller.vcf.gz \\
        -L ${interval} \\
        --exclude-intervals ${cpg_mask} \\
        --dragen-mode true \\
        --dragstr-params-path ${dragstr_model} \\
        --native-pair-hmm-threads ${hmm_threads} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${interval}"
    """
    touch ${prefix}.haplotypecaller.vcf.gz
    touch ${prefix}.haplotypecaller.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version | grep "Genome Analysis Toolkit" | sed 's/.*v//' || echo "stub")
    END_VERSIONS
    """
}
