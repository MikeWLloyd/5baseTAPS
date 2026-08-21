process SPLIT_FASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::seqkit=2.9.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(fastqs)
    val n_chunks

    output:
    tuple val(meta), path("split_out/R1.part_*.fastq.gz"), path("split_out/R2.part_*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p split_out
    ln -sf ${fastqs[0]} R1.fastq.gz
    ln -sf ${fastqs[1]} R2.fastq.gz

    seqkit split2 \\
        --by-part ${n_chunks} \\
        --read1 R1.fastq.gz \\
        --read2 R2.fastq.gz \\
        --out-dir split_out \\
        --threads ${task.cpus - 2} \\
        --extension .gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    def parts = (1..n_chunks).collect { String.format('%03d', it) }
    """
    mkdir -p split_out
    ${parts.collect { "touch split_out/R1.part_${it}.fastq.gz split_out/R2.part_${it}.fastq.gz" }.join('\n    ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.9.0
    END_VERSIONS
    """
}
