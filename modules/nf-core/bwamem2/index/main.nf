process BWAMEM2_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa' :
        'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa-mem2"), emit: index
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa-mem2
    bwa-mem2 \\
        index \\
        $args \\
        -p bwa-mem2/${prefix} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa-mem2
    touch bwa-mem2/${prefix}.0123
    touch bwa-mem2/${prefix}.amb
    touch bwa-mem2/${prefix}.ann
    touch bwa-mem2/${prefix}.bwt.2bit.64
    touch bwa-mem2/${prefix}.pac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
    END_VERSIONS
    """
}
