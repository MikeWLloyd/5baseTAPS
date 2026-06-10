process RASTAIR_CPG_MASK {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'
        : 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa' }"

    input:
    tuple val(meta), path(bed_gz)

    output:
    tuple val(meta), path("*.rastair_cpg_sites.bed.gz"),     emit: bed
    tuple val(meta), path("*.rastair_cpg_sites.bed.gz.tbi"), emit: tbi
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat ${bed_gz} \\
        | awk '!/^#/ { print \$1"\\t"\$2"\\t"\$3 }' \\
        | awk 'NR==1 { chr=\$1; s=\$2; e=\$3; next }
               \$1==chr && \$2==e { e=\$3; next }
               { print chr"\\t"s"\\t"e; chr=\$1; s=\$2; e=\$3 }
               END { print chr"\\t"s"\\t"e }' \\
        | bgzip -c > ${prefix}.rastair_cpg_sites.bed.gz

    tabix -p bed ${prefix}.rastair_cpg_sites.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version 2>&1 | head -1 | sed 's/bgzip (htslib) //')
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_cpg_sites.bed.gz
    touch ${prefix}.rastair_cpg_sites.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version 2>&1 | head -1 | sed 's/bgzip (htslib) //')
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/tabix (htslib) //')
    END_VERSIONS
    """
}
