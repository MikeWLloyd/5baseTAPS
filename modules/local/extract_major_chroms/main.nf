process EXTRACT_MAJOR_CHROMS {
    tag "${meta.id}"
    label 'process_single'

    // No container: only POSIX awk + sort + uniq are needed (host shell).

    input:
    tuple val(meta), path(meth_summary_tsv)

    output:
    tuple val(meta), path("*.major_chroms.txt"), emit: chroms
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract column 2 (chromosome names) from rastair methylation_summary.tsv,
    # keep only major chromosomes: numeric (chr1..chr22) and chrX/chrY.
    # Filters out 'ALL', chrM, chrEBV, and any chrName containing an underscore
    # (which excludes *_random, *_alt, chrUn_*, GL/KI alt contigs, etc.).
    awk -F'\\t' 'NR>1 && \$2 ~ /^chr([0-9]+|X|Y)\$/ { print \$2 }' \\
        ${meth_summary_tsv} \\
        | sort -u \\
        > ${prefix}.major_chroms.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/.*Awk //;s/,.*//' || echo "posix")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "chr1\\nchr2\\nchrX\\nchrY\\n" > ${prefix}.major_chroms.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: stub
    END_VERSIONS
    """
}
