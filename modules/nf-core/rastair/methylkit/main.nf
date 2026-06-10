process RASTAIR_METHYLKIT {
    label 'process_single'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.rastair_methylkit.txt.gz"), emit: methylkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Converts rastair call output to methylKit format.
    // Column layout: chr start end name beta_est strand unmod($7) mod($8) ... coverage ...
    // Follows rastair_call_to_methylkit.sh logic from 0.8.2:
    //   cov = unmod + mod; freqC = mod*100/cov; freqT = unmod*100/cov
    //   chrBase uses $3 (end = 1-based position)
    """
    zcat ${bed} \\
        | awk 'BEGIN {
                   print "chrBase\\tchr\\tbase\\tstrand\\tcoverage\\tfreqC\\tfreqT"
               }
               /^#/ { next }
               {
                   cov = \$7 + \$8
                   strand = (\$6 == "+") ? "F" : "R"
                   pct1 = (cov == 0) ? "0" : sprintf("%.2f", \$8 * 100 / cov)
                   pct2 = (cov == 0) ? "0" : sprintf("%.2f", \$7 * 100 / cov)
                   print \$1 ":" \$3 "\\t" \$1 "\\t" \$3 "\\t" strand "\\t" cov "\\t" pct1 "\\t" pct2
               }' \\
        | gzip -c > ${prefix}.rastair_methylkit.txt.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_methylkit.txt.gz
    """
}
