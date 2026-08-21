process TAPS_QC_REPORT {
    label 'process_low'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'biocontainers/multiqc:1.35--pyhdfd78af_1' }"

    input:
    path(mapping_csvs)
    path(methyl_csvs)
    path(wgs_csvs)
    path(duplex_csvs)
    path(vc_csvs)

    output:
    path("taps_qc_report.html"), emit: report
    path("taps_qc_report.tsv"),  emit: tsv
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Collect all CSV files passed as inputs into a single list
    all_csvs=""
    for f in ${mapping_csvs} ${methyl_csvs} ${wgs_csvs} ${duplex_csvs} ${vc_csvs}; do
        [ -f "\$f" ] && all_csvs="\$all_csvs \$f"
    done

    taps_qc_report.py \\
        --out taps_qc_report.html \\
        --tsv taps_qc_report.tsv \\
        \$all_csvs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch taps_qc_report.html
    touch taps_qc_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
