process TAPS_MULTIQC {
    tag "multiqc_report"
    label 'process_single'

    conda "bioconda::multiqc=1.35"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'biocontainers/multiqc:1.35--pyhdfd78af_1' }"

    input:
    val(trigger)   // ordering trigger — collects upstream channel items to ensure all
                   // publishDir writes are complete before this process runs
    val(genome)    // reference genome name for methylation column label (e.g. CHM13, GRCh38)

    output:
    path "multiqc_report.html",          emit: report
    path "multiqc_report_data",          emit: data
    path "multiqc_report_plots",         emit: plots,    optional: true
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // run_multiQC.py searches params.outdir for all metric files, regenerates *_mqc.tsv
    // via the bin/ scripts, then runs MultiQC. This avoids manually wiring every metric
    // file through Nextflow channels — new metrics are picked up automatically.
    // file() resolves relative outdir against workflow.launchDir so the absolute path is
    // valid from the task work directory.
    def abs_outdir = file(params.outdir)
    """
    python3 ${projectDir}/bin/run_multiQC.py \\
        ${abs_outdir} \\
        --genome "${genome}" \\
        --mqc-outdir ./mqc_work \\
        --multiqc multiqc

    cp    ${abs_outdir}/report/multiqc_report.html .
    cp -r ${abs_outdir}/report/multiqc_report_data/ .
    [ -d  ${abs_outdir}/report/multiqc_report_plots ] \\
        && cp -r ${abs_outdir}/report/multiqc_report_plots/ . \\
        || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.35
        python: 3.11
    END_VERSIONS
    """
}
