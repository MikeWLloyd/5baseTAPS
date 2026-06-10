process RASTAIR_VCF {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/tjongh/rastair:2.0.0-jaxGT' :
        'community.wave.seqera.io/library/rastair:2.0.0--5717eb9e4a70193d' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), val(parsed_trim_OT)
    tuple val(meta6), val(parsed_trim_OB)

    output:
    tuple val(meta), path("*.rastair_call.txt"),    emit: txt
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nt_OT_to_trim = meta.trim_OT ?: parsed_trim_OT
    def nt_OB_to_trim = meta.trim_OB ?: parsed_trim_OB

    """
    rastair call \\
        --nOT ${nt_OT_to_trim} \\
        --nOB ${nt_OB_to_trim} \\
        --fasta-file ${fasta} \\
        --vcf ${prefix}.rastair.vcf.gz \\
        --vcf-info-fields DP,AF,MQ,NAB,NOI,M5mC_Strands,CPG,CPGnovo,AS_SB \\
        --vcf-format-fields GT,DP,M5mC \\
        --cpg-novo-min-mapq 30 \\
        --m-vaf-min 0.1 \\
        --m-min-depth 5 \\
        --min-baseq 20 \\
        --threads ${task.cpus} \\
        ${bam} 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_call.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version 2>&1 || echo "stub")
    END_VERSIONS
    """
}
