process BCFTOOLS_SPLITVCF {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.17--h3cc50cf_1'
        : 'biocontainers/bcftools:1.17--h3cc50cf_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.filtered-SNP.vcf.gz"),      emit: snp_vcf
    tuple val(meta), path("*.filtered-SNP.vcf.gz.tbi"),  emit: snp_tbi
    tuple val(meta), path("*.filtered-INDEL.vcf.gz"),     emit: indel_vcf
    tuple val(meta), path("*.filtered-INDEL.vcf.gz.tbi"), emit: indel_tbi
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.haplotypecaller"
    """
    bcftools view --type snps   ${vcf} -O z -o ${prefix}.filtered-SNP.vcf.gz
    bcftools index --tbi ${prefix}.filtered-SNP.vcf.gz

    bcftools view --type indels ${vcf} -O z -o ${prefix}.filtered-INDEL.vcf.gz
    bcftools index --tbi ${prefix}.filtered-INDEL.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.haplotypecaller"
    """
    touch ${prefix}.filtered-SNP.vcf.gz   ${prefix}.filtered-SNP.vcf.gz.tbi
    touch ${prefix}.filtered-INDEL.vcf.gz ${prefix}.filtered-INDEL.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //' || echo "stub")
    END_VERSIONS
    """
}
