process ALIGN_BAM {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'
        : 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'}"

    input:
    tuple val(meta), path(unmapped_bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(bwa_dir)
    val sort_type

    output:
    tuple val(meta), path("*.mapped.bam"), emit: bam
    tuple val(meta), path("*.mapped.bam.bai"), emit: bai, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_fastq_args = task.ext.samtools_fastq_args ?: ''
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def bwa_args = task.ext.bwa_args ?: ''
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_threads      = task.ext.samtools_threads      ?: 2
    def samtools_sort_threads = task.ext.samtools_sort_threads ?: samtools_threads
    def samtools_sort_mem     = task.ext.samtools_sort_mem     ?: ''
    def bwa_cpu_reserve       = task.ext.bwa_cpu_reserve       ?: 6
    def filter_unmapped       = task.ext.filter_unmapped       ?: false
    def bwa_threads           = Math.max(1, task.cpus - bwa_cpu_reserve)
    def fgbio_threads         = Math.max(1, bwa_cpu_reserve - samtools_threads - samtools_sort_threads)
    def fgbio_mem_gb = 12
    def extra_command = ""

    if (!task.memory) {
        log.info('[fgbio ZipperBams] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.')
    }
    else if (fgbio_mem_gb > task.memory.giga) {
        fgbio_mem_gb = Math.max(1, (int)(task.memory.giga - 1))
    }

    if (sort_type == "none") {
        fgbio_zipper_bams_output = prefix + ".mapped.bam"
        fgbio_zipper_bams_compression = 1
    }
    else {
        fgbio_zipper_bams_output = "/dev/stdout"
        fgbio_zipper_bams_compression = 0
        def view_filter = filter_unmapped ? "| samtools view -@ ${samtools_threads} -bF 12 - " : ""
        extra_command = " ${view_filter}| samtools sort "
        extra_command += samtools_sort_args
        if (sort_type == "template-coordinate") {
            extra_command += " --template-coordinate"
        }
        else {
            if (sort_type != "coordinate") {
                log.info('[samtools sort] Unknown sort - defaulting to coordinate.')
            }
            extra_command += " --write-index"
        }
        extra_command += " -@ ${samtools_sort_threads}"
        if (samtools_sort_mem) extra_command += " -m ${samtools_sort_mem}"
        extra_command += " -o " + prefix + ".mapped.bam##idx##" + prefix + ".mapped.bam.bai"
        extra_command += " -"
    }

    """
    # The real path to the bwa-mem2 index prefix
    BWA_INDEX_PREFIX=`find -L ./ -name "*.bwt.2bit.64" | sed 's/.bwt.2bit.64//'`

    samtools fastq -@ ${samtools_threads} ${samtools_fastq_args} ${unmapped_bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${bwa_threads} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g -XX:ActiveProcessorCount=${fgbio_threads} \\
            --compression ${fgbio_zipper_bams_compression} \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${unmapped_bam} \\
            --ref ${fasta} \\
            --output ${fgbio_zipper_bams_output} \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            ${fgbio_args} \\
            ${extra_command};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_command = sort_type != "template-coordinate" ? "touch ${prefix}.mapped.bam.bai" : ""
    """
    touch ${prefix}.mapped.bam
    ${index_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
