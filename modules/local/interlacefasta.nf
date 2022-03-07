include { trimSuffix } from './custom_functions'

process INTERLACEFASTA {
    // tag "$meta.id"
    tag "${trimSuffix(forward_reads.simpleName, '_f_p')}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    tuple val(meta), path(forward_reads)
    tuple val(meta), path(reverse_reads)

    output:
    tuple val(meta), path("*.fasta"), emit: interlaced_reads
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = "${trimSuffix(forward_reads.simpleName, '_f_p')}"
    def prefix_forward   = "${forward_reads.simpleName}"
    def prefix_reverse   = "${reverse_reads.simpleName}"

    """
    gzip -d ${forward_reads}
    gzip -d ${reverse_reads}

    awk 'NR%4==1||NR%4==2' ${forward_reads.baseName} | \\
        tr '@' '>' > ${prefix_forward}.fa

    awk 'NR%4==1||NR%4==2' ${reverse_reads.baseName} | \\
        tr '@' '>' > ${prefix_reverse}.fa

    <${prefix_forward}.fa \\
        awk '{if (\$0 ~ /^>/) {printf (NR + 1) / 2 " f|"} else {print}}' > ${prefix_forward}.fa.txt

    <${prefix_reverse}.fa \\
        awk '{if (\$0 ~ /^>/) {printf (NR + 1) / 2 " r|"} else {print}}' > ${prefix_reverse}.fa.txt

    cat ${prefix_forward}.fa.txt ${prefix_reverse}.fa.txt | \\
        sort -k1,1n | \\
        awk -F' ' '{print ">" \$1 \$2}' | \\
        tr '|' '\\n' > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        INTERLACEFASTA: 1.0
    END_VERSIONS
    """
}
