process PARSESAM {
    tag '$bam'
    label 'process_low'

    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'quay.io/biocontainers/YOUR-TOOL-HERE' }"

    input:
    path bam

    output:
    path "*.bam", emit: bam
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsesam: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
