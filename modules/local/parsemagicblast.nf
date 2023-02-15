include { trimSuffix } from './custom_functions'

process PARSEMAGICBLAST {
    tag "${trimSuffix(magicblast_output.simpleName, '_output')}"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    path(magicblast_output)

    output:
    path("*_histogram.txt"), emit: mb_histogram
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def input_name  = "${trimSuffix(magicblast_output.simpleName, '_output')}"

    """
    READCOUNT=`<${magicblast_output} \\
        tail -n +4 | \\
        wc -l`
    
    <${magicblast_output} \\
        tail -n +4 | \\
        awk '{print \$2}' | \\
        sort | \\
        uniq -c | \\
        sort -k1,1nr | \\
        head -5 |
        awk -F" " -v var="\${READCOUNT}" '{print (\$1 / var * 100) "% " \$2}' > \\
        ${input_name}_histogram.txt \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PARSEMAGICBLAST: 1.0
    END_VERSIONS
    """
}
