include { trimSuffix } from './custom_functions'

process MAGICBLAST {
    tag "$meta.id"
    label 'process_long'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    path("*_output.txt"), emit: mb_results
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = "${trimSuffix(reads[0].simpleName, '_1')}_${db.simpleName}"
    //def prefix   = "${trimSuffix(reads[0].simpleName, '_R1')}_${db.simpleName}"
    
    """
    magicblast \\
        $args \\
        -num_threads $task.cpus \\
        -infmt fastq \\
        -outfmt tabular \\
        -query ${reads[0]} \\
        -query_mate ${reads[1]} \\
        -db ${db}/${db.simpleName} \\
        -out ${prefix}_output.txt \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magicblast: \$(magicblast -version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
