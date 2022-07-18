include { trimSuffix } from './custom_functions'

process REPEATEXPLORER {
    tag "$meta.id"
    label 'process_high'

    container 'sviatsidorov/tarean:1.1'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'shub://repeatexplorer/repex_tarean:0.3.8.dbaa07f':
    //     'kavonrtep/repeatexplorer:2.3.8' }"

    input:
    tuple val(meta), path(interlaced_fasta)

    output:
    tuple val(meta), path("output_*"), emit: repeat_contigs
    path("log_*.txt"), emit: tarean_log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export PATH=/opt/conda/bin:/home/toor/repex_tarean:\$PATH
    export PYTHONHASHSEED=0
    /home/toor/repex_tarean/seqclust \\
        -p ${interlaced_fasta} \\
        -l log_${meta.id}.txt \\
        -c $task.cpus \\
        -v output_${meta.id} \\
        -tax METAZOA3.0 \\
        -opt ILLUMINA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        REPEATEXPLORER: 0.3.7
    END_VERSIONS
    """
}