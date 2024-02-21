include { trimSuffix } from './custom_functions'

process REPEATEXPLORER {
    tag "$meta.id"
    label 'process_high'

    container 'sviatsidorov/tarean:1.1'

    input:
    tuple val(meta), path(interlaced_fasta)

    output:
    tuple val(meta), path("contigs_*.fasta"), emit: repeat_contigs
    tuple val(meta), path("output_*"), emit: repeat_output
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

    cp output_${meta.id}/contigs.fasta contigs_${meta.id}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        REPEATEXPLORER: 0.3.7
    END_VERSIONS
    """
}
