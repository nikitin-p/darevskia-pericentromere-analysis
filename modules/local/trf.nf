process TRF {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::trf=4.09.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trf:4.09.1--hec16e2b_2':
        'quay.io/biocontainers/trf:4.09.1--hec16e2b_2' }"

    input:
    tuple val(meta), path(contigs_fasta)
    // each path(contigs_fasta)

    output:
    tuple val(meta), path("*.dat"), emit: trf_dat
    // val meta, emit: ch_meta
    // path "*.dat",        emit: trf_dat
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    trf ${contigs_fasta} 2 5 7 80 10 50 2000 -d || echo "0"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        TRF: trf -v | head -2 | tail -1 | awk '{print \$5}'
    END_VERSIONS
    """
}