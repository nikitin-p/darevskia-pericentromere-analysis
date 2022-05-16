process RSCRIPTS {
    tag "$meta.id"
    label 'process_high'
    
    container 'rocker/tidyverse:3.6.3'

    input:
    tuple val(meta), path(repeat_units)

    output:
    tuple val(meta), path("*_tr_probes_by_distance_geq10.tsv"), emit: probes
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_probes.R
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RSCIPTS: 1.0
    END_VERSIONS
    """
}
