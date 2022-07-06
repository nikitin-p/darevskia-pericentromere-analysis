process RSCRIPTS {
    tag "probes"
    label 'process_high'
    
    container 'rocker/tidyverse:3.6.3'

    input:
    val repeat_units

    output:
    path "*_tr_probes_by_distance_geq10.tsv", emit: probes
    path "versions.yml", emit: versions

    script:
    """
    arg1=`echo ${repeat_units} | \\
        tr -d ',[]' | \\
        cut -d' ' -f1`
    arg2=`echo ${repeat_units} | \\
        tr -d ',[]' | \\
        cut -d' ' -f2`

    generate_probes.R \$arg1 \$arg2
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RSCRIPTS: 1.0
    END_VERSIONS
    """
}
