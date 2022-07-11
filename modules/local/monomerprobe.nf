process MONOMERPROBE {
    tag "monomer_probes"
    label 'process_high'
    
    container 'rocker/tidyverse:3.6.3'

    input:
    val repeat_units

    output:
    path "*_tr_units_by_distance_geq10.tsv", emit: units
    path "versions.yml", emit: versions

    script:
    """
    arg1=`echo ${repeat_units} | \\
        tr -d ',[]' | \\
        cut -d' ' -f1`
    arg2=`echo ${repeat_units} | \\
        tr -d ',[]' | \\
        cut -d' ' -f2`

    generate_probes_from_monomers.R \$arg1 \$arg2
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MONOMERPROBE: 1.0
    END_VERSIONS
    """
}
