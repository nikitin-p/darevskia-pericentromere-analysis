process PREPROCESSR {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    val meta
    each path(trf_dat)

    output:
    tuple val(meta), path("*.tsv"), emit: repeats_tsv
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #The code below is supposed to work with web-version. We should change it according to CLI version of TRF.
    #awk '{if ($0 ~ /^CL/) {printf $0 "\t"} else {print $0}}' N_6.txt > N_top10pc_tab.tsv
    #awk '{if ($0 ~ /^CL/) {printf $0 "\t"} else {print $0}}' V_6.txt > V_top10pc_tab.tsv
    #cat N_top10pc_tab.tsv | tr '(' '\t' | tr -d ')' | tr '-' '\t' | tr -d ' ' > N_top10pc_tab_bycol.tsv
    #cat V_top10pc_tab.tsv | tr '(' '\t' | tr -d ')' | tr '-' '\t' | tr -d ' ' > V_top10pc_tab_bycol.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PREPROCESSR: 1.0
    END_VERSIONS
    """
}
