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
    grep -v 'Parameters' ${trf_dat} | \\
        sed '/^\$/d' | \\
        tail -n+6 | \\
        awk 'BEGIN { prevseq="" } \\
            { \\
                if ( \$0 ~ /Sequence/ ) { \\
                    prevseq=\$0 \\
                } else { \\
                    print prevseq "\\t" \$0 \\
                } \\
            }' | \\
        awk '{print \$2 "\\t" \$3 "\\t" \$17}' | \\
        tr -d '(' | \\
        tr -d ')' | \\
        tr '-' '\\t' > \\
        ${meta.id}_top10pc_tab_bycol.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PREPROCESSR: 1.0
    END_VERSIONS
    """
}