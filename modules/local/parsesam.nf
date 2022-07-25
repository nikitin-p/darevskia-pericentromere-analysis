process PARSESAM {
    tag '$sam'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    tuple val(meta), path(sam)

    output:
    path "*.tsv", emit: probe_mapping
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    <${sam} \\
        grep -v '@' | \\
        awk '{ \\
            if (\$2 == 4) { \\
                print \$1 "\tunmapped" \\
            } else { \\
                print \$1 "\tmapped" \\
                } \\
            }' > probe_mapping_on_${meta.id}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PARSESAM: 1.0
    END_VERSIONS
    """
}
