process EXTRACTCONTIG {
    tag "$species_name"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    tuple val(species_name), path(sam), path(contigs)

    output:
    path "*.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:

    """
    grep -A1 \$(tail -1 ${sam} | cut -f3) ${contigs} > ${species_name}.tmp

    contig_name=\$(grep '^>' ${species_name}.tmp | cut -d" " -f1 | tr -d ">")

    mv ${species_name}.tmp \${contig_name}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EXTRACTCONTIG: 1.0
    END_VERSIONS
    """
}
