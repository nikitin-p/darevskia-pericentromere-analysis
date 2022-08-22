process EXTRACTCONTIG {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    path sam
    tuple val(meta), path(contigs)

    output:
    path "*.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:

    """
    grep -A1 $(tail -1 clsat36_mapped_on_${meta.id}.sam | cut -f3) ../preprocesstrf/contigs_N_merged_all.fasta > 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EXTRACTCONTIG: 1.0
    END_VERSIONS
    """
}
