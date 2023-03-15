process DOWNLOADREADS {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    tuple val(meta), val(srrs)

    output:
    tuple val(meta), tuple path("*_1.fastq.gz"), path("*_2.fastq.gz"), emit: fastq
    path "versions.yml"           , emit: versions

    script:
    """
    wget ${srrs[0]} -O ${meta.id}_1.fastq.gz
    wget ${srrs[1]} -O ${meta.id}_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget -V | head -1| cut -d" " -f3)
    END_VERSIONS
    """
}
