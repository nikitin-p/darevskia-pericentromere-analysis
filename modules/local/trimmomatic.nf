process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::trimmomatic=0.39" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    
    tuple val(meta), path(reads)
    each path(primer)

    output:
    tuple val(meta), path("*_trimmed_f_p.fastq.gz"), emit: trimmed_reads_f_p
    tuple val(meta), path("*_trimmed_r_p.fastq.gz"), emit: trimmed_reads_r_p
    tuple val(meta), path("*_trimmed_*_u.fastq.gz"), emit: trimmed_reads_u
    path "*_trimlog.txt" , emit: trimlog
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_trimmed"
    """
    trimmomatic \\
        PE \\
        -phred33 \\
        -threads $task.cpus \\
        -trimlog ${prefix}_trimlog.txt \\
        ${reads[0]} \\
        ${reads[1]} \\
        ${prefix}_f_p.fastq.gz \\
        ${prefix}_f_u.fastq.gz \\
        ${prefix}_r_p.fastq.gz \\
        ${prefix}_r_u.fastq.gz \\
        HEADCROP:25 \\
        ILLUMINACLIP:${primer}:8:30:10 \\
        ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:20

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version 2>&1 | tail -1)
    END_VERSIONS
    """
}
