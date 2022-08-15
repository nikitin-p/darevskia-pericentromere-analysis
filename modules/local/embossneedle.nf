process EMBOSSNEEDLE {
    tag '$bam'
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf657eab_5':
        'quay.io/biocontainers/emboss:6.6.0--h440b012_4' }"

    input:
    path bam

    output:
    path "*.bam", emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    
    """
    revseq cl107contig1_V.fasta cl107contig1_revcomp_V.fasta

    needle cl107contig1_revcomp_V.fasta cl1contig21_N.fasta -outfile alignment_for_probes.needle -gapopen 10 -gapextend 0.5

    needle cl118contig1_V.fasta cl1contig21_N.fasta -outfile alignment_for_probes.needle -gapopen 10 -gapextend 0.5
    
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        embossneedle: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
