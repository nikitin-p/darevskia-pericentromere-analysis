process EMBOSSNEEDLE {
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf657eab_5':
        'quay.io/biocontainers/emboss:6.6.0--h440b012_4' }"

    input:
    tuple path(cl1contig21), path(cl107contig1)

    output:
    path "*revcomp.fasta", emit: revcomp_fasta
    path "*.needle", emit: needle
    path "versions.yml", emit: versions

    script:
    
    """
    revseq ${cl107contig1} CL107Contig1_revcomp.fasta

    needle ${cl1contig21} CL107Contig1_revcomp.fasta -outfile alignment_for_probes.needle -gapopen 10 -gapextend 0.5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EMBOSS: \$(needle --version 2>&1 | cut -d":" -f2)
    END_VERSIONS
    """
}