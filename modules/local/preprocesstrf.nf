process PREPROCESSTRF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    tuple val(meta), path(contigs_dir)

    output:
    tuple val(meta), path("contigs_*_top10pc_seq.fasta"), emit: top10pc_contigs
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk '{if (\$0 ~ /^>/) {print "\n" \$0} else {printf \$0}}' ${contigs_dir}/contigs.fasta > contigs_${meta.id}_merged.fasta

    grep '^>' ${contigs_dir}/contigs.fasta > contigs_${meta.id}_headers.txt

    cat contigs_${meta.id}_headers.txt | \\
        tr '(' '|' | \\
        tr -d ')' | \\
        tr '-' '|' | \\
        tr -d ' ' | \\
        tr '.' '|' | \\
        sort -t '|' -k3,3nr > \\
        contigs_${meta.id}_headers_sorted.txt

    x=$(bc <<< "`wc -l contigs_${meta.id}_headers_sorted.txt | cut -d " " -f1` / 10 + 1")

    head -\${x} contigs_${meta.id}_headers_sorted.txt > contigs_${meta.id}_headers_sorted_top10pc.txt 

    for contig in `awk -F"|" '{print \$1}' contigs_${meta.id}_headers_sorted_top10pc.txt`; do \\
        grep -w -A 1 "^\$contig" contigs_${meta.id}_merged.fasta >> contigs_${meta.id}_top10pc_seq.fasta; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PREPROCESSTRF: 1.0
    END_VERSIONS
    """
}