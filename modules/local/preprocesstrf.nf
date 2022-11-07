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
    // tuple val(meta), path("*_seq.fasta"), emit: fixed_contigs
    tuple val(meta), path("contigs_*_top10pc.fasta"), emit: top10pc_contigs
    tuple val(meta), path("contigs_*_merged_all.fasta"), emit: all_contigs
    tuple val(meta), path("contigs_*_tab.tsv"), emit: all_contigs_tab
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk '{if (\$0 ~ /^>/) {if (NR != 1) {printf "\\n"} print \$0} else {printf \$0}}' ${contigs_dir}/contigs.fasta > contigs_${meta.id}_merged_all.fasta

    cat contigs_${meta.id}_merged_all.fasta | \\
        tr '(' '|' | \\
        tr -d ')' | \\
        tr '-' '|' | \\
        tr -d ' ' | \\
        awk 'BEGIN{RS=">"}{print \$1"\t"\$2;}' | \\
        tr '|' '\t' > \\
        contigs_${meta.id}_tab.tsv

    grep '^>' ${contigs_dir}/contigs.fasta > contigs_${meta.id}_headers.txt

    cat contigs_${meta.id}_headers.txt | \\
        tr '(' '|' | \\
        tr -d ')' | \\
        tr '-' '|' | \\
        tr -d ' ' | \\
        tr '.' '|' | \\
        sort -t '|' -k3,3nr > \\
        contigs_${meta.id}_headers_sorted.txt

    w=\$(wc -l contigs_${meta.id}_headers_sorted.txt | cut -d " " -f1)
    x=\$(awk -v var=\$w '{printf "%.0f", var / 10}' <(echo "1"))

    head -\${x} contigs_${meta.id}_headers_sorted.txt > contigs_${meta.id}_headers_sorted_top10pc.txt 

    for contig in `awk -F"|" '{print \$1}' contigs_${meta.id}_headers_sorted_top10pc.txt`; do \\
        grep -w -A 1 "^\$contig" contigs_${meta.id}_merged_all.fasta >> contigs_${meta.id}_top10pc.fasta; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PREPROCESSTRF: 1.0
    END_VERSIONS
    """
}