
process QUAST {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::quast=5.0.2" : null)
    container 'quay.io/biocontainers/quast:5.0.2--py37pl5321h09c1ff4_7'

    input:
    tuple val(meta), path(contigs) // path(contigs_dir)

    output:
    tuple val(meta), path("contigs_*"), emit: contigs_qc
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // ${contigs_dir}/contigs.fasta
    """
    quast ${contigs} \\
        -o contigs_${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
