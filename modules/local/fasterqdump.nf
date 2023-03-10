process SRATOOLS_FASTERQDUMP {
    tag "$meta.id"
    label 'process_medium'

    // conda "bioconda::sra-tools=2.11.0 conda-forge::pigz=2.6"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6a9ff0e76ec016c3d0d27e0c0d362339f2d787e6-0' :
    //     'quay.io/biocontainers/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6a9ff0e76ec016c3d0d27e0c0d362339f2d787e6-0' }"
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6a9ff0e76ec016c3d0d27e0c0d362339f2d787e6-0'

    input:
    val meta

    output:
    tuple val(meta), path(fastq), emit: reads
    path "versions.yml"         , emit: versions

    script:

    // WARNING: Paired-end data extracted by fasterq-dump (--split-3 the default)
    // always creates *_1.fastq *_2.fastq files but sometimes also
    // an additional *.fastq file for unpaired reads which we ignore here.
    fastq = meta.single_end ? '*.fastq.gz' : '*_R{1,2}.fastq.gz'
    def outfile = meta.single_end ? "${meta.id}.fastq" : "${meta.id}"

    """
    fasterq-dump \\
        --threads $task.cpus \\
        --outfile $outfile \\
        ${meta.srr}

    mv ${outfile}_1.fastq ${outfile}_R1.fastq

    mv ${outfile}_2.fastq ${outfile}_R2.fastq

    pigz \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
