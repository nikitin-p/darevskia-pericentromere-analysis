process BOWTIE2_CLSAT_ALIGN {
    label "process_low"

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6" : null)
    container "${ workflow.containerEngine == "singularity" && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0" :
        "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0" }"

    input:
    tuple val(species_name), path(index)
    each path(clsat36)
    
    output:
    path "*.sam"        , emit: sam
    path  "versions.yml", emit: versions

    script:
    def reads_args = "-U ${clsat36}"

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/.rev.1.bt2//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/.rev.1.bt2l//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 \\
        -f \\
        -x \$INDEX \\
        $reads_args \\
        --threads $task.cpus \\
        -S clsat36_mapped_on_${species_name}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
