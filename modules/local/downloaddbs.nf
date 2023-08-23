process DOWNLOADDBS {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    input:
    val(db)

    output:
    tuple val(db), emit: db_gz
    tuple, emit: json
    tuple, emit: md5
    path "versions.yml"           , emit: versions

    script:
    """
    wget ${srrs[0]} -O ${meta.id}_1.fastq.gz
    wget ${srrs[1]} -O ${meta.id}_2.fastq.gz

    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_viroids_rep_genomes*"
    md5sum -c ref_viroids_rep_genomes.tar.gz.md5

    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_viruses_rep_genomes*"
    md5sum -c ref_viruses_rep_genomes.tar.gz.md5

    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "16S_ribosomal_RNA*"
    md5sum -c 16S_ribosomal_RNA.tar.gz.md5

    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_prok_rep_genomes*"
    md5sum -c ref_prok_rep_genomes.*.tar.gz.md5

    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_euk_rep_genomes*"
    md5sum -c ref_euk_rep_genomes.*.tar.gz.md5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget 2>&1 | head -1 | cut -d" " -f1,2)
    END_VERSIONS
    """
}