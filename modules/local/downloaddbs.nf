process DOWNLOADDBS {
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::magicblast=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magicblast:1.6.0--h95f258a_0':
        'quay.io/biocontainers/magicblast:1.6.0--h95f258a_0' }"

    output:
    path "ref_*", emit: db
    path "versions.yml", emit: versions

    script:
    """
    #mkdir ref_viroids_rep_genomes
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/ref_viroids_rep_genomes*" -P ref_viroids_rep_genomes/
    cd ref_viroids_rep_genomes/
    md5sum -c ref_viroids_rep_genomes.tar.gz.md5
    tar -xzvf *.tar.gz
    cd ../

    #mkdir ref_viruses_rep_genomes
    #wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_viruses_rep_genomes*" -P ref_viruses_rep_genomes/
    #md5sum -c ref_viruses_rep_genomes/ref_viruses_rep_genomes.tar.gz.md5

    #mkdir ref_16S_ribosomal_RNA
    #wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "16S_ribosomal_RNA*"
    #md5sum -c 16S_ribosomal_RNA.tar.gz.md5

    #mkdir ref_prok_rep_genomes
    #wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_prok_rep_genomes*"
    #md5sum -c ref_prok_rep_genomes.*.tar.gz.md5

    #mkdir ref_euk_rep_genomes
    #wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/*" -A "ref_euk_rep_genomes*"
    #md5sum -c ref_euk_rep_genomes.*.tar.gz.md5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget 2>&1 | head -1 | cut -d" " -f1,2)
    END_VERSIONS
    """
}