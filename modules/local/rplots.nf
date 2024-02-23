process RPLOTS {
    label 'process_low'
    stageInMode 'copy'

    container 'sviatsidorov/r_machine:1.1'

    input:
    path rmd
    path contigs_n
    path contigs_v
    path units_n
    path units_v

    output:
    path "contig_gc_content.pdf"    , emit: contigsgc
    path "units_gc_content.pdf"     , emit: unitsgc
    path "contigs_length.pdf"       , emit: contigslength
    path "units_length.pdf"         , emit: unitslength
    path "plot_GC_length_distr.html", emit: knitted_html
    path "versions.yml"             , emit: versions

    script:
    """
    render_rmd.R $rmd ${contigs_n} ${contigs_v} ${units_n} ${units_v}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
        R_tidyr: \$(Rscript -e 'packageVersion("tidyr")' | awk '{print \$2}' | tr -d "‘’")
        R_stringr: \$(Rscript -e 'packageVersion("stringr")' | awk '{print \$2}' | tr -d "‘’")
        
    END_VERSIONS
    """
}
