include { extract_species } from '../modules/local/custom_functions.nf'
include { extract_reverse_species } from '../modules/local/custom_functions.nf'

include { DOWNLOADREADS } from '../modules/local/downloadreads.nf'
include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
include { TRIMMOMATIC } from '../modules/local/trimmomatic.nf'
include { MAGICBLAST } from '../modules/local/magicblast.nf'
include { PARSEMAGICBLAST } from '../modules/local/parsemagicblast.nf'
include { INTERLACEFASTA } from '../modules/local/interlacefasta.nf'
include { REPEATEXPLORER } from '../modules/local/repeatexplorer.nf'
include { PREPROCESSTRF } from '../modules/local/preprocesstrf.nf'
include { QUAST } from '../modules/local/quast.nf'
include { TRF } from '../modules/local/trf.nf'
include { PREPROCESSR } from '../modules/local/preprocessr.nf'
include { MONOMERPROBE } from '../modules/local/monomerprobe.nf'
include { RPLOTS } from '../modules/local/rplots.nf'
include { BOWTIE2_BUILD } from '../modules/nf-core/modules/bowtie2/build/main.nf'
include { BOWTIE2_CROSS_ALIGN } from '../modules/local/crossalign.nf'
include { PARSESAM } from '../modules/local/parsesam.nf'
include { BOWTIE2_CLSAT_ALIGN } from '../modules/local/bowtie2clsatalign.nf'
include { EXTRACTCONTIG } from '../modules/local/extractcontig.nf'

// Replace test SRR numbers below with the real ones and uncomment the definition of the srr structure.
// Real SRR numbers are SRR25825523 for D. raddei nairensis and SRR25825522 for D. valentini.

/*srr = [
    [
    [
        id: "N"
    ],
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/070/SRR20851170/SRR20851170_1.fastq.gz",
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/070/SRR20851170/SRR20851170_2.fastq.gz"
    ],
    [
    [
        id: "V"
    ],
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/071/SRR20851171/SRR20851171_1.fastq.gz",
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR208/071/SRR20851171/SRR20851171_2.fastq.gz"
    ]
]*/

Channel
    .from( srr )
    .map{ row -> [ row[0], [ row[1], row[2] ] ] }
    .set{ ch_srr }

Channel
    .fromPath('./darevskia-pericentromere-analysis/primer/primer.fasta')
    .set{ primer }

Channel
    .fromPath('./darevskia-pericentromere-analysis/clsat36/clsat36.fasta')
    .set{ clsat36 }

rmd_handler = file( "./darevskia-pericentromere-analysis/rmd/plot_GC_length_distr.Rmd" )

params.from_fastq = false
params.enable_tarean = false 
params.enable_magicblast = false
params.db_dir = false

workflow DAREVSKIA {

    if ( params.enable_tarean && !params.from_fastq) {
        exit 1, "ERROR: TAREAN can be enabled only when --from_fastq is set."
    }

    if ( params.enable_magicblast && !params.from_fastq ) {
        exit 1, "ERROR: Magic-BLAST can be enabled only when --from_fastq is set."
    }

    if ( params.enable_magicblast && !params.db_dir ) {
        exit 1, "ERROR: Path to Magic-BLAST databases is not specified."
    }

    if ( params.db_dir && !params.enable_magicblast ) {
        exit 1, "ERROR: --enable_magicblast is not set."
    }

    if ( params.db_dir && !params.from_fastq) {
        exit 1, "ERROR: --from_fastq is not set."
    }

    if ( params.db_dir ) {
        if (params.db_dir.toString() == "true") {
            exit 1, "ERROR: Specify a path to a folder with Magic-BLAST databases after --db_dir. See README."
        }
        file(params.db_dir, checkIfExists: true)
    }

    if ( params.from_fastq ) {
        DOWNLOADREADS(
            ch_srr
        )

        DOWNLOADREADS.out.fastq
            .map { it -> [it[0], [it[1], it[2]]]}
            .set{ ch_reads }

        FASTQC ( 
            ch_reads 
        )

        TRIMMOMATIC (
            ch_reads,
            primer
        )

        if ( params.enable_magicblast ) {
            if (params.db_dir) {
                Channel
                    // This fromPath(...) factory may not
                    // list all directories within the params.db_dir directory,
                    // but just the params.db_dir directory itself,
                    // if run with Nextflow 23.04.1 (or, possibly, higher).
                    // Needs addtional testing.
                    .fromPath(params.db_dir, type: 'dir' )
                    .set{ db_dir }

                db_dir.view()

                MAGICBLAST (
                    ch_reads,
                    db_dir
                )
        
                PARSEMAGICBLAST (
                    MAGICBLAST.out.mb_results
                )
            } else {
               exit 1, "ERROR: Specify a path to a folder with Magic-BLAST databases after --db_dir. See README."
            }
        }
        
        if ( params.enable_tarean ) {
            INTERLACEFASTA (
                TRIMMOMATIC.out.trimmed_reads_f_p,
                TRIMMOMATIC.out.trimmed_reads_r_p
            )

            REPEATEXPLORER (
                INTERLACEFASTA.out.interlaced_reads
            )
        }
        
    }

    if ( params.enable_tarean ) {
        ch_contigs = REPEATEXPLORER.out.repeat_contigs
    } else {
        contigs = [
            [
                [
                    id: "N"
                ],
                "./darevskia-pericentromere-analysis/contigs/contigs_N.fasta",
            ],
            [
                [
                    id: "V"
                ],
                "./darevskia-pericentromere-analysis/contigs/contigs_V.fasta"
            ]
        ]

        Channel
            .from( contigs )
            .map{ row -> [ row[0], file(row[1]) ] }
            .set{ ch_contigs }
    }

    PREPROCESSTRF (
        ch_contigs
    )

    QUAST (
        ch_contigs
    )

    TRF (
        PREPROCESSTRF.out.top10pc_contigs.concat(PREPROCESSTRF.out.all_contigs)
    )

    PREPROCESSR (
        TRF.out.trf_dat
    )

    MONOMERPROBE (
        PREPROCESSR.out.repeats_tsv.filter( ~/.*top10pc.*/ ).collect()
    )

    RPLOTS (
        rmd_handler,
        PREPROCESSTRF.out.all_contigs_tab.filter( ~/.*_N_.*/ ),
        PREPROCESSTRF.out.all_contigs_tab.filter( ~/.*_V_.*/ ),
        PREPROCESSR.out.repeats_tsv.filter( ~/.*N_all_.*/ ),
        PREPROCESSR.out.repeats_tsv.filter( ~/.*V_all_.*/ )
    )

    BOWTIE2_BUILD (
        PREPROCESSTRF.out.all_contigs.map {it -> it[1]}
    )

    BOWTIE2_BUILD.out.contigs_index
        .map {it -> [extract_species(it), it]}
        .set{ ch_contigs_index }
        
    Channel
        .fromPath('./darevskia-pericentromere-analysis/probes/probes_*.fasta')
        .map { it -> [extract_reverse_species(it), it] }
        .set{ ch_selected_probes }

    ch_contigs_index.join(ch_selected_probes)
        .map { it -> [ [id: it[0], single_end: true] , it[1], it[2] ] }
        .set{ ch_index_probes }

    BOWTIE2_CROSS_ALIGN (
        ch_index_probes,
        false,
        false
    )

    PARSESAM (
        BOWTIE2_CROSS_ALIGN.out.sam
    )

    BOWTIE2_CLSAT_ALIGN (
        ch_contigs_index,
        clsat36
    )

    EXTRACTCONTIG (
        BOWTIE2_CLSAT_ALIGN.out.sam
            .map {it -> [extract_species(it), it]}
            .join(PREPROCESSTRF.out.all_contigs
                .map {it -> [extract_species(it[1]), it[1]]}
            )
    )

}
