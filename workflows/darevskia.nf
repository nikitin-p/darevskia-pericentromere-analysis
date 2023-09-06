include { extract_species } from '../modules/local/custom_functions.nf'
include { extract_reverse_species } from '../modules/local/custom_functions.nf'

include { DOWNLOADREADS } from '../modules/local/downloadreads.nf'
include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
include { TRIMMOMATIC } from '../modules/local/trimmomatic.nf'
// include { MAGICBLAST } from '../modules/local/magicblast.nf'
// include { PARSEMAGICBLAST } from '../modules/local/parsemagicblast.nf'
// include { INTERLACEFASTA } from '../modules/local/interlacefasta.nf'
// include { REPEATEXPLORER } from '../modules/local/repeatexplorer.nf'
// include { PREPROCESSTRF } from '../modules/local/preprocesstrf.nf'
// include { QUAST } from '../modules/local/quast.nf'
// include { TRF } from '../modules/local/trf.nf'
// include { PREPROCESSR } from '../modules/local/preprocessr.nf'
// include { MONOMERPROBE } from '../modules/local/monomerprobe.nf'
// include { RPLOTS } from '../modules/local/rplots.nf'
// include { BOWTIE2_BUILD } from '../modules/nf-core/modules/bowtie2/build/main.nf'
// include { BOWTIE2_CROSS_ALIGN } from '../modules/local/crossalign.nf'
// include { PARSESAM } from '../modules/local/parsesam.nf'
// include { BOWTIE2_CLSAT_ALIGN } from '../modules/local/bowtie2clsatalign.nf'

// Do not forget to replace SRR with the real ones

srr = [
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
]

Channel
    .from( srr )
    .map{ row -> [ row[0], [ row[1], row[2] ] ] }
    .set{ ch_srr }

Channel
    .fromPath('./primer/primer.fasta')
    .set{ primer }

Channel
    .fromPath('./clsat36/clsat36.fasta')
    .set{ clsat36 }

rmd_handler = file( "./rmd/plot_GC_length_distr.Rmd" )

// This should allow to skip everything before if --from_contigs is given
params.from_contigs = false
// This should allow to enable tarean if --enable_tarean is given
params.enable_tarean = false 
// This should allow to enable magicblast if --enable_magicblast is given
params.enable_magicblast = false 

workflow DAREVSKIA {

    // The idea:
    // no flags:
    //     run from contigs to end
    // --from_fastq:
    //     downloads fastq automatically from SRA
    // --enable_tarean:
    //     requires --from_fastq, almost the same as no flags but enables tarean
    // --enable_magicblast:
    //     requires --from_fastq, almost the same as no flags but enables magicblast
    // --dbdir 

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

        // if ( params.enable_magicblast ) {
        //     if (params.db_dir) {
        //         Channel
        //             // path example: '/home/nikitinp/lizards/pipeline/magicblast_db_test/*'
        //             .fromPath(params.db_dir, type: 'dir' )
        //             .set{ db_dir }

        //         MAGICBLAST (
        //             ch_reads,
        //             db_dir
        //         )

        //         PARSEMAGICBLAST (
        //             MAGICBLAST.out.mb_results
        //         )
        //     } else {
        //         exit 1, "Path to MAGIC-Blast databases is not specified!"
        //     }
        // }
        
        // if ( params.enable_tarean ) {
        //     INTERLACEFASTA (
        //         TRIMMOMATIC.out.trimmed_reads_f_p,
        //         TRIMMOMATIC.out.trimmed_reads_r_p
        //     )

        //     REPEATEXPLORER (
        //         INTERLACEFASTA.out.interlaced_reads
        //     )
        // }
        
    }

    // if ( params.enable_tarean ) {
    //     ch_contigs = REPEATEXPLORER.out.repeat_contigs
    // } else {
    //     contigs = [
    //         [
    //             [
    //                 id: "N"
    //             ],
    //             "./contigs/contigs_N.fasta",
    //         ],
    //         [
    //             [
    //                 id: "V"
    //             ],
    //             "./contigs/contigs_V.fasta"
    //         ]
    //     ]

    //     Channel
    //         .from( contigs )
    //         .map{ row -> [ row[0], file(row[1]) ] }
    //         .set{ ch_contigs }
    // }

    // PREPROCESSTRF (
    //     ch_contigs
    // )

    // QUAST (
    //     ch_contigs
    // )

    // TRF (
    //     PREPROCESSTRF.out.top10pc_contigs.concat(PREPROCESSTRF.out.all_contigs)
    // )

    // PREPROCESSR (
    //     TRF.out.trf_dat
    // )

    // MONOMERPROBE (
    //     PREPROCESSR.out.repeats_tsv.filter( ~/.*top10pc.*/ ).collect()
    // )

    // RPLOTS (
    //     rmd_handler,
    //     PREPROCESSTRF.out.all_contigs_tab.filter( ~/.*_N_.*/ ),
    //     PREPROCESSTRF.out.all_contigs_tab.filter( ~/.*_V_.*/ ),
    //     // contigs_n_tab,
    //     // contigs_v_tab,
    //     PREPROCESSR.out.repeats_tsv.filter( ~/N_all.*/ ),
    //     PREPROCESSR.out.repeats_tsv.filter( ~/V_all.*/ )
    //     // units_n_tab,
    //     // units_v_tab
    //     // REPEATEXPLORER.out.repeat_contigs,
    //     // PREPROCESSR.out.all_repeats_tab
    //     //     !do not forget to add length plots to this
    // )

    // BOWTIE2_BUILD (
    //     PREPROCESSTRF.out.all_contigs.map {it -> it[1]}
    // )
        
    // Channel
    //     .fromPath('./darevskia-pericentromere-analysis/probes/probes_*.fasta')
    //     .map { it -> [extract_reverse_species(it), it] }
    //     .set{ ch_selected_probes }

    // ch_contigs_index.join(ch_selected_probes)
    //     .map { it -> [ [id: it[0], single_end: true] , it[1], it[2] ] }
    //     .set{ ch_index_probes }

    // BOWTIE2_CROSS_ALIGN (
    //     ch_index_probes,
    //     false,
    //     false
    // )

    // PARSESAM (
    //     BOWTIE2_ALIGN.out.sam
    // )

    // BOWTIE2_CLSAT_ALIGN (
    //     BOWTIE2_BUILD.out.contigs_index
    //         .map {it -> [extract_species(it), it]},
    //     clsat36
    // )

}