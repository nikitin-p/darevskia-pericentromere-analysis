include { extract_species } from '../modules/local/custom_functions.nf'
include { extract_reverse_species } from '../modules/local/custom_functions.nf'

// decide which to turn on
include { DOWNLOADREADS } from '../modules/local/downloadreads.nf'


include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
// make it optional
// include { MAGICBLAST } from '../modules/local/magicblast.nf'
// include { PARSEMAGICBLAST } from '../modules/local/parsemagicblast.nf'
//

include { TRIMMOMATIC } from '../modules/local/trimmomatic.nf'
// include { INTERLACEFASTA } from '../modules/local/interlacefasta.nf'

// turn off below, but not remove, make it optional
// include { REPEATEXPLORER } from '../modules/local/repeatexplorer.nf'
//

// by default begin from here (--from_fastq)
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

// srr_n_meta = [id: "N", srr: "SRR20851170", single_end: false]
// srr_v_meta = [id: "V", srr: "SRR20851171", single_end: false]

// Use this to specify SRRs

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

// Don't forget to add them in full pipeline!

// ref_prok_rep_genomes
// ref_euk_rep_genomes
// 16S_ribosomal_RNA

// Channel
//     .of( "ref_viroids_rep_genomes", "ref_viruses_rep_genomes" )
//     .set{ ch_magicblast_dbs }

// contigs = [
//     [
//     [
//         id: "N"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/preprocesstrf/contigs_N_top10pc_bt2.fasta",
//     ],
//     [
//     [
//         id: "V"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/preprocesstrf/contigs_V_top10pc_bt2.fasta"
//     ]
// ]

// contigs = [
//     [
//     [
//         id: "N"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/repeatexplorer/output_N",
//     ],
//     [
//     [
//         id: "V"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/repeatexplorer/output_V"
//     ]
// ]

// contigs = [
//     [
//     [
//         id: "N"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/repex_test/output_N",
//     ],
//     [
//     [
//         id: "V"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/repex_test/output_V"
//     ]
// ]

// trf = [
//     [
//         "/home/nikitinp/lizards/pipeline/results/trf/N_contigs_top10pc.fasta.2.5.7.80.10.50.2000.dat"
//     ],
//     [
//         "/home/nikitinp/lizards/pipeline/results/trf/V_contigs_top10pc.fasta.2.5.7.80.10.50.2000.dat"
//     ]
// ]

// trf_meta = [
//     [
//         id: "N"
//     ],
//     [
//         id: "V"
//     ]
// ]

// trf = [
//     [
//     [
//         id: "N"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/trf/contigs_N_top10pc.fasta.2.7.7.80.10.50.500.dat",
//     ],
//     [
//     [
//         id: "V"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/trf/contigs_V_top10pc.fasta.2.7.7.80.10.50.500.dat"
//     ]
// ]

// rtables = [
//     [
//     [
//         id: "repeat_units"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/preprocessr/N_top10pc_tab_bycol.tsv",
//     "/home/nikitinp/lizards/pipeline/results/preprocessr/V_top10pc_tab_bycol.tsv"
//     ]
// ]

// rtables = [
//     [
//     [
//         id: "N"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/preprocessr/N_top10pc_tab_bycol.tsv",
//     ],
//     [
//     [
//         id: "V"
//     ],
//     "/home/nikitinp/lizards/pipeline/results/preprocessr/V_top10pc_tab_bycol.tsv"
//     ]
// ]

// reads = [
//     [
//     [
//         id: "N"
//     ],
//     "/home/nikitinp/lizards/pipeline/reads/N_R1.fastq.gz",
//     "/home/nikitinp/lizards/pipeline/reads/N_R2.fastq.gz"
//     ],
//     [
//     [
//         id: "V"
//     ],
//     "/home/nikitinp/lizards/pipeline/reads/V_R1.fastq.gz",
//     "/home/nikitinp/lizards/pipeline/reads/V_R2.fastq.gz"
//     ]
// ]

// reads = [
//     [
//     [
//         id: "N_sample"
//     ],
//     "/home/nikitinp/lizards/pipeline/subsample/N_sample_R1.fastq.gz",
//     "/home/nikitinp/lizards/pipeline/subsample/N_sample_R2.fastq.gz"
//     ]
// ]

// genome_valentini = [
//     [
//     [
//         id: "V"
//     ], 
//     "/home/nikitinp/lizards/valentini/fasta/ncbi_dataset/data/GCA_024498535.1/GCA_024498535.1_Dval_245_genomic.fna"
//     ]
// ]

// Channel
//     .from( [srr_n_meta, srr_v_meta] )
//     .set{ ch_srr_meta }

// Channel
//     .from( contigs )
//     .map{ row -> [ row[0], file(row[1]) ] }
//     .set{ ch_contigs }

// Channel
//     .from( rtables )
//     .map{ row -> [ row[0], file(row[1]) ] }
//     .set{ ch_rtables }

// Channel
//     .from( trf )
//     .map{ row -> [ row[0], file(row[1]) ] }
//     .set{ ch_trf }

// Channel
//     .from( trf )
//     .map{ row -> file(row[0]) }
//     .set{ ch_trf }

// Channel
//     .from( trf_meta )
//     .set{ ch_trf_meta }

// Channel
//     .from( rtables )
//     .map{ row -> [ row[0], [ file(row[1]), file(row[2]) ] ] }
//     .set{ ch_rtables }

// Channel
//     .fromPath('/home/nikitinp/lizards/pipeline/magicblast_db_test/*', type: 'dir' )
//     .set{ db_dir }

// Channel
//     .from( reads )
//     .map{ row -> [ row[0], [ file(row[1]), file(row[2]) ] ] }
//     .set{ ch_reads }

// Use this to specify SRRs

Channel
    .from( srr )
    .map{ row -> [ row[0], [ row[1], row[2] ] ] }
    .set{ ch_srr }

// Channel
//     .from( genome_valentini )
//     .map{ row -> [ row[0], [ file(row[1]) ] ] }
//     .set{ ch_genome_valentini }

// Reminder: move primer.fasta into repo

// Channel
//     .fromPath('/home/nikitinp/lizards/pipeline/primers/primer.fasta')
//     .set{ primer }

// Channel
//     .fromPath('/home/nikitinp/lizards/pipeline/darevskia-pericentromere-analysis/clsat36/clsat36.fasta')
//     .set{ clsat36 }

// rmd_handler = file( "/home/nikitinp/lizards/pipeline/darevskia-pericentromere-analysis/rmd/plot_GC_length_distr.Rmd" )

// contigs_n_tab = file( "/home/nikitinp/lizards/pipeline/darevskia-pericentromere-analysis/rmd/contigs_N_tab.tsv" )
// contigs_v_tab = file( "/home/nikitinp/lizards/pipeline/darevskia-pericentromere-analysis/rmd/contigs_V_tab.tsv" )
// units_n_tab = file( "/home/nikitinp/lizards/pipeline/darevskia-pericentromere-analysis/rmd/N_all_tab_bycol.tsv" )
// units_v_tab = file( "/home/nikitinp/lizards/pipeline/darevskia-pericentromere-analysis/rmd/V_all_tab_bycol.tsv" )

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
        //     MAGICBLAST (
        //         ch_reads,
        //         db_dir
        //     )

        //     PARSEMAGICBLAST (
        //         MAGICBLAST.out.mb_results
        //     )
        // }
        // else {
            
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
        // else {
            
        // }
        
    }
    else {
        
    }



    // PREPROCESSTRF (
    //     ch_contigs
    //     // REPEATEXPLORER.out.repeat_contigs
    // )

    // QUAST (
    //     // ch_contigs
    //     REPEATEXPLORER.out.repeat_contigs
    // )

    // TRF (
    //     // PREPROCESSTRF.out.fixed_contigs
    //     // PREPROCESSTRF.out.ch_meta,
    //     PREPROCESSTRF.out.top10pc_contigs.concat(PREPROCESSTRF.out.all_contigs)
    // )

    // PREPROCESSR (
    //     // TRF.out.ch_meta,
    //     // TRF.out.trf_dat
    //     // // TRF.out.top10pc_repeats,
    //     // // TRF.out.all_repeats
    //     // ch_trf_meta,
    //     // ch_trf.filter( ~/.*top10pc.*/ )
    //     ch_trf
    // )

    // MONOMERPROBE (
    //     // PREPROCESSR.out.top10pc_repeats_tab
    //     PREPROCESSR.out.repeats_tsv.filter( ~/.*top10pc.*/ ).collect()
    //     // ch_rtables
    // )

    // RPLOTS (
    //     rmd_handler,
    //     contigs_n_tab,
    //     contigs_v_tab,
    //     units_n_tab,
    //     units_v_tab
    //     // REPEATEXPLORER.out.repeat_contigs,
    //     // PREPROCESSR.out.all_repeats_tab
    //     //     !do not forget to add length plots to this
    // )

    // BOWTIE2_BUILD (
    //     // ch_contigs.map {it -> it[1]}
    //     PREPROCESSTRF.out.all_contigs.map {it -> it[1]}
    // )

    // contigs_*_merged_all.fasta
    // PREPROCESSTRF.out.all_contigs
        
    // Channel
    //     .fromPath('./darevskia-pericentromere-analysis/probes/probes_*.fasta')
    //     .map { it -> [extract_reverse_species(it), it] }
    //     .set{ ch_selected_probes }

    // ch_contigs_index.join(ch_selected_probes)
    //     .map { it -> [ [id: it[0], single_end: true] , it[1], it[2] ] }
    //     .set{ ch_index_probes }

    // BOWTIE2_CROSS_ALIGN (
    //     // ch_contigs
    //     // REPEATEXPLORER.out.repeat_contigs
    //     // BOWTIE2BUILD.out.contigs_index
    //     // selected_probes
    //     ch_index_probes,
    //     false,
    //     false
    // )

    // PARSESAM (
    //     BOWTIE2_ALIGN.out.sam
    // )

    // BOWTIE2_CLSAT_ALIGN (
    //     // BOWTIE2_BUILD.out.contigs_index
    //     BOWTIE2_BUILD.out.contigs_index
    //         .map {it -> [extract_species(it), it]},
    //     clsat36
    // )

}