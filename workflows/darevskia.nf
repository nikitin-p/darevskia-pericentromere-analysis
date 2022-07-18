// include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
// include { MAGICBLAST } from '../modules/local/magicblast.nf'
// include { PARSEMAGICBLAST } from '../modules/local/parsemagicblast.nf'
// include { TRIMMOMATIC } from '../modules/local/trimmomatic.nf'
// include { INTERLACEFASTA } from '../modules/local/interlacefasta.nf'
// include { REPEATEXPLORER } from '../modules/local/repeatexplorer.nf'
// include { PREPROCESSTRF } from '../modules/local/preprocesstrf.nf'
// include { QUAST } from '../modules/local/quast.nf'
// include { TRF } from '../modules/local/trf.nf'
// include { PREPROCESSR } from '../modules/local/preprocessr.nf'
// include { MONOMERPROBE } from '../modules/local/monomerprobe.nf'
// include { KMERPROBE } from '../modules/local/kmerprobe.nf'
// include { PYSCRIPTS } from '../modules/local/pyscripts.nf'
include { BOWTIE2BUILD } from '../modules/nf-core/modules/bowtie2/build/main.nf'
// include { BOWTIE2ALIGN } from '../modules/nf-core/modules/bowtie2/align/main.nf'


contigs = [
    [
    [
        id: "N"
    ],
    "/home/nikitinp/lizards/pipeline/results/repeatexplorer/output_N",
    ],
    [
    [
        id: "V"
    ],
    "/home/nikitinp/lizards/pipeline/results/repeatexplorer/output_V"
    ]
]

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

Channel
    .from( contigs )
    .map{ row -> [ row[0], file(row[1]) ] }
    .set{ ch_contigs }

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

Channel
    .fromPath('./darevskia-pericentromere-analysis/probes/*_probes.fasta')
    .set{ selected_probes }

// Reminder: move primer.fasta into repo

// Channel
//     .fromPath('/home/nikitinp/lizards/pipeline/primers/primer.fasta')
//     .set{ primer }

workflow DAREVSKIA {

    // FASTQC ( 
    //     ch_reads 
    // )

    // MAGICBLAST (
    //     ch_reads,
    //     db_dir
    // )

    // PARSEMAGICBLAST (
    //     MAGICBLAST.out.mb_results
    // )

    // TRIMMOMATIC (
    //     ch_reads,
    //     primer
    // )

    // INTERLACEFASTA (
    //     TRIMMOMATIC.out.trimmed_reads_f_p,
    //     TRIMMOMATIC.out.trimmed_reads_r_p
    // )

    // REPEATEXPLORER (
    //     INTERLACEFASTA.out.interlaced_reads
    // )

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

    // KMERPROBE (
    //     // PREPROCESSR.out.top10pc_repeats_tab
    //     PREPROCESSR.out.repeats_tsv.filter( ~/.*top10pc.*/ ).collect()
    //     // ch_rtables
    // )

    // PYSCRIPTS (
    //     REPEATEXPLORER.out.repeat_contigs,
    //     PREPROCESSR.out.all_repeats_tab
            // !do not forget to add length plots to this
    // )

    BOWTIE2BUILD (
        ch_contigs
        // PREPROCESSTRF.out.all_contigs
    )

    // contigs_*_merged_all.fasta
    // PREPROCESSTRF.out.all_contigs

    // BOWTIE2BUILD.out.contigs_index
    // .map {it -> [extract_species(it), it]}

    // def contigs_name = "contigs_N_merged_all.fasta";
    // def m = contigs_name =~ /contigs_([NV])/;
    // print m[0][1]​

    // BOWTIE2ALIGN (
    //     // ch_contigs
    //     // REPEATEXPLORER.out.repeat_contigs
    //     BOWTIE2BUILD.out.contigs_index
    //     selected_probes
    // )

}
