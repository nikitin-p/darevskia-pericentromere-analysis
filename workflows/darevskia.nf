// include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
// include { MAGICBLAST } from '../modules/local/magicblast.nf'
// include { PARSEMAGICBLAST } from '../modules/local/parsemagicblast.nf'
// include { TRIMMOMATIC } from '../modules/local/trimmomatic.nf'
// include { INTERLACEFASTA } from '../modules/local/interlacefasta.nf'
// include { REPEATEXPLORER } from '../modules/local/repeatexplorer.nf'
// include { PREPROCESSTRF } from '../modules/local/preprocesstrf.nf'
include { QUAST } from '../modules/local/quast.nf'
// include { TRF } from '../modules/local/trf.nf'

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
//     .fromPath('/home/nikitinp/lizards/pipeline/magicblast_db_test/*', type: 'dir' )
//     .set{ db_dir }

// Channel
//     .from( reads )
//     .map{ row -> [ row[0], [ file(row[1]), file(row[2]) ] ] }
//     .set{ ch_reads }

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
    //     REPEATEXPLORER.out.repeat_contigs
    // )

    QUAST (
        ch_contigs
        // REPEATEXPLORER.out.repeat_contigs
    )

    // TRF (
        // PREPROCESSTRF.out.top10pc_contigs
    // )
}
