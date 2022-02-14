include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'
// include { READ_NAMES } from '../modules/local/read_names/main.nf'

reads = [
    [
        id: "N_sample"
    ],
    "/home/nikitinp/lizards/pipeline/subsample/N_sample_R1.fastq.gz",
    "/home/nikitinp/lizards/pipeline/subsample/N_sample_R2.fastq.gz"
]

// reads = [
//     [
//         [
//             id: "testx"
//         ],
//         "https://github.com/hartwigmedical/testdata/raw/master/100k_reads_hiseq/TESTX/TESTX_H7YRLADXX_S1_L001_R1_001.fastq.gz"
//     ],
//     [
//         [
//             id: "testy"
//         ],
//         "https://github.com/hartwigmedical/testdata/raw/master/100k_reads_hiseq/TESTY/TESTY_H7YRLADXX_S1_L001_R1_001.fastq.gz"
//     ]
// ]

// Channel
//     .from( reads )
//     .map{ row -> [ row[0], [ file(row[1]) ] ] }
//     .view()
    // .map{ row -> [ row[0], [ file(row[1]), file(row[2]) ] ] }
    // .set{ ch_reads }

// Channel
//     .from( reads )
//     .map{ row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set{ ch_reads }

workflow DAREVSKIA {
//     FASTQC( ch_reads )
//     // READ_NAMES( ch_reads )
Channel
    .from( reads )
    .map{ row -> [ row[0], [ file(row[1]) ] ] }
    .view()
}
