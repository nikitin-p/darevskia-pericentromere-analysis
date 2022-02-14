include { FASTQC } from '../modules/nf-core/modules/fastqc/main.nf'

reads = [
    [
    [
        id: "N_sample"
    ],
    "/home/nikitinp/lizards/pipeline/subsample/N_sample_R1.fastq.gz",
    "/home/nikitinp/lizards/pipeline/subsample/N_sample_R2.fastq.gz"
    ]
]

Channel
    .from( reads )
    .map{ row -> [ row[0], [ file(row[1]), file(row[2]) ] ] }
    .set{ ch_reads }

workflow DAREVSKIA {
    FASTQC( ch_reads )
}
