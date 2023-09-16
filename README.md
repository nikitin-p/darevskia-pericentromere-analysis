# darevskia-pericentromere-analysis

The pipeline was developed using Nextflow DSL2 by [Pavel Nikitin](https://github.com/nikitin-p) and [Sviatoslav Sidorov](https://github.com/sidorov-si).

## Introduction

Here we describe the pipeline used for the analysis of the pericentromeric sequences of _Darevskia raddei nairensis_ and _D. valentini_, parental species of a hybrid parthenogenetic lizard _D. unisexualis_. In targeted sequencing data obtained from the pericentromeres of the parental species, we search for tandem repeat monomers and predict species-specific pericentromeric DNA FISH probes to differentially stain the subgenomes in the hybrid karyotype.

## Requirements

* `Nextflow v21.10.6` (we developed the pipeline with this version of Nextflow; we did not test it with later versions).

* `Singularity` or `Docker`. We developed the pipeline using `Singularity v3.8.7` and `Docker v23.0.1` but did not test it with other versions of these containerisation systems.

## Usage

You can run the pipeline using:

```bash
nextflow run darevskia-pericentromere-analysis/main.nf \
    -profile <docker/singularity/.../institute> \
    [--from_fastq [--enable_magicblast --db_dir <dir>] [--enable_tarean]]
```

### Options

* `--from_fastq` Start the pipeline from raw reads. The reads in the `fastq` format will be downloaded automatically from SRA. By deafult, the FastQC analysis and read trimming will be performed.

    * `--enable_magicblast` Assess contamination with Magic-BLAST. **Warning!** This step is resource-heavy and may require several days. This option can be specified only with the `--from_fastq` option and requires the `--db_dir` option.

        * `--db_dir <magicblast_dir>` Directory with Magic-BLAST databases (see the **Input** section for details). This option can be specified only with the `--enable_magicblast` option.

    * `--enable_tarean` Assemble contigs with TAREAN. **Warning!** This step is resource-heavy. This option can be specified only with the `--from_fastq` option.

## Input

The pipeline requires no input data. By default, it will start from the pre-assembled contigs. Otherwise, if the `--from_fastq` option is specified, the pipeline will start from raw reads that it will download.

We provide contigs that we assembled before the development of this pipeline with TAREAN and used for our analysis. In our pipeline, we also include the contig assembly step using TAREAN. However, it is switched off by default because it does not exactly reproduce the contigs that we assembled and used.

For the contamination assessment, we used the following Magic-BLAST databases version 5: `16S_ribosomal_RNA`, `ref_viroids_rep_genomes`, `ref_viruses_rep_genomes`, `ref_prok_rep_genomes`, `ref_euk_rep_genomes`. They can be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/db/v5 using, for example, `lftp` client (see [man lftp](https://linux.die.net/man/1/lftp)). The directory with the Magic-BLAST databases must have the following structure:

```
magicblast_dir
├── ref_viroids_rep_genomes
│   ├── ref_viroids_rep_genomes.ndb
│   ├── ...
│   ├── ref_viroids_rep_genomes.nto
│   ├── taxdb.btd
│   └── taxdb.bti
├── ref_prok_rep_genomes
│   ├── ref_prok_rep_genomes.01.nhr
│   ├── ...
│   ├── ref_prok_rep_genomes.01.nsq
│   ├── ref_prok_rep_genomes.02.nhr
│   ├── ...
│   ├── ref_prok_rep_genomes.02.nsq
│   ├── ...
│   ├── taxdb.btd
│   └── taxdb.bti
...
```

## Pipeline flowchart

![flowchart](https://github.com/nikitin-p/darevskia-pericentromere-analysis/blob/master/pipeline.png)

## Output

```
.
├── bowtie2clsatalign
│   ├── clsat36_mapped_on_N.sam
│   ├── clsat36_mapped_on_V.sam
│   └── versions.yml
├── embossneedle
│   ├── alignment_for_probes.needle
│   ├── CL107Contig1_V_revcomp.fasta
│   └── versions.yml
├── extractcontig
│   ├── CL107Contig1_V.fasta
│   ├── CL1Contig21_N.fasta
│   └── versions.yml
├── fastqc
│   ├── N_1_fastqc.html
│   ├── N_1_fastqc.zip
│   ├── N_2_fastqc.html
│   ├── N_2_fastqc.zip
│   ├── V_1_fastqc.html
│   ├── V_1_fastqc.zip
│   ├── V_2_fastqc.html
│   ├── V_2_fastqc.zip
│   └── versions.yml
├── interlacefasta
│   ├── N_trimmed.fasta
│   ├── V_trimmed.fasta
│   └── versions.yml
├── magicblast
│   ├── N_16S_ribosomal_RNA_output.txt
│   ├── N_ref_prok_rep_genomes_output.txt
│   ├── N_ref_viroids_rep_genomes_output.txt
│   ├── N_ref_viruses_rep_genomes_output.txt
│   ├── V_16S_ribosomal_RNA_output.txt
│   ├── V_ref_viruses_rep_genomes_output.txt
│   ├── V_ref_prok_rep_genomes_output.txt
│   ├── V_ref_viroids_rep_genomes_output.txt
│   └── versions.yml
├── monomerprobe
│   ├── N_tr_units_by_distance_geq10.tsv
│   ├── N_tr_units_by_distance.tsv
│   ├── versions.yml
│   ├── V_tr_units_by_distance_geq10.tsv
│   └── V_tr_units_by_distance.tsv
├── parsemagicblast
│   ├── N_16S_ribosomal_RNA_histogram.txt
│   ├── N_ref_prok_rep_genomes_histogram.txt
│   ├── N_ref_viroids_rep_genomes_histogram.txt
│   ├── N_ref_viruses_rep_genomes_histogram.txt
│   ├── V_16S_ribosomal_RNA_histogram.txt
│   ├── versions.yml
│   ├── V_ref_prok_rep_genomes_histogram.txt
│   ├── V_ref_viroids_rep_genomes_histogram.txt
│   └── V_ref_viruses_rep_genomes_histogram.txt
├── parsesam
│   ├── probe_mapping_on_N.tsv
│   ├── probe_mapping_on_V.tsv
│   └── versions.yml
├── preprocessr
│   ├── N_all_tab_bycol.tsv
│   ├── N_top10pc_tab_bycol.tsv
│   ├── V_all_tab_bycol.tsv
│   ├── versions.yml
│   └── V_top10pc_tab_bycol.tsv
├── preprocesstrf
│   ├── contigs_N_all_seq.fasta
│   ├── contigs_N_merged_all.fasta
│   ├── contigs_N_tab.tsv
│   ├── contigs_N_top10pc_bt2.fasta
│   ├── contigs_N_top10pc.fasta
│   ├── contigs_N_top10pc_seq.fasta
│   ├── contigs_V_all_seq.fasta
│   ├── contigs_V_merged_all.fasta
│   ├── contigs_V_tab.tsv
│   ├── contigs_V_top10pc_bt2.fasta
│   ├── contigs_V_top10pc.fasta
│   ├── contigs_V_top10pc_seq.fasta
│   └── versions.yml
├── quast
│   ├── contigs_N
│   │   ├── basic_stats
│   │   │   ├── contigs_GC_content_plot.pdf
│   │   │   ├── cumulative_plot.pdf
│   │   │   ├── GC_content_plot.pdf
│   │   │   └── Nx_plot.pdf
│   │   ├── icarus.html
│   │   ├── icarus_viewers
│   │   │   └── contig_size_viewer.html
│   │   ├── quast.log
│   │   ├── report.html
│   │   ├── report.pdf
│   │   ├── report.tex
│   │   ├── report.tsv
│   │   ├── report.txt
│   │   ├── transposed_report.tex
│   │   ├── transposed_report.tsv
│   │   └── transposed_report.txt
│   ├── contigs_V
│   │   ├── basic_stats
│   │   │   ├── contigs_GC_content_plot.pdf
│   │   │   ├── cumulative_plot.pdf
│   │   │   ├── GC_content_plot.pdf
│   │   │   └── Nx_plot.pdf
│   │   ├── icarus.html
│   │   ├── icarus_viewers
│   │   │   └── contig_size_viewer.html
│   │   ├── quast.log
│   │   ├── report.html
│   │   ├── report.pdf
│   │   ├── report.tex
│   │   ├── report.tsv
│   │   ├── report.txt
│   │   ├── transposed_report.tex
│   │   ├── transposed_report.tsv
│   │   └── transposed_report.txt
│   └── versions.yml
├── reads
│   ├── N_1.fastq.gz
│   ├── N_2.fastq.gz
│   ├── V_1.fastq.gz
│   ├── V_2.fastq.gz
│   └── versions.yml
├── repeatexplorer
│   ├── output_N
│   │   └── contigs.fasta
│   └── output_V
│       └── contigs.fasta
├── rplots
│   ├── contig_gc_content.pdf
│   ├── contigs_length.pdf
│   ├── plot_GC_length_distr.html
│   ├── units_gc_content.pdf
│   ├── units_length.pdf
│   └── versions.yml
├── rscipts
│   ├── N_tr_probes_by_distance_geq10.tsv
│   ├── V_tr_probes_by_distance_geq10.tsv
│   └── versions.yml
├── trf
│   ├── contigs_N_all_seq.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_N_merged_all.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_N_merged_all.fasta.2.7.7.80.10.50.500.dat
│   ├── contigs_N_top10pc.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_N_top10pc.fasta.2.7.7.80.10.50.500.dat
│   ├── contigs_N_top10pc_seq.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_V_all_seq.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_V_merged_all.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_V_merged_all.fasta.2.7.7.80.10.50.500.dat
│   ├── contigs_V_top10pc.fasta.2.5.7.80.10.50.2000.dat
│   ├── contigs_V_top10pc.fasta.2.7.7.80.10.50.500.dat
│   ├── contigs_V_top10pc_seq.fasta.2.5.7.80.10.50.2000.dat
│   ├── N_contigs_top10pc.fasta.2.5.7.80.10.50.2000.dat
│   ├── V_contigs_top10pc.fasta.2.5.7.80.10.50.2000.dat
│   └── versions.yml
└── trimmomatic
    ├── N_trimmed_f_p.fastq
    ├── N_trimmed_f_p.fastq.gz
    ├── N_trimmed_f_u.fastq
    ├── N_trimmed_f_u.fastq.gz
    ├── N_trimmed_r_p.fastq
    ├── N_trimmed_r_p.fastq.gz
    ├── N_trimmed_r_u.fastq
    ├── N_trimmed_r_u.fastq.gz
    ├── N_trimmed_trimlog.txt
    ├── V_trimmed_trimlog.txt
    ├── V_trimmed_f_p.fastq
    ├── V_trimmed_f_p.fastq.gz
    ├── V_trimmed_f_u.fastq
    ├── V_trimmed_f_u.fastq.gz
    ├── V_trimmed_r_p.fastq
    ├── V_trimmed_r_p.fastq.gz
    ├── V_trimmed_r_u.fastq
    ├── V_trimmed_r_u.fastq.gz
    └── versions.yml
```

## Repository structure

```
.
├── bin
│   ├── generate_probes_from_kmers.R
│   ├── generate_probes_from_monomers.R
│   └── render_rmd.R
├── CITATIONS.md
├── clsat36
│   └── clsat36.fasta
├── conf
│   ├── base.config
│   ├── igenomes.config
│   └── modules.config
├── contigs
│   ├── contigs_N.fasta
│   └── contigs_V.fasta
├── contigs_for_probes
│   ├── cl107contig1_V.fasta
│   └── cl1contig21_N.fasta
├── LICENSE
├── main.nf
├── modules
│   ├── local
│   │   ├── bowtie2clsatalign.nf
│   │   ├── crossalign.nf
│   │   ├── custom_functions.nf
│   │   ├── downloaddbs.nf
│   │   ├── downloadreads.nf
│   │   ├── embossneedle.nf
│   │   ├── extractcontig.nf
│   │   ├── fasterqdump.nf
│   │   ├── interlacefasta.nf
│   │   ├── kmerprobe.nf
│   │   ├── magicblast.nf
│   │   ├── monomerprobe.nf
│   │   ├── parsemagicblast.nf
│   │   ├── parsesam.nf
│   │   ├── preprocessr.nf
│   │   ├── preprocesstrf.nf
│   │   ├── quast.nf
│   │   ├── repeatexplorer.nf
│   │   ├── rplots.nf
│   │   ├── trf.nf
│   │   └── trimmomatic.nf
│   └── nf-core
│       └── modules
│           ├── bowtie2
│           │   └── build
│           │       ├── main.nf
│           │       └── meta.yml
│           ├── bwa
│           │   ├── index
│           │   │   ├── main.nf
│           │   │   └── meta.yml
│           │   └── mem
│           │       ├── main.nf
│           │       └── meta.yml
│           └── fastqc
│               ├── main.nf
│               └── meta.yml
├── modules.json
├── nextflow.config
├── nf-core
│   └── modules
│       ├── bowtie2
│       │   └── align
│       │       ├── main.nf
│       │       └── meta.yml
│       └── nf-core
│           └── bwa
│               ├── index
│               │   ├── main.nf
│               │   └── meta.yml
│               └── mem
│                   ├── main.nf
│                   └── meta.yml
├── pipeline.png
├── primer
│   └── primer.fasta
├── probes
│   ├── probes_N.fasta
│   └── probes_V.fasta
├── README.md
├── rmd
│   ├── contigs_N_tab.tsv
│   ├── contigs_V_tab.tsv
│   ├── N_all_tab_bycol.tsv
│   ├── plot_GC_length_distr.Rmd
│   └── V_all_tab_bycol.tsv
├── temp.nf
└── workflows
    └── darevskia.nf
```

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## To-do

1. Clarify if the user needs to unpack the Magic-BLAST database *.tar.gz files.
2. Write up the Output section.
3. Write up the Repository structure section.
