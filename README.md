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

## Output

After execution of the whole pipeline (i.e. with all parameters), all the output files will be located in `results` folder, which consists of the folowwing subfolders:
* `reads` - compressed (gz) reads in `FASTQ` format
* `fastqc` - FASTQC raw reads quality assessment in `HTML` and `ZIP` format
* `trimmomatic` - trimmed forward/reverse paired/unpaired reads in `FASTQ` format
* `magicblast` - Magic-BLAST raw output in `TXT` format
* `parsemagicblast` - processed Magic-BLAST output in `TXT` format
* `interlacefasta` - interlaced reads for TAREAN assembly in `FASTA` format
* `repeatexplorer` - contigs assembled with TAREAN in two separate folders in `FASTA` format
* `quast` - assembly quality assesment with QUAST in two separate folders in `PDF`, `TEX`, `TSV` and `TXT` formats
* `preprocesstrf` - all contigs and top 10% highly-covered contigs in `FASTA` format and in `TSV` format with metadata
* `trf` - Tandem Repeat Finder (TRF) output in `DAT` format
* `preprocessr` - all tandem repeat monomers and tandem repeat monomers from top 10% highly-covered contigs in `TSV` format
* `monomerprobe` - suggested FISH probes from top 10% highly-covered contigs in `TSV` format
* `rplots` - GC-content and sequences length distributions in `PDF` format and knitted R notebook in `HTML` format
* `bowtie2clsatalign` - CLsat36 alignment on contigs in `SAM` format
* `parsesam` - processed contigs on which CLsat36 aligned in `TSV` format
* `extractcontig` - contigs from which selected FISH probes originate from in `FASTA` format

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
