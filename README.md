# darevskia-pericentromere-analysis

The pipeline was developed using Nextflow DSL2 by [Pavel Nikitin](https://github.com/nikitin-p) and [Sviatoslav Sidorov](https://github.com/sidorov-si).

## Introduction

Here we describe the pipeline used for the analysis of the pericentromeric sequences of _Darevskia raddei nairensis_ and _D. valentini_, parental species of a hybrid parthenogenetic lizard _D. unisexualis_. In targeted sequencing data obtained from the pericentromeres of the parental species, we search for tandem repeat monomers and, based on them, predict species-specific pericentromeric DNA FISH probes to differentially stain parental subgenomes in the hybrid karyotype.

## Requirements

* `Nextflow v21.10.6` (we developed the pipeline with this Nextflow version) or `Nextflow v23.04.1` (the only other version with which we tested the pipeline).

* `Singularity v3.8.7` or `Docker v23.0.1`. We developed the pipeline using these versions and additionally tested it with `Singularity v3.11.3`.

## Usage

You can run the pipeline using:

```bash
nextflow run darevskia-pericentromere-analysis/main.nf \
    -profile <docker/singularity/.../institute> \
    [--from_fastq [--enable_magicblast --db_dir <dir>] [--enable_tarean]]
```

### Options
Without optional arguments, the pipeline will start from pre-assembled contigs included in this repository.

* `--from_fastq` Start the pipeline from raw reads. The reads in the `fastq` format will be downloaded automatically from SRA. By deafult, the FastQC analysis and read trimming will be performed.

    * `--enable_magicblast` Assess contamination with Magic-BLAST. **Warning!** This step is resource-heavy and may require several days. This option can be specified only with the `--from_fastq` option and requires the `--db_dir` option.

        * `--db_dir <magicblast_dir>` Directory with Magic-BLAST databases (see the **Input** section for details). This option can be specified only with the `--enable_magicblast` option.

    * `--enable_tarean` Assemble contigs with TAREAN. **Warning!** This step is resource-heavy. This option can be specified only with the `--from_fastq` option.

## Input

The pipeline requires no input data. By default, it will start from pre-assembled contigs. Otherwise, if the `--from_fastq` option is specified, the pipeline will start from raw reads that it will download.

We provide contigs that we assembled before the development of this pipeline with TAREAN and used for our analysis. In our pipeline, we also include the contig assembly step using TAREAN. However, it is switched off by default because, when run in this pipeline within a container, TAREAN does not exactly reproduce the contigs that we assembled and used.

For the contamination assessment, we used the following Magic-BLAST databases version 5: `16S_ribosomal_RNA`, `ref_viroids_rep_genomes`, `ref_viruses_rep_genomes`, `ref_prok_rep_genomes`, `ref_euk_rep_genomes`. They can be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/db/v5 using, for example, `lftp` client (see [man lftp](https://linux.die.net/man/1/lftp)). Unpack the downloaded parts of the databases. The directory with the Magic-BLAST databases must have the following structure:

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

The output is placed in the **results** folder. If run without parameters (i.e., starting from the pre-assembed contigs), the pipeline will generate the following subfolders in the **results** folder:

* **quast** - Quality assesment of contigs with QUAST, including `PDF` reports.
* **preprocesstrf** - All contigs and the top 10% highly covered contigs in the `FASTA` format (for Tandem Repeat Finder) and, with stats, in the `TSV` format (for the analysis in R).
* **trf** - Tandem Repeat Finder (TRF) output in its `DAT` format.
* **preprocessr** - TRF output in a tabular format: all tandem repeat monomers and tandem repeat monomers from the top 10% highly covered contigs, annotated with the source contigs.
* **monomerprobe** - Tables of pairwise edit distances between all tandem repeat monomers found in the top 10% highly covered contigs.
* **rplots** - Plots of the GC-content and sequence length distributions, in `PDF`, and the corresponding knitted R notebook, in `HTML`.
* **parsesam** - Tables of predicted DNA FISH probes annotated as "mapped" or "unmapped" to the contigs of the opposite species.
* **bowtie2build** - Bowtie2 index files for the full sets of contigs.
* **bowtie2crossalign** - Alignment of manually selected candidate probes to the contigs of the opposite species, in the `SAM` format.
* **bowtie2clsatalign** - Alignment of the CLsat36 sequence to contigs, in the `SAM` format.
* **extractcontig** - Contigs, in the `FASTA` format, to which the CLsat36 sequence aligned. The predicted DNA FISH probes originate from these contigs.

If `--from_fastq` is set, then, depending on additional options, the pipeline will generate the following additional subfolders:

* **reads** - Gzipped raw reads, in the `FASTQ` format.
* **fastqc** - FASTQC reports on the raw reads, in `HTML` and zipped.
* **trimmomatic** - Trimmed forward/reverse paired/unpaired reads, in the gzipped `FASTQ` format.
* **magicblast** - Contamination assessment report produced by Magic-BLAST as a table, in the `TXT` format.
* **parsemagicblast** - Summarised top predicted contaminants from the Magic-BLAST report, in a space-delimited table, in the `TXT` format.
* **interlacefasta** - Interlaced reads prepared as input for TAREAN, in the `FASTA` format.
* **repeatexplorer** - Contigs assembled with TAREAN, in the `FASTA` format. Importantly, the output of the TAREAN module in the pipeline does not exactly match the pre-assembled contigs.

## Repository structure

* **bin** - R scripts run in Nextflow modules.
* **clsat36/clsat36.fasta** - CLsat36 sequence.
* **conf** - Nextflow configs.
* **contigs** - Pre-assembled contigs.
* **contigs_for_probes** - Contigs with the CLsat monomers from which probes were manually selected.
* **modules** - Nextflow modules implementing the steps of the pipeline.
* **primer/primer.fasta** - DOP-PCR primer used in library preparations.
* **probes** - FASTA files with validated species-specific DNA FISH probes for _D. raddei nairensis_ and _D. valentini_.
* **rmd/plot_GC_length_distr.Rmd** - R Markdown notebook with the analysis of the GC content and length of contigs and tandem repeat monomers.
* **workflows/darevskia.nf** - The workflow that implements the pipeline.
* **CITATIONS.md** - References to the software tools and R packages used in the pipeline, with versions.
* **main.nf** - A wrapper workflow that runs `workflows/darevskia.nf`.
* **nextflow.config** - The main Nextflow config file that includes the config files from the `conf` folder.

## Citations

An extensive list of references for the tools used in the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## To-do

1. **[In progress]** Run the whole pipeline (try all options and different versions of Singularity and Nextflow).
2. Add the real SRR IDs when the paper is published.
