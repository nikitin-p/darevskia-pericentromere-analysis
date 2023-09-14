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
    [--from_fastq [--enable_magicblast --dbdir <dir>] [--enable_tarean]]
```

### Options

* `--from_fastq` Start the pipeline from raw reads. The reads in the `fastq` format will be downloaded automatically from SRA. By deafult, the FastQC analysis and read trimming will be performed.

    * `--enable_magicblast` Assess contamination with Magic-BLAST. **Warning!** This step is resource-heavy and may require several days. This option can be specified only with the `--from_fastq` option and requires the `--dbdir` option.

        * `--dbdir <magicblast_dir>` Directory with Magic-BLAST databases (see the **Input** section for details). This option can be specified only with the `--enable_magicblast` option.

    * `--enable_tarean` Assemble contigs with TAREAN. **Warning!** This step is resource-heavy. This option can be specified only with the `--from_fastq` option.

## Input

The pipeline requires no input data. By default, it will start from the pre-assembled contigs. Otherwise, if the `--from_fastq` option is specified, the pipeline will start from raw reads that it will download.

We provide contigs that we assembled before the development of this pipeline with TAREAN and used for our analysis. In our pipeline, we also include the contig assembly step using TAREAN. However, it is switched off by default because it does not exactly reproduce the contigs that we assembled and used.

For the contamination assessment, we used the following Magic-BLAST databases version 5: `ref_euk_rep_genomes`, `16S_ribosomal_RNA`, `ref_prok_rep_genomes`, `ref_viruses_rep_genomes`, `ref_viroids_rep_genomes`. They can be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/db/v5 using, for example, `lftp` client (see [man lftp](https://linux.die.net/man/1/lftp)). The directory with the Magic-BLAST databases must have the following structure:

```
magicblast_dir
├── ref_viroids_rep_genomes
│   └── ref_viroids_rep_genomes.tar.gz 
├── ref_prok_rep_genomes
│   ├── ref_prok_rep_genomes.00.tar.gz
│   ├── ref_prok_rep_genomes.01.tar.gz 
│   ├── ...
...
```

## Output

...

## Repository structure

...

## Pipeline flowchart

![flowchart](https://github.com/nikitin-p/darevskia-pericentromere-analysis/blob/master/pipeline.png)

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## To-do

1. Clarify if the user needs to unpack the Magic-BLAST database *.tar.gz files.
2. Write up the Output section.
3. Write up the Repository structure section.
