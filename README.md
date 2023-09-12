# darevskia-pericentromere-analysis

## Introduction

**darevskia-pericentromere-analysis** is a bioinformatics pipeline that can be used to analyse pericentromeric sequences of _Darevskia_ lizards. By default it requires no input and starts from assembled contigs in contigs directory. It performs tandem repeat monomers search and creates FISH probe candidates for subgenomes staining as an output.

However, you can run pipeline from reads, assemble them and run contamination assessment, however it requires a plenty of computional resources. Reads will be downloaded automatically from SRA.

## Usage

You can run the pipeline using:

```bash
nextflow run darevskia-pericentromere-analysis/main.nf \
    -profile <singularity/docker> \
        [--from_fastq] \
            [--enable_magicblast] \
                [--dbdir] \
            [--enable_tarean]
```

<b>Options</b>

    * `--from_fastq` Begin pipeline from reads. Reads in fastq format will be downloaded automatically from SRA. By deafult it provides reads trimming and FastQC analysis.

        * `--enable_magicblast` Run Magic-BLAST contamination assessment. Warning! This step requires a plenty of computional resources. This option can be specified only with `--from_fastq` option.

            * `--dbdir` Allocate Magic-BLAST databases.

        * `--enable_tarean` Enable contigs assembling with TAREAN.Warning! This step requires a plenty of computional resources. This option can be specified only with `--from_fastq` option.

### Pipeline flowchart

![flowchart](https://github.com/nikitin-p/darevskia-pericentromere-analysis/blob/master/pipeline.png)

## Pipeline discussion

Originally, contigs were assembled with TAREAN locally and during re-writing pipeline in Nextflow DSL2, we noticed that there were backward changes in TAREAN commit we used. However, contigs assembled in pipeline are do not differ much from what we had, we decided to provide original pre-assembled contigs for reproducibility.

Also, we developed an R script which is looking for most distant tandem repeat monomers from the whole family of the other subgenome using Levenshtein distance. We suppose it is not the best way to compare sequences of different lengths, so it can be improved in the way of comparing tandem repeat k-mers. However we use the whole monomer approach for reproducibility, as we originally got FISH probes with this method.

## Pipeline output

## Credits

The pipeline was written in Nextflow DSL2 by:

- [Sviatoslav Sidorov](https://github.com/sidorov-si)
- [Pavel Nikitin](https://github.com/nikitin-p)

## Requirements

The pipeline was working fine with Nextflow v21.10.6.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.