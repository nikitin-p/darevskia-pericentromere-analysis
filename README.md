# darevskia-pericentromere-analysis
A Nextflow pipeline for analysis of pericentromeric sequences of 
Darevskia lizards.

To run pipeline use the following:
```
nextflow run darevskia-pericentromere-analysis/main.nf -profile <singularity/docker>
    [--from_fastq] # Download fastq automatically from SRA
        [--enable_magicblast] # Run Magic-BLAST
            [--dbdir] # Allocate Magic-BLAST databases
        [--enable_tarean] # Enable contigs assembling with TAREAN
```

By default pipeline begins from assembled contigs as TAREAN assembly and Magic-BLAST contamination analysis requires many comuptional resources.

The pipeline was working fine with Nextflow v21.10.6.

### Pipeline flowchart

![flowchart](https://github.com/nikitin-p/darevskia-pericentromere-analysis/blob/master/pipeline.png)
