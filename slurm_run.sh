#!/bin/bash
#SBATCH --job-name=nextflow-example
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

module purge
module load Nextflow
module load Singularity

srun nextflow run main.nf \
    -profile crick
