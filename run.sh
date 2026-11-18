#!/bin/bash
#SBATCH --job-name=snakemake_pipeline
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda --cores $SLURM_CPUS_PER_TASK