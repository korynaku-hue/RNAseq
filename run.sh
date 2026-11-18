#!/bin/bash
#SBATCH --job-name=snakemake_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake --use-conda --cores 8