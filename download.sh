#!/bin/bash

#SBATCH --job-name=downloading raw data
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00

export PATH=$PATH:/scratch/lustre/home/koku1049/2ND/tools/sratoolkit.3.4.1-ubuntu64/bin

ACCESSION_FILE="accessions.txt"

mkdir -p raw_data

while IFS= read -r id || [[ -n "$id" ]]; do
        [[ -z "$id" ]] && continue
        echo "Downloading sample: $id"
        fasterq-dump --split-files --outdir raw_data "$id"

done < "$ACCESSION_FILE"