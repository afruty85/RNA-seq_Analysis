#!/bin/bash
#SBATCH --job-name=ferret_download
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# ============================================================
# RNA-seq Final Project 2025
# Download raw FASTQ data for Kim et al. (2022) ferret SARS-CoV-2 dataset
# BioProject: PRJNA783897
# ============================================================

BASE_DIR="/home/FCAM/afrutos/ISG/Final_project2025"
ACC_LIST="${BASE_DIR}/SRR_Acc_ferret_List.txt"

mkdir -p ${BASE_DIR}/data/raw
cd ${BASE_DIR}/data/raw || exit 1

module load sratoolkit

if [[ ! -f "$ACC_LIST" ]]; then
    echo "ERROR: Accession list not found at $ACC_LIST"
    exit 1
fi

echo "Starting download from NCBI SRA..."
date

# Loop over accessions (compatible with SRA Toolkit 3.x)
for SRR in $(cat "$ACC_LIST"); do
    echo "Downloading $SRR..."
    fasterq-dump --split-files --threads 10 --outdir ${BASE_DIR}/data/raw $SRR
done

echo "Download complete."
date

echo "Compressing FASTQ files..."
gzip ${BASE_DIR}/data/raw/*.fastq

echo "All downloads and compression complete."
echo "Files saved in: ${BASE_DIR}/data/raw"
