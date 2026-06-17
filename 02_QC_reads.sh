#!/bin/bash
#SBATCH --job-name=ferret_qc
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# ============================================================
# RNA-seq Final Project 2025
# Step 2: Quality Control (FastQC + MultiQC)
# ============================================================

BASE_DIR="/home/FCAM/afrutos/ISG/Final_project2025"
RAW_DIR="${BASE_DIR}/data/raw"
QC_DIR="${BASE_DIR}/results/qc"

mkdir -p ${QC_DIR}

module load fastqc
module load multiqc

echo "Starting FastQC analysis..."
date

# Run FastQC on all FASTQ.gz files
fastqc ${RAW_DIR}/*.fastq.gz --outdir ${QC_DIR} --threads 8

echo "FastQC complete."
date

echo "Running MultiQC summary..."
multiqc ${QC_DIR} --outdir ${QC_DIR}

echo "QC pipeline complete."
date
echo "FastQC + MultiQC reports saved in: ${QC_DIR}"