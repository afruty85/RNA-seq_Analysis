#!/bin/bash
#SBATCH --job-name=download_assembly_reads
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH --array=1-27%3
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

# The data are from this study:
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189619
    # https://www.ncbi.nlm.nih.gov/bioproject/PRJNA783897

# DIRECTORIES
wkdir="/home/FCAM/afrutos/ISG/Final_project2025"
ACC_LIST="${wkdir}/SRR_Acc_ferret_List.txt"
raw="${wkdir}/data/raw"
sra_dir="${wkdir}/data/sra"

module load sratoolkit

mkdir -p "$raw" "$sra_dir"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_LIST")

# STEP 1: PREFETCH
prefetch "$SRR" --output-directory "$sra_dir"

# STEP 2: CONVERT TO FASTQ
fasterq-dump \
    "$sra_dir/$SRR/$SRR.sra" \
    --skip-technical \
    --split-files \
    --threads 8 \
    --progress \
    --outdir "$raw"

# STEP 3: COMPRESS IF EXISTS
if [ -f "$raw/${SRR}_1.fastq" ]; then gzip "$raw/${SRR}_1.fastq"; fi
if [ -f "$raw/${SRR}_2.fastq" ]; then gzip "$raw/${SRR}_2.fastq"; fi
if [ -f "$raw/${SRR}.fastq" ]; then gzip "$raw/${SRR}.fastq"; fi

echo "Done, reads are inside $raw"

