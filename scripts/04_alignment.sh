#!/bin/bash
#SBATCH --job-name=ferretp4_align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# RNA-seq Final Project 2025
# Align all ferret FASTQ files with HISAT2

# Load modules
module load hisat2
module load samtools

# Define directories
BASE_DIR="/home/FCAM/afrutos/ISG/Final_project2025"
RAW_DIR="${BASE_DIR}/data/raw"
RESULTS_DIR="${BASE_DIR}/results/alignments"
INDEX_DIR="${BASE_DIR}/ferret/MusPutFur1.0_index"

# Make output directory
mkdir -p $RESULTS_DIR

# Align all paired-end samples
for READ1 in ${RAW_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename "${READ1}" _1.fastq.gz)
    echo "Processing sample: ${SAMPLE}"

    READ2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"
    SAM_OUT="${RESULTS_DIR}/${SAMPLE}.sam"
    BAM_OUT="${RESULTS_DIR}/${SAMPLE}_sorted.bam"
    LOG_OUT="${RESULTS_DIR}/${SAMPLE}_hisat2.log"

    # Align reads
    hisat2 -p 8 \
        --dta \
        --rna-strandness RF \
        -x ${INDEX_DIR} \
        -1 ${READ1} \
        -2 ${READ2} \
        -S ${SAM_OUT} \
        2> ${LOG_OUT}

    # Convert to sorted BAM
    samtools sort -@ 8 -o ${BAM_OUT} ${SAM_OUT}
    samtools index ${BAM_OUT}

    # Remove SAM file to save space
    rm ${SAM_OUT}

    echo "Completed sample: $SAMPLE"
done

echo "Alignments complete."
