#!/bin/bash
#SBATCH --job-name=ferretp5_htseq
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Quantify gene expression using HTseq-count

module load htseq

# Define directories
BASE_DIR="/home/FCAM/afrutos/ISG/Final_project2025"
ALIGN_DIR="${BASE_DIR}/results/alignments"
ANNOT_GTF="${BASE_DIR}/ferret/Mustela_putorius_furo.MusPutFur1.0.115.gtf"
COUNT_DIR="${BASE_DIR}/results/counts"

mkdir -p ${COUNT_DIR}

# Loop through BAM files and generate count tables
for BAM in ${ALIGN_DIR}/*_sorted.bam; do
    SAMPLE=$(basename "${BAM}" _sorted.bam)
    echo "Counting reads for sample: ${SAMPLE}"

    htseq-count \
        -f bam \
        -r pos \
        -s reverse \
        -t exon \
        -i gene_id \
        ${BAM} \
        ${ANNOT_GTF} \
        > ${COUNT_DIR}/${SAMPLE}.counts.txt
done

echo "HTSeq-count finished for all samples."
