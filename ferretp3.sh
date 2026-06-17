#!/bin/bash
#SBATCH --job-name=ferretp3_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-user=alex.frutos@uconn.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Build HISAT2 index for ferret genome

# Load HISAT2
module load hisat2

# Set paths
REF_DIR="/home/FCAM/afrutos/ISG/Final_project2025/ferret"
GENOME_FA="$REF_DIR/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa"
OUT_PREFIX="$REF_DIR/MusPutFur1.0_index"

# Build index
hisat2-build -p 8 $GENOME_FA $OUT_PREFIX
