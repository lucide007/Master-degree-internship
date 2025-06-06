#!/bin/bash

#submit : sbatch job_merge_fastq.sh

#SBATCH --job-name=workflow_brb_tise_step1

#SBATCH --account=chronobiologie
#SBATCH --partition=cpucourt

#SBATCH --time=0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64GB

#SBATCH --output=/workspace/lbolelli/logs/1_%j.out
#SBATCH --error=/workspace/lbolelli/logs/1_%j.err

#SBATCH --mail-user=lucie.bolelli@etu.unice.fr
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniconda

# commands to run job
mkdir -p /workspace/lbolelli/logs
mkdir -p /workspace/lbolelli/raw_data
mkdir -p /workspace/lbolelli/data

conda activate /home/lbolelli/.conda/envs/env_brb_tise

FASTQ_DIR="/workspace/lbolelli/raw_data"
FASTQ_DIR_MERGED="/workspace/lbolelli/data"

echo "Merging FASTQ files..."
cat $FASTQ_DIR/*R1*.fastq.gz > $FASTQ_DIR_MERGED/merged_R1.fastq.gz
cat $FASTQ_DIR/*R2*.fastq.gz > $FASTQ_DIR_MERGED/merged_R2.fastq.gz
echo "FASTQ files merged"

conda deactivate

