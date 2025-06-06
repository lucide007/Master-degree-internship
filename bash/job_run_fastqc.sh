#!/bin/bash

#submit : sbatch job_run_fastqc.sh

#SBATCH --job-name=workflow_brb_tise_step2

#SBATCH --account=chronobiologie
#SBATCH --partition=cpucourt

#SBATCH --time=0-10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64GB

#SBATCH --output=/workspace/lbolelli/logs/2_%j.out
#SBATCH --error=/workspace/lbolelli/logs/2_%j.err

#SBATCH --mail-user=lucie.bolelli@etu.unice.fr
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniconda

# commands to run job
mkdir -p /workspace/lbolelli/logs
mkdir -p /workspace/lbolelli/data/fastqc_report/reads_1
mkdir -p /workspace/lbolelli/data/fastqc_report/reads_1

conda activate /home/lbolelli/.conda/envs/env_brb_tise

QC_DIR_1="/workspace/lbolelli/data/fastqc_report/reads_1"
QC_DIR_2="/workspace/lbolelli/data/fastqc_report/reads_2"

cd /workspace/lbolelli/data
fastqc --outdir $QC_DIR_1 merged_R1.fastq.gz
fastqc --outdir $QC_DIR_2 merged_R2.fastq.gz
echo "QC done"

conda deactivate
