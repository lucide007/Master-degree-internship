#!/bin/bash

#submit : sbatch job_align_reads.sh

#SBATCH --job-name=workflow_brb_tise_step3

#SBATCH --account=chronobiologie
#SBATCH --partition=cpucourt

#SBATCH --time=0-10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64GB

#SBATCH --output=/workspace/lbolelli/logs/3_%j.out
#SBATCH --error=/workspace/lbolelli/logs/3_%j.err

#SBATCH --mail-user=lucie.bolelli@etu.unice.fr
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniconda

# commands to run job
mkdir -p /workspace/lbolelli/logs
mkdir -p /workspace/lbolelli/data/alignment_results

conda activate /home/lbolelli/.conda/envs/env_brb_tise

BAM_DIR="/workspace/lbolelli/data/alignment_results"
STAR_INDEX_DIR="/workspace/lbolelli/genome/star_index"
FASTQ_DIR_MERGED="/workspace/lbolelli/data"
BARCODES_FILE="/workspace/lbolelli/genome/barcode_V5D_MERCURIUS.txt"

STAR --runMode alignReads --outSAMmapqUnique 60 --runThreadN 16 --outSAMunmapped Within \
    --soloStrand Forward --quantMode GeneCounts --outBAMsortingThreadN 8 \
    --genomeDir $STAR_INDEX_DIR --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 14 \
    --soloUMIstart 15 --soloUMIlen 14 --soloUMIdedup NoDedup 1MM_All --soloCellFilter None \
    --soloCBwhitelist $BARCODES_FILE --soloBarcodeReadLength 0 --soloFeatures Gene \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $BAM_DIR/ --readFilesIn $FASTQ_DIR_MERGED/merged_R2.fastq.gz $FASTQ_DIR_MERGED/merged_R1.fastq.gz

conda deactivate
