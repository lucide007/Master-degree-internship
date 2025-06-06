#!/bin/bash

# submit : sbatch job_prepare_genome.sh

# SBATCH --job-name=prepare_genome

# SBATCH --account=chronobiologie
# SBATCH --partition=cpucourt

# SBATCH --time=0-5:00:00
# SBATCH --nodes=1
# SBATCH --ntasks=1
# SBATCH --cpus-per-task=12
# SBATCH --mem=64GB

# SBATCH --output=/workspace/lbolelli/logs/job_%j.out
# SBATCH --error=/workspace/lbolelli/logs/job_%j.err

# SBATCH --mail-user=lucie.bolelli@etu.unice.fr
# SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load miniconda

# commands to run job
mkdir -p /workspace/lbolelli/logs
conda activate /home/lbolelli/.conda/envs/env_brb_tise
bash prepare_genome.sh
conda deactivate
