#!/bin/bash
#SBATCH --job-name=GS_CV0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --nodes=1
#SBATCH --mem=250GB
#SBATCH -t 96:00:00
#SBATCH --output=logs/GS_CV0_%j.out
#SBATCH --error=logs/GS_CV0_%j.err


# ── Load modules ──
module purge
module load R

# ── Run ──
Rscript GS_CV0_LOEO.R
