#!/bin/bash
#SBATCH --job-name=GS_CV1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --nodes=1
#SBATCH --mem=250GB
#SBATCH -t 96:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your-email>
#SBATCH --qos=<your-qos>
#SBATCH --account=<your-account>
#SBATCH --output=logs/GS_CV1_%j.out
#SBATCH --error=logs/GS_CV1_%j.err

# ── Load modules ──
module purge
module load R

# ── Run ──
Rscript GS_CV1_GBLUP.R
