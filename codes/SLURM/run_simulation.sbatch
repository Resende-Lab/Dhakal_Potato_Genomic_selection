#!/bin/bash
#SBATCH --job-name=PotatoSim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10GB
#SBATCH --time=96:00:00
#SBATCH --mail-type=END
#SBATCH --output=logs/sim_%a_%A.out
#SBATCH --error=logs/sim_%a_%A.err
#SBATCH --array=1-30

module purge
module load R

# Each line of rep.txt: <rep_number> <ddMeanSG>
INPUT=$(head -n $SLURM_ARRAY_TASK_ID rep.txt | tail -n 1)

echo "Running rep: $INPUT"
Rscript RUNME.R $INPUT
