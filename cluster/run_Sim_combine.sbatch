#!/bin/bash
# ---------------------------------------------------------------
# High-carriage simulations: combine the 20 fits (see run_Sim_fits.sbatch).
# ---------------------------------------------------------------
#SBATCH --job-name=Sim_combine
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Sim_combine_%j.out
#SBATCH --error=logs/Sim_combine_%j.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

module load R/4.5.0

Rscript "06_sensitivity/simulations/Sim_01_fit_20261015.R" combine
