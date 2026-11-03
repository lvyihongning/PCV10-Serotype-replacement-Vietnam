#!/bin/bash
# ---------------------------------------------------------------
# Tables and cross-approach checks of the 8 fits; run after ALL fits
# have succeeded (see run_M1to4_fits.sbatch). Can also be run on the
# Mac after copying results/02_approaches_1_4/fits/ back.
# ---------------------------------------------------------------
#SBATCH --job-name=M1to4_sum
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/M1to4_sum_%j.out
#SBATCH --error=logs/M1to4_sum_%j.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

module load R/4.5.0

Rscript "02_approaches_1_4/M1to4_02_summarise_20261015.R"
