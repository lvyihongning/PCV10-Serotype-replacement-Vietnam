#!/bin/bash
# ---------------------------------------------------------------
# 6A sensitivity analysis (Fig S7): computed from the 8 primary fits
# (same posterior; only the VT/NVT sums change), so no fits of its own.
# Run after ALL primary fits have succeeded:
#   sbatch --dependency=afterok:$fid "cluster/run_Sensi_6A.sbatch"
# ---------------------------------------------------------------
#SBATCH --job-name=Sensi_6A
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Sensi_6A_%j.out
#SBATCH --error=logs/Sensi_6A_%j.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

module load R/4.5.0

Rscript "06_sensitivity/6A/Sensi_6A_20261015.R"
