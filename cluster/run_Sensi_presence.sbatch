#!/bin/bash
# ---------------------------------------------------------------
# Presence-based sensitivity analysis (Fig S8): unadjusted and
# adjusted fits as a 2-task Slurm array, 4 chains in parallel each.
#   task 1  unadj
#   task 2  adj
# Each task writes
#   results/06_sensitivity/presence_based/fits/Sensi_presence_based_fit_<run>_20261015.rds
#   results/06_sensitivity/presence_based/Sensi_presence_based_<run>_20261015.xlsx
#
# Needs on the cluster: results/00_setup/analysis_data_20261015.rds
# and pneumoalloc installed (setup_pneumoalloc.sbatch).
#
# Submit from the project directory:
#   sbatch "cluster/run_Sensi_presence.sbatch"
# Progress: bash "cluster/progress.sh" <jobid> Sensi_presence
# ---------------------------------------------------------------
#SBATCH --job-name=Sensi_presence
#SBATCH --array=1-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Sensi_presence_%A_%a.out
#SBATCH --error=logs/Sensi_presence_%A_%a.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

RUN_TYPES=(unadj adj)
RUN=${RUN_TYPES[$((SLURM_ARRAY_TASK_ID - 1))]}

module load R/4.5.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "Task $SLURM_ARRAY_TASK_ID: presence-based, $RUN, $SLURM_CPUS_PER_TASK CPUs on $(hostname)"

Rscript "06_sensitivity/presence_based/Sensi_presence_based_20261015.R" "$RUN"
