#!/bin/bash
# ---------------------------------------------------------------
# Approaches 1-4 x (unadj, adj): 8 independent fits as a Slurm array.
# Each array task runs ONE fit (02_approaches_1_4/M1to4_01_fit) with 4
# chains in parallel, one core per chain.
#
# Needs on the cluster: results/00_setup/analysis_data_20261015.rds
# (built on the Mac by 00_setup/01_prepare_data_20261015.R).
#
# Submit from the project directory, after setup_pneumoalloc.sbatch:
#   fid=$(sbatch --parsable --dependency=afterok:$sid "cluster/run_M1to4_fits.sbatch")
#   sbatch --dependency=afterok:$fid "cluster/run_M1to4_summarise.sbatch"
# Re-run some fits only, e.g. tasks 5 and 8:
#   sbatch --array=5,8 "cluster/run_M1to4_fits.sbatch"
# Progress: bash "cluster/progress.sh" $fid
# ---------------------------------------------------------------
#SBATCH --job-name=M1to4
#SBATCH --array=1-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/M1to4_%A_%a.out
#SBATCH --error=logs/M1to4_%A_%a.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

APPROACHES=(1 2 3 4 1 2 3 4)
RUN_TYPES=(unadj unadj unadj unadj adj adj adj adj)
i=$((SLURM_ARRAY_TASK_ID - 1))

module load R/4.5.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "Task $SLURM_ARRAY_TASK_ID: approach ${APPROACHES[$i]}, ${RUN_TYPES[$i]}, $SLURM_CPUS_PER_TASK CPUs on $(hostname)"

Rscript "02_approaches_1_4/M1to4_01_fit_20261015.R" "${APPROACHES[$i]}" "${RUN_TYPES[$i]}"
