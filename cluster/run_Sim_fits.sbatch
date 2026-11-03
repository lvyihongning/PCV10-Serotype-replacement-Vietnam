#!/bin/bash
# ---------------------------------------------------------------
# High-carriage simulations (optional on the cluster; a few minutes
# per fit, so the Mac also works): the 20 fits as a Slurm array.
#   tasks  1-5   sim1, unadj, approaches 1-4 + presence-based
#   tasks  6-10  sim1, adj
#   tasks 11-15  sim2, unadj
#   tasks 16-20  sim2, adj
# Each task saves
# results/06_sensitivity/simulations/fits/Sim_fit_<sim>_<run>_a<k>_20261015.rds.
#
# Needs on the cluster: results/06_sensitivity/simulations/Sim_data_20261015.rds
# (built on the Mac by Sim_00_simulate_data_20261015.R).
#
# Submit after setup_pneumoalloc.sbatch:
#   fid=$(sbatch --parsable "cluster/run_Sim_fits.sbatch")
#   sbatch --dependency=afterok:$fid "cluster/run_Sim_combine.sbatch"
# Progress: bash "cluster/progress.sh" $fid Sim_fit
# ---------------------------------------------------------------
#SBATCH --job-name=Sim_fit
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/Sim_fit_%A_%a.out
#SBATCH --error=logs/Sim_fit_%A_%a.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

# Same order as `runs` in Sim_01: sim, then run type, then approach
SIMS=(sim1 sim2)
RUN_TYPES=(unadj adj)
METHODS=("dominant serotype" "equal share" "relative abundance"
         "serotype-specific density" "presence-based")
i=$((SLURM_ARRAY_TASK_ID - 1))
SIM=${SIMS[$((i / 10))]}
RUN=${RUN_TYPES[$(((i / 5) % 2))]}
METHOD=${METHODS[$((i % 5))]}

module load R/4.5.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

echo "Task $SLURM_ARRAY_TASK_ID: $SIM, $RUN, $METHOD, $SLURM_CPUS_PER_TASK CPUs on $(hostname)"

Rscript "06_sensitivity/simulations/Sim_01_fit_20261015.R" "$SLURM_ARRAY_TASK_ID"
