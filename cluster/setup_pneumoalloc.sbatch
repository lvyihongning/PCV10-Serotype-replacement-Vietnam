#!/bin/bash
# ---------------------------------------------------------------
# Step 1: install pneumoalloc from the source tarball and compile its
# Stan models ONCE, so the array tasks do not all compile at the same
# time. Run again only when the package changes.
#
# Needs on the cluster (copied from the Mac; see README):
#   cluster/pneumoalloc_0.1.0.tar.gz
#
# Submit from the project directory:
#   mkdir -p logs
#   sid=$(sbatch --parsable "cluster/setup_pneumoalloc.sbatch")
# ---------------------------------------------------------------
#SBATCH --job-name=pa_setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/pa_setup_%j.out
#SBATCH --error=logs/pa_setup_%j.err
#SBATCH --partition=compute

set -eo pipefail   # no -u: the module command uses unset variables
cd "$SLURM_SUBMIT_DIR"

module load R/4.5.0
export MAKEFLAGS="-j${SLURM_CPUS_PER_TASK:-4}"

rm -rf "$HOME"/R/x86_64-pc-linux-gnu-library/4.5/00LOCK*
R CMD INSTALL "cluster/pneumoalloc_0.1.0.tar.gz"
Rscript -e 'for (m in c("unified", "unified_multiple", "categorical", "presence", "density")) pneumoalloc::carriage_model(m)' \
        -e 'cat("pneumoalloc", format(packageVersion("pneumoalloc")), "ready\n")'
