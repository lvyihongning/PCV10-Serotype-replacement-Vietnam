############################################################
# Approaches i-iv, unified framework — step 1: fit
#
# pneumoalloc::fit_carriage(): carriage (logistic, all children),
# serotype composition among carriers (fractional multinomial logit)
# and density among carriers (lognormal with multiplicity sub-model),
# fitted jointly, for each approach x run type:
#   approach 1  dominant serotype
#   approach 2  equal share
#   approach 3  relative abundance
#   approach 4  serotype-specific density (carrier weight
#               Theta_i / Theta_bar_t, Theta_bar_t = mean density among
#               carriers of the same survey round)
# All fits share the population, design, priors and sampler; only the
# allocation weights (and the carrier weight of approach 4) differ.
#
# Writes: results/02_approaches_1_4/fits/M1to4_fit_M<a>_<run>_<STAMP>.rds
# Next:   M1to4_02_summarise_<STAMP>.R (tables), 06_sensitivity/6A
#
# About 20 min per adjusted fit on the Mac (M4); recommended on the
# cluster, 8 fits in parallel (cluster/run_M1to4_fits.sbatch).
#
# Run from the project directory:
#   Rscript "02_approaches_1_4/M1to4_01_fit_20261015.R"                   all 8 fits
#   Rscript "02_approaches_1_4/M1to4_01_fit_20261015.R" <approach> <run>  one fit
#   (approach 1-4, run unadj or adj)
############################################################

rm(list = ls())

source("00_setup/settings_20261015.R")

approaches <- 1:4
run_types  <- RUN_TYPES
cl_args <- commandArgs(trailingOnly = TRUE)
if (length(cl_args) == 2) {
  approaches <- as.integer(cl_args[1])
  run_types  <- cl_args[2]
  stopifnot(approaches %in% 1:4, run_types %in% RUN_TYPES)
}

fit_dir <- dir_create("results/02_approaches_1_4/fits")
fit_file <- function(a, run_type) {
  file.path(fit_dir, paste0("M1to4_fit_M", a, "_", run_type, "_", STAMP, ".rds"))
}

dat <- load_analysis_data()$dat
mod <- carriage_model("unified")    # compiled once, then cached

for (run_type in run_types) {
  for (a in approaches) {
    message("\n================ M", a, " (", APPROACHES[a], "), ", run_type, " ================")
    fit <- do.call(fit_carriage, c(
      list(dat, approach = APPROACHES[a], covariates = covariates_for(run_type),
           priors = PRIORS, model = mod),
      SAMPLER
    ))
    print(check_fit(fit)$summary)
    save_fit(fit, fit_file(a, run_type))
    message("Written: ", fit_file(a, run_type))
  }
}
