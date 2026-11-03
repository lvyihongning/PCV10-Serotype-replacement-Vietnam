############################################################
# 06_sensitivity/presence_based/Sensi_presence_based_20261015.R
#
# Presence-based sensitivity analysis (Supplement; Fig S8).
#
# pneumoalloc::fit_presence(): each carriage outcome is a child-level
# binary indicator, defined independently of dominance, relative
# abundance or density, with its own logistic model:
#   overall   >= 1 serotype detected
#   VT        >= 1 PCV10 serotype detected
#   NVT       >= 1 non-PCV10 serotype detected
#   multiple  >= 2 serotypes detected
#   serotype  that serotype detected ("Other NVT": any non-retained NVT)
# VT and NVT are therefore not constrained to sum to overall carriage.
#
# Same as the primary analysis: analysis population, covariates,
# serotype categories, design ~ 0 + survey (+ covariates), prior sd 2.5
# on the logit scale, unadjusted and adjusted runs.
#
# Output (results/06_sensitivity/presence_based/):
#   fits/Sensi_presence_based_fit_<run>_<STAMP>.rds
#   Sensi_presence_based_<unadj|adj>_<STAMP>.xlsx  (layout read by 07_plots)
#
# Roughly 10-20 min per fit on the Mac. Run from the project directory:
#   Rscript "06_sensitivity/presence_based/Sensi_presence_based_20261015.R" [unadj|adj]
############################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(openxlsx)

source("00_setup/settings_20261015.R")
source("00_setup/helpers_20261015.R")

cl_args <- commandArgs(trailingOnly = TRUE)
run_types <- if (length(cl_args) >= 1) cl_args else RUN_TYPES
stopifnot(run_types %in% RUN_TYPES)

method_lab <- "Presence-based"
out_dir <- dir_create("results/06_sensitivity/presence_based")
fit_dir <- dir_create(out_dir, "fits")

dat <- load_analysis_data()$dat
vt  <- dat$vaccine_types
class_of <- c(Overall = "overall", VT = "VT", NVT = "NVT", `Multiple carriage` = "multiple")

# Crude proportions of children with each outcome
crude <- describe_prevalence(dat, "presence")

for (run_type in run_types) {
  message("\n================ Presence-based, ", run_type, " ================")
  fit <- do.call(fit_presence,
                 c(list(dat, covariates = covariates_for(run_type), priors = PRIORS), SAMPLER))
  save_fit(fit, file.path(fit_dir, paste0("Sensi_presence_based_fit_", run_type, "_", STAMP, ".rds")))
  chk <- check_fit(fit)
  print(chk$summary)

  mp <- marginal_prevalence(fit)
  pr <- prevalence_ratios(fit)

  preva_class <- mp |>
    filter(group != "serotype") |>
    legacy("class", "P") |>
    mutate(class = unname(class_of[class]), variable = paste0("P_", class),
           prevalence = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5), method = method_lab)

  preva_sero <- mp |>
    filter(group == "serotype") |>
    legacy("serotype", "P") |>
    mutate(variable = paste0("P_", serotype),
           prevalence = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5),
           sero_group = if_else(serotype %in% vt, "VT", "NVT"), method = method_lab)

  pr_class <- pr |>
    filter(group != "serotype") |>
    legacy("class", "PR") |>
    mutate(class = unname(class_of[class]), variable = paste0("PR_", class),
           PR = fmt_ci(median, q2.5, q97.5), method = method_lab)

  # Flag serotypes whose PR is unstable because baseline is near zero,
  # with the same 0.1% threshold as the primary analysis
  baseline_prev <- preva_sero |>
    filter(year == BASELINE) |>
    select(serotype, baseline_med = median, baseline_q97.5 = q97.5)
  pr_sero <- pr |>
    filter(group == "serotype") |>
    legacy("serotype", "PR") |>
    mutate(variable = paste0("PR_", serotype),
           PR = fmt_ci(median, q2.5, q97.5),
           sero_group = if_else(serotype %in% vt, "VT", "NVT"), method = method_lab) |>
    left_join(baseline_prev, by = "serotype") |>
    mutate(PR_reportable = baseline_q97.5 >= 0.001,
           PR_note = if_else(PR_reportable, "",
                             "baseline near zero; report absolute difference instead"))

  write.xlsx(
    list(
      Run_metadata     = fit_metadata(fit, run_type),
      Preva_Class_long = preva_class,
      Preva_Sero_long  = preva_sero,
      PR_Class_long    = pr_class,
      PR_Sero_long     = pr_sero,
      Crude_check      = chk$crude,
      Presence_check   = crude,
      Checks           = chk$summary
    ),
    file = file.path(out_dir, paste0("Sensi_presence_based_", run_type, "_", STAMP, ".xlsx")),
    overwrite = TRUE
  )

  # Key manuscript numbers: VT and NVT, October 2016 and October 2019
  print(preva_class |> filter(year %in% c(BASELINE, "Oct. 2019")) |> select(year, class, prevalence))
  print(pr_class |> filter(year == "Oct. 2019") |> select(year, class, PR))
}

message("\nDone. Written to ", out_dir, "/.")

