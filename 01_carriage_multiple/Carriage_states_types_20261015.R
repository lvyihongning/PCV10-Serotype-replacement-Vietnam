############################################################
# Overall and multiple serotype carriage, and types of multiple
# carriage — unadjusted and adjusted (Fig 1; Tables S2, S6, S7)
#
# 1. Carriage states: number of detected serotypes (0, 1, 2, 3),
#    multinomial logit, pneumoalloc::fit_carriage_states(). Gives
#    overall carriage (>= 1) and multiple carriage (>= 2) prevalence
#    and prevalence ratios versus Oct. 2016.
# 2. Types of multiple carriage among children with >= 2 serotypes
#    (VT-only, NVT-only, mixed VT-dominant, mixed NVT-dominant),
#    multinomial logit, pneumoalloc::fit_multiple_types(). The type is
#    built from the detected serotypes and the dominant serotype
#    (highest relative abundance); 00_setup checks it against the raw
#    domin_subdomin_sero_1 pattern.
# Survey-specific intercepts (+ covariates in the adjusted model);
# g-computation over the analysis population.
#
# Writes (results/01_carriage_multiple/), same layout as before:
#   0_co-colonisation_prevalence_PR_<run>_<STAMP>.xlsx
#   0_co-colonisation_proportion_PR_<run>_<STAMP>.xlsx
#   fits/Carriage_states_<run>_<STAMP>.rds, fits/Multiple_types_<run>_<STAMP>.rds
# Figure 1: 07_plots/Overall_Multi_Preva_Prop_Plot_<STAMP>.R
#
# Run from the project directory (a few minutes per fit on the Mac):
#   Rscript "01_carriage_multiple/Carriage_states_types_20261015.R" [unadj|adj]
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

out_dir <- dir_create("results/01_carriage_multiple")
fit_dir <- dir_create(out_dir, "fits")

dat <- load_analysis_data()$dat

# Category labels of the earlier result files
CAT_LABEL <- c("No carriage" = "m0", "1 serotype" = "m1", "2 serotypes" = "m2",
               "3 serotypes" = "m3", "3+ serotypes" = "m3",
               "Carriage" = "m>=1 (carriage)", "Multiple carriage" = "m>=2 (co-colonisation)")
TYPE_LEVELS <- c("VT-only", "NVT-only", "Mixed, VT-dominant", "Mixed, NVT-dominant")

# One decimal, "-" separator, as in the earlier tables
fmt_pct <- function(med, lo, hi) fmt_ci(100 * med, 100 * lo, 100 * hi, 1, sep = "-")

crude <- describe_carriage_states(dat) |>
  mutate(year = as.character(survey))

# ==========================================================
# 1) CARRIAGE STATES
# ==========================================================
run_states <- function(run_type) {
  message("\n================ Carriage states, ", run_type, " ================")
  fit <- do.call(fit_carriage_states,
                 c(list(dat, covariates = covariates_for(run_type), priors = PRIORS), SAMPLER))
  save_fit(fit, file.path(fit_dir, paste0("Carriage_states_", run_type, "_", STAMP, ".rds")))
  chk <- check_fit(fit)
  print(chk$summary)

  prev_all <- marginal_prevalence(fit) |>
    legacy("category", "P_class") |>
    mutate(category = unname(CAT_LABEL[category]),
           group = RUN_LABEL[[run_type]],
           prevalence = fmt_pct(median, q2.5, q97.5)) |>
    select(-variable)

  pr_all <- prevalence_ratios(fit) |>
    legacy("category", "PR") |>
    mutate(category = unname(CAT_LABEL[category]),
           group = RUN_LABEL[[run_type]],
           prevalence_ratio = fmt_ci(median, q2.5, q97.5, 2, sep = "-")) |>
    select(-variable)

  crude_vs_model <- prev_all |>
    select(year, category, prevalence, median) |>
    left_join(crude |>
                filter(measure == "prevalence") |>
                transmute(year, category = unname(CAT_LABEL[outcome]), crude_prev = 100 * value),
              by = c("year", "category")) |>
    mutate(diff_pp = 100 * median - crude_prev)
  message("Max |model - crude| = ", round(max(abs(crude_vs_model$diff_pp)), 3), " pp")

  write.xlsx(
    list(preva_colonistation = prev_all,
         pr_all              = pr_all,
         preva_table         = prev_all |> select(year, group, category, prevalence) |>
                                 pivot_wider(names_from = year, values_from = prevalence),
         pr_table            = pr_all |> select(year, group, category, prevalence_ratio) |>
                                 pivot_wider(names_from = year, values_from = prevalence_ratio),
         crude_vs_model      = crude_vs_model,
         diagnostics         = fit_metadata(fit, run_type),
         checks              = chk$summary),
    file = file.path(out_dir, paste0("0_co-colonisation_prevalence_PR_", run_type, "_", STAMP, ".xlsx")),
    overwrite = TRUE
  )
}

# ==========================================================
# 2) TYPES OF MULTIPLE CARRIAGE
# ==========================================================
run_types_model <- function(run_type) {
  message("\n================ Types of multiple carriage, ", run_type, " ================")
  # adapt_delta 0.95 (package default for this model), as before
  fit <- do.call(fit_multiple_types,
                 c(list(dat, covariates = covariates_for(run_type), priors = PRIORS), SAMPLER))
  save_fit(fit, file.path(fit_dir, paste0("Multiple_types_", run_type, "_", STAMP, ".rds")))
  chk <- check_fit(fit)
  print(chk$summary)

  prev_class <- marginal_prevalence(fit) |>
    legacy("category", "prop") |>
    mutate(category = factor(category, levels = TYPE_LEVELS),
           group = RUN_LABEL[[run_type]],
           prevalence = fmt_pct(median, q2.5, q97.5)) |>
    arrange(year, category) |>
    mutate(category = as.character(category)) |>
    select(-variable)

  pr <- prevalence_ratios(fit) |> legacy("category", "PR")
  pr_class <- pr |>
    mutate(group = RUN_LABEL[[run_type]],
           prevalence_ratio = fmt_ci(median, q2.5, q97.5, 2, sep = "-")) |>
    select(-variable)
  # log ratio and percent change: transformed posterior quantiles
  logpr_class <- pr |>
    transmute(year, category, median = log(median), q2.5 = log(q2.5), q97.5 = log(q97.5),
              group = RUN_LABEL[[run_type]])
  pct_class <- pr |>
    transmute(year, category, median = 100 * (median - 1), q2.5 = 100 * (q2.5 - 1),
              q97.5 = 100 * (q97.5 - 1), group = RUN_LABEL[[run_type]]) |>
    mutate(pct_change_fmt = paste0(fmt_ci(median, q2.5, q97.5, 1), "%"))

  prev_comparison <- prev_class |>
    select(year, category, median, q2.5, q97.5, prevalence, group) |>
    left_join(crude |>
                filter(measure == "proportion among multiple carriers") |>
                transmute(year, category = outcome, n = count, prop_crude = value),
              by = c("year", "category")) |>
    mutate(diff = median - prop_crude)
  message("Max |model - crude| = ", round(100 * max(abs(prev_comparison$diff)), 2), " pp")

  write.xlsx(
    list(prev_class      = prev_class,
         pr_class        = pr_class,
         logpr_class     = logpr_class,
         pct_class       = pct_class,
         prev_table      = prev_class |> select(year, group, category, prevalence) |>
                             pivot_wider(names_from = year, values_from = prevalence),
         pr_table        = pr_class |> select(year, group, category, prevalence_ratio) |>
                             pivot_wider(names_from = year, values_from = prevalence_ratio),
         prev_comparison = prev_comparison,
         diagnostics     = fit_metadata(fit, run_type),
         checks          = chk$summary),
    file = file.path(out_dir, paste0("0_co-colonisation_proportion_PR_", run_type, "_", STAMP, ".xlsx")),
    overwrite = TRUE
  )
}

for (rt in run_types) {
  run_states(rt)
  run_types_model(rt)
}
message("Done: ", paste(run_types, collapse = ", "))
