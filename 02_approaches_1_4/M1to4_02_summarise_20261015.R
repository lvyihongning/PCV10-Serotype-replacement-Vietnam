############################################################
# Approaches i-iv, unified framework — step 2: tables and checks
#
# Reads the fits of M1to4_01_fit_<STAMP>.R and writes, per fit, the
# workbook layout of the earlier analysis (read by 07_plots):
#   results/02_approaches_1_4/All_M<a>_<run>_re0_<STAMP>.xlsx
# and across approaches:
#   results/02_approaches_1_4/M1to4_cross_approach_re0_<STAMP>.xlsx
#     Shared_check_<run>  overall carriage, multiplicity and mean density
#                         must agree across approaches within Monte
#                         Carlo error (pneumoalloc::compare_approaches)
#     Class_all_<run>, Sero_all_<run>, Ratios_all_<run>
#
# "_re0" in the file names is kept for the figure scripts; there are
# no commune effects.
#
# Runs in a few minutes on the Mac (or on the cluster after the fits:
# cluster/run_M1to4_summarise.sbatch). Run from the project directory:
#   Rscript "02_approaches_1_4/M1to4_02_summarise_20261015.R"
############################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(openxlsx)

source("00_setup/settings_20261015.R")
source("00_setup/helpers_20261015.R")

out_dir <- "results/02_approaches_1_4"
fit_dir <- file.path(out_dir, "fits")
fit_file <- function(a, run_type) {
  file.path(fit_dir, paste0("M1to4_fit_M", a, "_", run_type, "_", STAMP, ".rds"))
}

ad <- load_analysis_data()

# ==========================================================
# 1) ONE FIT -> WORKBOOK
# ==========================================================
summarise_one <- function(fit, run_type) {
  method <- method_of(fit)
  vt <- fit$data$vaccine_types
  cats <- fit$labels$categories

  mp <- marginal_prevalence(fit)
  pr <- prevalence_ratios(fit)
  class_of <- c(Overall = "overall", VT = "VT", NVT = "NVT")

  # ---- prevalence ----
  P_class <- mp |>
    filter(group %in% c("overall", "VT", "NVT"), measure == "prevalence") |>
    legacy("class", "P") |>
    mutate(class = unname(class_of[class]), variable = paste0("P_", class),
           prevalence = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 2),
           method = method)

  P_sero <- mp |>
    filter(group == "serotype", measure == "prevalence") |>
    legacy("serotype", "P_sero") |>
    mutate(prevalence = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 2),
           sero_group = if_else(serotype %in% vt, "VT", "NVT"),
           method = method)

  # ---- prevalence ratios versus Oct. 2016 ----
  PR_class <- pr |>
    filter(group %in% c("overall", "VT", "NVT"), measure == "prevalence") |>
    legacy("class", "PR") |>
    mutate(class = unname(class_of[class]), variable = paste0("PR_", class),
           PR = fmt_ci(median, q2.5, q97.5, 2), method = method)

  baseline_prev <- P_sero |>
    filter(year == BASELINE) |>
    select(serotype, baseline_med = median, baseline_q97.5 = q97.5)
  PR_sero <- pr |>
    filter(group == "serotype", measure == "prevalence") |>
    legacy("serotype", "PR_sero") |>
    mutate(PR = fmt_ci(median, q2.5, q97.5, 2),
           sero_group = if_else(serotype %in% vt, "VT", "NVT"),
           method = method) |>
    left_join(baseline_prev, by = "serotype") |>
    mutate(PR_reportable = baseline_q97.5 >= 0.001,   # 0.1% threshold
           PR_note = if_else(PR_reportable, "",
                             "baseline near zero; report absolute difference instead"))

  PD_sero <- summarise_stan(fit, "PD_sero", cats) |>
    filter(year != BASELINE) |>
    mutate(PD_pp = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 2), method = method)

  # ---- composition among carriers ----
  share_class <- mp |>
    filter(group %in% c("VT", "NVT"), measure == "proportion among carriers") |>
    legacy("class", "share") |>
    mutate(variable = paste0("share_", class),
           Prop_among_carriers = fmt_ci(median, q2.5, q97.5, 3), method = method)
  share_sero <- mp |>
    filter(group == "serotype", measure == "proportion among carriers") |>
    legacy("serotype", "share") |>
    mutate(Prop_among_carriers = fmt_ci(median, q2.5, q97.5, 3),
           sero_group = if_else(serotype %in% vt, "VT", "NVT"),
           method = method)

  # ---- 15B/C + 15A + 23A among NVT (percent) ----
  gp <- group_proportions(fit, TARGET_GROUP, among = "NVT")
  prop_3 <- gp |>
    filter(measure == "proportion among NVT") |>
    legacy("serotype_group", "prop_15BC_15A_23A_among_NVT") |>
    mutate(across(c(mean, median, sd, q2.5, q97.5), ~ 100 * .x),
           proportion = fmt_ci(median, q2.5, q97.5, 1), method = method)
  prop_3_ratio <- gp |>
    filter(measure == "proportion ratio") |>
    legacy("serotype_group", "proportion_ratio") |>
    mutate(ratio = fmt_ci(median, q2.5, q97.5, 2), method = method)

  # ---- density and multiplicity ----
  md <- marginal_density(fit)
  meanD <- md |>
    filter(outcome == "Density") |>
    legacy("outcome", "meanD_pos") |>
    mutate(mean_density_GEml  = fmt_ci(median, q2.5, q97.5, 1),
           mean_density_log10 = fmt_ci(log10(median), log10(q2.5), log10(q97.5), 2))
  meanD_by_m <- md |>
    filter(outcome != "Density") |>
    legacy("n_serotypes", "meanD_pos_by_m") |>
    mutate(n_serotypes = sub("All carriers with ", "", n_serotypes),
           mean_density_GEml = fmt_ci(median, q2.5, q97.5, 1))
  dr <- density_ratios(fit, multiplicity = TRUE)
  PR_meanD <- dr |>
    filter(outcome == "Density") |>
    legacy("outcome", "PR_meanD_pos") |>
    mutate(density_ratio = fmt_ci(median, q2.5, q97.5, 2))
  dens_folds <- dr |>
    filter(outcome != "Density") |>
    transmute(comparison = outcome, across(any_of(SUMMARY_COLS)),
              estimate = fmt_ci(median, q2.5, q97.5, 2))

  P_mult <- mp |>
    filter(group == "multiplicity") |>
    legacy("n_serotypes", "P_mult") |>
    mutate(proportion = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 2))
  # observed proportions with 1, 2, 3 serotypes among carriers
  states <- describe_carriage_states(fit$data) |>
    filter(outcome %in% unique(P_mult$n_serotypes)) |>
    group_by(survey) |>
    mutate(obs_prop = 100 * count / sum(count)) |>
    ungroup() |>
    transmute(year = as.character(survey), n_serotypes = outcome, obs_prop)
  mult_check <- P_mult |>
    left_join(states, by = c("year", "n_serotypes")) |>
    mutate(model_prop = 100 * median, diff = model_prop - obs_prop)

  # ---- checks ----
  chk <- check_fit(fit)
  print(chk$summary)

  list(
    Run_metadata     = fit_metadata(fit, run_type,
                                    n_dropped_missing_data = ad$n_dropped,
                                    n_income_capped = ad$n_income_capped),
    Preva_Class_long = P_class,
    Preva_Sero_long  = P_sero,
    PR_Class_long    = PR_class,
    PR_Sero_long     = PR_sero,
    PD_Sero_long     = PD_sero,
    Share_Class_long = share_class,
    Share_Sero_long  = share_sero,
    Prop_3_among_NVT = prop_3,
    Prop_3_ratio     = prop_3_ratio,
    MeanD_pos_long   = meanD,
    PR_meanD_pos     = PR_meanD,
    MeanD_by_m       = meanD_by_m,
    Density_folds    = dens_folds,
    P_mult           = P_mult,
    Mult_check       = mult_check,
    Checks           = chk$summary,
    Identity_checks  = chk$identities,
    Crude_check      = chk$crude
  )
}

# ==========================================================
# 2) ALL FITS PRESENT
# ==========================================================
cross <- list()
for (run_type in RUN_TYPES) {
  files <- setNames(fit_file(1:4, run_type), APPROACHES)
  files <- files[file.exists(files)]
  if (length(files) == 0) next
  fits <- lapply(files, readRDS)

  for (nm in names(fits)) {
    a <- fits[[nm]]$approach$approach
    message("\n================ M", a, " (", nm, "), ", run_type, " ================")
    xlsx_file <- file.path(out_dir, paste0("All_M", a, "_", run_type, "_re0_", STAMP, ".xlsx"))
    write.xlsx(summarise_one(fits[[nm]], run_type), xlsx_file, overwrite = TRUE)
    message("Written: ", xlsx_file)
  }

  if (length(fits) >= 2) {
    ref <- if ("density" %in% names(fits)) "density" else names(fits)[1]
    cmp <- compare_approaches(fits, reference_approach = ref)
    message("\n[", run_type, "] shared quantities vs ", ref, ": max |z| = ",
            round(max(abs(cmp$shared$z), na.rm = TRUE), 2), "; flagged = ",
            sum(cmp$shared$flag, na.rm = TRUE), " of ", nrow(cmp$shared))
    cross[[paste0("Shared_check_", run_type)]] <- cmp$shared
    cross[[paste0("Class_all_", run_type)]] <- cmp$prevalence |> filter(group != "serotype")
    cross[[paste0("Sero_all_", run_type)]] <- cmp$prevalence |> filter(group == "serotype")
    cross[[paste0("Ratios_all_", run_type)]] <- cmp$ratios
  }
}

if (length(cross) > 0) {
  xlsx_file <- file.path(out_dir, paste0("M1to4_cross_approach_re0_", STAMP, ".xlsx"))
  write.xlsx(cross, xlsx_file, overwrite = TRUE)
  message("Written: ", xlsx_file)
}
