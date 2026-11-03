############################################################
# 05_density/Density_20261015.R
#
# Pneumococcal density among carriage-positive samples
# (Figure 3, Tables S6 and S9).
#
# pneumoalloc::fit_density(): lognormal regression on the natural-log
# density scale with nested multiplicity effects (>= 2 and >= 3
# detected serotypes), fitted separately for all pneumococci, PCV10
# (VT) and non-PCV10 (NVT) serotypes. VT and NVT density is the sum of
# the serotype-specific densities (total density x relative abundance)
# of the group's serotypes; their multiplicity effect counts the
# serotypes within the group (multiplicity = "group"). Multiplicity
# effects are unconstrained.
#
# Four configurations:
#   unadj, adj                          primary density model
#                                       (Table S9, Fig 3)
#   interaction_unadj, interaction_adj  survey x multiplicity interaction
#                                       (Table S6)
#
# Output (results/05_density/):
#   fits/Density_fit_<config>_<STAMP>.rds     the three group fits
#   Density_PR_<config>_<STAMP>.xlsx          one workbook per configuration
#   Density_summaries_<config>_<STAMP>.rds    summaries used for the tables
#                                             (and by the simulations)
#   TableS6_Density ratios_<STAMP>.xlsx
#   TableS9_Median density ratios_<STAMP>.xlsx
#
# 07_plots/All_density_plot_<STAMP>.R reads Density_PR_adj_<STAMP>.xlsx.
# 12 fits; roughly 5-10 min each on the Mac.
#
# Run from the project directory:
#   Rscript "05_density/Density_20261015.R"                 all configurations
#   Rscript "05_density/Density_20261015.R" <config> ...    some, e.g. adj
#   Rscript "05_density/Density_20261015.R" tables          Tables S6/S9 only
############################################################

rm(list = ls())
library(tidyverse)
library(posterior)
library(openxlsx)

source("00_setup/settings_20261015.R")


############################################################
#### 0) Settings ###########################################
############################################################

configs <- list(
  unadj             = list(run_type = "unadj", interaction = FALSE),
  adj               = list(run_type = "adj",   interaction = FALSE),
  interaction_unadj = list(run_type = "unadj", interaction = TRUE),
  interaction_adj   = list(run_type = "adj",   interaction = TRUE)
)
configs_to_run <- names(configs)
cl_args <- commandArgs(trailingOnly = TRUE)
if (length(cl_args) >= 1) {
  configs_to_run <- if (identical(cl_args, "tables")) character() else cl_args
  stopifnot(configs_to_run %in% names(configs))
}

mult_source <- "group"       # VT/NVT: serotypes counted within the group
m_rep <- 4000L               # predictive draws per round (median, quantiles)
baseline_level <- BASELINE

group_labels <- c(total = "Overall", VT = "VT", NVT = "NVT")

out_dir <- dir_create("results/05_density")
fit_dir <- dir_create(out_dir, "fits")


############################################################
#### 1) Functions ##########################################
############################################################

# ---------------------------------------------------------------------
# Small formatting helpers
# ---------------------------------------------------------------------

fmt_ci <- function(med, lo, hi, digits = 2) {
  paste0(
    formatC(med, digits = digits, format = "f"), " (",
    formatC(lo, digits = digits, format = "f"), "-",
    formatC(hi, digits = digits, format = "f"), ")"
  )
}

safe_log10 <- function(x) {
  ifelse(x > 0, log10(x), NA_real_)
}

fmt_ci_log10 <- function(med, lo, hi, digits = 2) {
  fmt_ci(safe_log10(med), safe_log10(lo), safe_log10(hi), digits = digits)
}


# ---------------------------------------------------------------------
# Tidy a posterior rvar into a data frame with named dimensions
#
# Indices are taken from the bracketed part of the name only; parsing
# every digit in the name would also match the digits inside names such
# as q025D and mislabel the survey.
# ---------------------------------------------------------------------

tidy_rvar <- function(rv, name, dim_labels = NULL) {
  summ <- posterior::summarise_draws(
    posterior::draws_rvars(x = rv),
    posterior::default_summary_measures(),
    posterior::default_convergence_measures(),
    q2.5  = ~ posterior::quantile2(.x, probs = 0.025),
    q97.5 = ~ posterior::quantile2(.x, probs = 0.975)
  )

  idx_txt <- stringr::str_match(summ$variable, "^x\\[(.*)\\]$")[, 2]

  if (!all(is.na(idx_txt))) {
    idx_list <- stringr::str_split(idx_txt, ",")
    n_dim <- max(lengths(idx_list))
    idx_mat <- do.call(
      rbind,
      purrr::map(idx_list, \(z) as.integer(trimws(z)))
    )
    colnames(idx_mat) <- paste0("i", seq_len(n_dim))
    summ <- dplyr::bind_cols(summ, tibble::as_tibble(idx_mat))

    if (!is.null(dim_labels)) {
      summ <- purrr::reduce(
        seq_along(dim_labels),
        function(acc, d) {
          acc[[names(dim_labels)[d]]] <- dim_labels[[d]][acc[[paste0("i", d)]]]
          acc
        },
        .init = summ
      )
    }
  }

  summ |>
    dplyr::mutate(variable = name) |>
    tibble::as_tibble()
}


# ---------------------------------------------------------------------
# What a pneumoalloc density fit holds about its data: survey rounds,
# multiplicity levels, the observed density and number of serotypes
# ---------------------------------------------------------------------

fit_prep <- function(fit) {
  list(
    time_levels = fit$data$survey_levels,
    baseline_level = fit$data$reference_level,
    baseline = fit$data$reference,
    m_max = fit$labels$m_levels,
    dat = tibble::tibble(
      survey_yearmon = fit$data$survey[fit$rows],
      m_sh = fit$labels$m,
      dens_tot = fit$labels$density
    )
  )
}


# ---------------------------------------------------------------------
# Post-processing (generated quantities of the package's density.stan)
# ---------------------------------------------------------------------

density_vars_to_pull <- c(
  "meanD", "medianD", "q025D", "q25D", "q75D", "q975D", "IQRD",
  "meanD_by_m", "medianD_by_m", "q025D_by_m", "q25D_by_m",
  "q75D_by_m", "q975D_by_m", "IQRD_by_m",
  "PR_meanD", "PR_medianD", "PR_shiftD", "shift_maxdev",
  "PR_meanD_by_m", "PR_medianD_by_m",
  "sigma_dens", "gamma",
  "dens_fold_vs1", "dens_fold_stepwise", "RR_dens_fold_vs1"
)

m_labels_for <- function(m_max) {
  c("One serotype", "Two serotypes", "Three serotypes")[seq_len(m_max)]
}


summarise_density_fit <- function(fit, prep, group_label) {

  time_levels <- prep$time_levels
  baseline_level <- prep$baseline_level
  baseline <- prep$baseline
  m_labels <- m_labels_for(prep$m_max)

  d_rv <- posterior::as_draws_rvars(fit$draws(density_vars_to_pull))

  dims_t <- list(year = time_levels)
  dims_tm <- list(year = time_levels, m = m_labels)

  tidy_t <- function(nm) tidy_rvar(d_rv[[nm]], nm, dims_t)
  tidy_tm <- function(nm) {
    tidy_rvar(d_rv[[nm]], nm, dims_tm) |>
      dplyr::rename(m_id = i2)
  }

  # ---- overall density by survey round ----
  overall_density <- purrr::map(
    c("medianD", "q025D", "q25D", "q75D", "q975D", "IQRD", "meanD"),
    tidy_t
  ) |>
    purrr::list_rbind() |>
    dplyr::select(variable, year, median, q2.5, q97.5) |>
    tidyr::pivot_wider(
      names_from = variable,
      values_from = c(median, q2.5, q97.5)
    ) |>
    dplyr::mutate(
      year = factor(year, levels = time_levels),
      group = group_label,
      median_log10 = safe_log10(median_medianD),
      q25_log10 = safe_log10(median_q25D),
      q75_log10 = safe_log10(median_q75D),
      q025_log10 = safe_log10(median_q025D),
      q975_log10 = safe_log10(median_q975D),
      median_IQR_log10 = fmt_ci(median_log10, q25_log10, q75_log10),
      median_95range_log10 = fmt_ci(median_log10, q025_log10, q975_log10),
      median_95CrI_log10 = fmt_ci_log10(median_medianD, q2.5_medianD, q97.5_medianD),
      mean_ci_GEml = fmt_ci(median_meanD, q2.5_meanD, q97.5_meanD, digits = 1)
    ) |>
    dplyr::arrange(year)

  # ---- density by counterfactual multiplicity ----
  density_by_m <- purrr::map(
    c("medianD_by_m", "q025D_by_m", "q25D_by_m", "q75D_by_m",
      "q975D_by_m", "IQRD_by_m", "meanD_by_m"),
    tidy_tm
  ) |>
    purrr::list_rbind() |>
    dplyr::select(variable, year, m_id, m, median, q2.5, q97.5) |>
    tidyr::pivot_wider(
      names_from = variable,
      values_from = c(median, q2.5, q97.5)
    ) |>
    dplyr::mutate(
      year = factor(year, levels = time_levels),
      m = factor(m, levels = m_labels),
      group = group_label,
      median_log10 = safe_log10(median_medianD_by_m),
      q25_log10 = safe_log10(median_q25D_by_m),
      q75_log10 = safe_log10(median_q75D_by_m),
      median_IQR_log10 = fmt_ci(median_log10, q25_log10, q75_log10),
      median_95CrI_log10 = fmt_ci_log10(
        median_medianD_by_m, q2.5_medianD_by_m, q97.5_medianD_by_m
      )
    ) |>
    dplyr::arrange(m_id, year)

  # ---- density ratios versus baseline ----
  pr_overall <- purrr::map(
    c("PR_meanD", "PR_medianD", "PR_shiftD"),
    tidy_t
  ) |>
    purrr::list_rbind() |>
    dplyr::mutate(m = "Overall", m_id = NA_integer_)

  pr_by_m <- purrr::map(
    c("PR_meanD_by_m", "PR_medianD_by_m"),
    tidy_tm
  ) |>
    purrr::list_rbind()

  pr_density <- dplyr::bind_rows(pr_overall, pr_by_m) |>
    dplyr::filter(as.character(year) != baseline_level) |>
    dplyr::mutate(
      group = group_label,
      year = factor(year, levels = time_levels),
      PR = fmt_ci(median, q2.5, q97.5, digits = 2)
    ) |>
    dplyr::select(group, variable, year, m, m_id, median, q2.5, q97.5,
                  PR, rhat, ess_bulk, ess_tail)

  # ---- shift diagnostic ----
  # shift_maxdev near zero confirms that the covariate part of the survey
  # effect is a pure additive shift on the log scale, in which case
  # PR_shiftD is the exact quantile ratio and PR_medianD should agree
  # with it. With the interaction the multiplicity effects also change
  # by survey, so PR_shiftD no longer equals the quantile ratio.
  shift_diag <- tidy_t("shift_maxdev") |>
    dplyr::mutate(
      group = group_label,
      year = factor(year, levels = time_levels)
    ) |>
    dplyr::select(group, year, median, q97.5)

  # ---- multiplicity effects (Table S6 contrasts) ----
  # Constant across survey rounds without the interaction. Groups whose
  # carriers never have three serotypes of the group have only the
  # two-versus-one contrast.
  folds <- list("Two versus single serotype" = d_rv$dens_fold_vs1[, 2])
  if (prep$m_max >= 3) {
    folds[["Three versus single serotype"]] <- d_rv$dens_fold_vs1[, 3]
    folds[["Three versus two serotypes"]] <- d_rv$dens_fold_stepwise[, 3]
  }
  ratio_of_ratios <- lapply(folds, function(f) f / f[baseline])

  multiplicity_effects <- dplyr::bind_rows(
    purrr::imap(folds, \(f, nm) tidy_rvar(f, nm, dims_t)) |>
      purrr::list_rbind() |>
      dplyr::mutate(estimand = "Density ratio"),
    purrr::imap(ratio_of_ratios, \(f, nm) tidy_rvar(f, nm, dims_t)) |>
      purrr::list_rbind() |>
      dplyr::mutate(estimand = "Ratio of density ratios")
  ) |>
    dplyr::mutate(
      group = group_label,
      year = factor(year, levels = time_levels),
      estimate = fmt_ci(median, q2.5, q97.5, digits = 2)
    ) |>
    dplyr::select(group, estimand, contrast = variable, year,
                  median, q2.5, q97.5, estimate, rhat, ess_bulk, ess_tail) |>
    dplyr::arrange(estimand, contrast, year)

  # ---- raw parameters ----
  parameters <- posterior::summarise_draws(
    fit$draws(c("sigma_dens", "gamma")),
    posterior::default_summary_measures(),
    posterior::default_convergence_measures(),
    q2.5 = ~ posterior::quantile2(.x, probs = 0.025),
    q97.5 = ~ posterior::quantile2(.x, probs = 0.975)
  ) |>
    dplyr::mutate(
      group = group_label,
      estimate = fmt_ci(median, q2.5, q97.5, digits = 3)
    ) |>
    tibble::as_tibble()

  # ---- observed densities, for plotting ----
  observed <- tibble::tibble(
    group = group_label,
    year = factor(prep$dat$survey_yearmon, levels = time_levels),
    m_id = prep$dat$m_sh,
    m = factor(m_labels[prep$dat$m_sh], levels = m_labels),
    density = prep$dat$dens_tot,
    log10_density = log10(prep$dat$dens_tot)
  )

  list(
    overall_density = overall_density,
    density_by_m = density_by_m,
    pr_density = pr_density,
    shift_diagnostics = shift_diag,
    multiplicity_effects = multiplicity_effects,
    parameters = parameters,
    observed_density = observed
  )
}


# ---------------------------------------------------------------------
# Sampler diagnostics, from pneumoalloc::check_fit()
# ---------------------------------------------------------------------

sampler_diagnostics <- function(fit, group_label) {
  chk <- check_fit(fit)$summary
  v <- function(name) chk$value[chk$check == name]
  tibble::tibble(
    group = group_label,
    divergences = v("divergent transitions"),
    max_treedepth = v("max treedepth hits"),
    min_ebfmi = v("min E-BFMI"),
    max_rhat = v("max R-hat"),
    min_ess_bulk = v("min bulk ESS"),
    min_ess_tail = v("min tail ESS")
  )
}


# ---------------------------------------------------------------------
# Stack one summary table across the three serotype groups
# ---------------------------------------------------------------------

collect <- function(results, field) {
  purrr::map(results, \(r) r$summaries[[field]]) |>
    purrr::list_rbind() |>
    dplyr::mutate(group = factor(group, levels = c("Overall", "VT", "NVT")))
}


############################################################
#### 2) Data ###############################################
############################################################

dat <- load_analysis_data()$dat
print(describe_density(dat))


############################################################
#### 3) Fit and export each configuration ##################
############################################################

if (length(configs_to_run) > 0) mod <- carriage_model("density")

for (cfg_name in configs_to_run) {

  cfg <- configs[[cfg_name]]
  message("\n##### Configuration: ", cfg_name, " #####")

  fits <- do.call(fit_density, c(
    list(dat, group = "all", covariates = covariates_for(cfg$run_type),
         interaction = cfg$interaction, multiplicity = mult_source,
         m_rep = m_rep, priors = PRIORS, model = mod),
    SAMPLER
  ))
  save_fit(fits, file.path(fit_dir, paste0("Density_fit_", cfg_name, "_", STAMP, ".rds")))

  results <- purrr::imap(fits, \(f, g) {
    list(
      diagnostics = sampler_diagnostics(f, group_labels[[g]]),
      summaries = summarise_density_fit(f$fit, fit_prep(f), group_labels[[g]])
    )
  })

  summaries <- list(
    overall_density = collect(results, "overall_density"),
    density_by_m = collect(results, "density_by_m"),
    pr_density = collect(results, "pr_density"),
    shift_diagnostics = collect(results, "shift_diagnostics"),
    multiplicity_effects = collect(results, "multiplicity_effects"),
    parameters = collect(results, "parameters"),
    observed_density = collect(results, "observed_density"),
    sampler_diagnostics = purrr::map(results, \(r) r$diagnostics) |>
      purrr::list_rbind()
  )
  print(summaries$sampler_diagnostics)

  # The covariate part of the survey effect is a pure additive shift on
  # the log scale whenever the design has no survey-by-covariate
  # interaction. A non-negligible shift_maxdev means only PR_meanD and
  # PR_shiftD stay exact.
  if (max(summaries$shift_diagnostics$median, na.rm = TRUE) > 1e-6) {
    warning(
      "[", cfg_name, "] shift_maxdev is non-negligible: the design contains ",
      "survey-by-covariate variation, so PR_medianD is not an exact shift. ",
      "Report PR_meanD."
    )
  }

  # Median density ratios versus baseline, by group (Table S9 layout)
  pr_median_wide <- summaries$pr_density |>
    dplyr::filter(variable == "PR_medianD") |>
    dplyr::select(group, year, PR) |>
    dplyr::mutate(analysis = cfg$run_type) |>
    tidyr::pivot_wider(names_from = year, values_from = PR)

  saveRDS(
    summaries,
    file.path(out_dir, paste0("Density_summaries_", cfg_name, "_", STAMP, ".rds"))
  )

  openxlsx::write.xlsx(
    list(
      Overall_density = summaries$overall_density,
      Density_by_multiplicity = summaries$density_by_m,
      PR_density = summaries$pr_density,
      PR_median_wide = pr_median_wide,
      Multiplicity_effects = summaries$multiplicity_effects,
      Parameters = summaries$parameters,
      Shift_diagnostics = summaries$shift_diagnostics,
      Sampler_diagnostics = summaries$sampler_diagnostics,
      Observed_density = summaries$observed_density
    ),
    file = file.path(out_dir, paste0("Density_PR_", cfg_name, "_", STAMP, ".xlsx")),
    overwrite = TRUE
  )
}


############################################################
#### 4) Tables S6 and S9 ###################################
############################################################
# Built from the saved summaries, so this section also runs on its own
# after the configurations have been fitted.

read_summaries <- function(cfg_name) {
  readRDS(file.path(out_dir, paste0("Density_summaries_", cfg_name, "_", STAMP, ".rds")))
}
have <- function(cfg) file.exists(file.path(out_dir, paste0("Density_summaries_", cfg, "_", STAMP, ".rds")))

analysis_labels <- c(unadj = "Unadjusted", adj = "Adjusted")


# ---- Table S6: survey-specific multiplicity density ratios, all pneumococci ----
# Density ratios for every survey; ratios of density ratios for the
# post-PCV10 surveys relative to October 2016.

contrast_labels <- c(
  "Two versus single serotype" = "Two versus single",
  "Three versus single serotype" = "Three versus single",
  "Three versus two serotypes" = "Three versus two"
)

if (have("interaction_unadj") && have("interaction_adj")) {
  table_s6 <- purrr::map(
    c(unadj = "interaction_unadj", adj = "interaction_adj"),
    \(cfg_name) read_summaries(cfg_name)$multiplicity_effects
  ) |>
    purrr::list_rbind(names_to = "run_type") |>
    dplyr::filter(
      group == "Overall",
      !(estimand == "Ratio of density ratios" & as.character(year) == baseline_level)
    ) |>
    dplyr::mutate(
      class = factor(analysis_labels[run_type], levels = analysis_labels),
      group1 = factor(
        dplyr::if_else(estimand == "Density ratio", "Density ratios", estimand),
        levels = c("Density ratios", "Ratio of density ratios")
      ),
      contrast = factor(contrast_labels[contrast], levels = contrast_labels)
    ) |>
    dplyr::select(group1, class, contrast, year, estimate) |>
    dplyr::arrange(class, group1, contrast, year) |>
    tidyr::pivot_wider(names_from = year, values_from = estimate)

  openxlsx::write.xlsx(
    table_s6,
    file.path(out_dir, paste0("TableS6_Density ratios_", STAMP, ".xlsx")),
    overwrite = TRUE
  )
}


# ---- Table S9: median density ratios versus October 2016, by group ----

if (have("unadj") && have("adj")) {
  table_s9 <- purrr::map(
    c(unadj = "unadj", adj = "adj"),
    \(cfg_name) read_summaries(cfg_name)$pr_density
  ) |>
    purrr::list_rbind(names_to = "run_type") |>
    dplyr::filter(variable == "PR_medianD") |>
    dplyr::mutate(class = factor(analysis_labels[run_type], levels = analysis_labels)) |>
    dplyr::select(class, group, year, PR) |>
    dplyr::arrange(class, group, year) |>
    tidyr::pivot_wider(names_from = year, values_from = PR)

  openxlsx::write.xlsx(
    table_s9,
    file.path(out_dir, paste0("TableS9_Median density ratios_", STAMP, ".xlsx")),
    overwrite = TRUE
  )
}

message("\nDone. Written to ", out_dir, "/.")
