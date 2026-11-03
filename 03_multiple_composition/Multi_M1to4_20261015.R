############################################################
# Approaches i-iv among children with multiple serotype carriage
# Unified framework — fit, summarise and compare approaches
# (Table S4; Fig S2)
#
# pneumoalloc::fit_carriage(population = "multiple"): serotype
# composition of multiple carriage (fractional multinomial logit),
# density (lognormal) and three versus two serotypes (logistic),
# among the children of the analysis population with >= 2 detected
# serotypes. Serotype categories: those of the primary analysis that
# are detected among multiple carriers ("Other NVT" last).
#
# Estimand: share_s(t), standardised over the multiple carriers;
# VT/NVT shares; ratios and differences versus Oct. 2016.
#
# Writes (results/03_multiple_composition/):
#   fits/Multi_M1to4_fit_M<a>_<run>_<STAMP>.rds
#   Multi_M1to4_results_re0_<STAMP>.xlsx   all approaches x run types,
#                                          layout read by 07_plots
#
# Small population: a few minutes per fit on the Mac.
# Run from the project directory:
#   Rscript "03_multiple_composition/Multi_M1to4_20261015.R"                  all 8 fits
#   Rscript "03_multiple_composition/Multi_M1to4_20261015.R" <approach> <run> one fit
#   Rscript "03_multiple_composition/Multi_M1to4_20261015.R" combine          tables only
############################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(openxlsx)

source("00_setup/settings_20261015.R")
source("00_setup/helpers_20261015.R")

approaches <- 1:4
run_types  <- RUN_TYPES
do_fit     <- TRUE
cl_args <- commandArgs(trailingOnly = TRUE)
if (length(cl_args) == 1 && cl_args == "combine") {
  do_fit <- FALSE
} else if (length(cl_args) == 2) {
  approaches <- as.integer(cl_args[1])
  run_types  <- cl_args[2]
  stopifnot(approaches %in% 1:4, run_types %in% RUN_TYPES)
}

out_dir <- dir_create("results/03_multiple_composition")
fit_dir <- dir_create(out_dir, "fits")
fit_file <- function(a, run_type) {
  file.path(fit_dir, paste0("Multi_M1to4_fit_M", a, "_", run_type, "_", STAMP, ".rds"))
}

dat <- load_analysis_data()$dat
vt  <- dat$vaccine_types

# ==========================================================
# 1) FIT
# ==========================================================
if (do_fit) {
  mod <- carriage_model("unified_multiple")
  for (run_type in run_types) {
    for (a in approaches) {
      message("\n================ M", a, " (", APPROACHES[a], "), ", run_type, " ================")
      fit <- do.call(fit_carriage, c(
        list(dat, approach = APPROACHES[a], population = "multiple",
             covariates = covariates_for(run_type), priors = PRIORS, model = mod),
        SAMPLER
      ))
      print(check_fit(fit)$summary)
      save_fit(fit, fit_file(a, run_type))
    }
  }
}

# ==========================================================
# 2) COMBINE ALL FITS PRESENT
# ==========================================================
runs <- expand_grid(run_type = RUN_TYPES, a = 1:4) |>
  mutate(file = fit_file(a, run_type)) |>
  filter(file.exists(file))
stopifnot(nrow(runs) > 0)
message("\nCombining ", nrow(runs), " fits")
fits <- lapply(runs$file, readRDS)

# Labels of the earlier workbook: method "Method <a>", group
# unadjusted/adjusted, approach name and run type
label_rows <- function(tab, fit, run_type) {
  # look up the label outside mutate(): inside it, `run_type` would refer to
  # the column just created rather than the argument
  run_label <- RUN_LABEL[[run_type]]
  tab |>
    mutate(approach = fit$approach$approach, method = method_of(fit),
           approach_name = fit$approach$label, run_type = run_type,
           group = run_label, .before = 1)
}


per_fit <- lapply(seq_along(fits), function(i) {
  fit <- fits[[i]]
  rt  <- runs$run_type[i]
  mp  <- marginal_prevalence(fit)
  pr  <- prevalence_ratios(fit)
  chk <- check_fit(fit)
  list(
    share_class = mp |> filter(group %in% c("VT", "NVT")) |>
      legacy("class", "share") |> label_rows(fit, rt),
    share_sero = mp |> filter(group == "serotype") |>
      legacy("serotype", "share") |> label_rows(fit, rt),
    pr_class = pr |> filter(group %in% c("VT", "NVT")) |>
      legacy("class", "PR") |> label_rows(fit, rt),
    # PRs of categories near zero at baseline are flagged (0.1% threshold)
    pr_sero = pr |> filter(group == "serotype") |>
      mutate(PR_reportable = reportable) |>
      legacy("serotype", "PR", keep = "PR_reportable") |>
      label_rows(fit, rt),
    pd_sero = summarise_stan(fit, "PD_sero", fit$labels$categories) |>
      filter(year != BASELINE) |> label_rows(fit, rt),
    dens_mult = bind_rows(
      marginal_density(fit) |> legacy("outcome", "meanD"),
      density_ratios(fit) |> legacy("outcome", "PR_meanD"),
      mp |> filter(group == "multiplicity") |> legacy("outcome", "P_m3"),
      density_ratios(fit, multiplicity = TRUE) |> filter(is.na(survey)) |>
        legacy("outcome", "dens_fold_3vs2")
    ) |> label_rows(fit, rt),
    crude = chk$crude |> mutate(year = as.character(survey), .keep = "unused") |>
      label_rows(fit, rt),
    meta = fit_metadata(fit, rt) |>
      pivot_wider(names_from = item, values_from = value) |>
      mutate(approach = fit$approach$approach, .before = 1)
  )
})
grab <- function(nm) bind_rows(lapply(per_fit, `[[`, nm))

share_class <- grab("share_class") |>
  mutate(estimate = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 1))
share_sero <- grab("share_sero") |>
  mutate(sero_group = if_else(serotype %in% vt, "VT", "NVT"),
         estimate = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 1))
pr_class <- grab("pr_class") |> mutate(estimate = fmt_ci(median, q2.5, q97.5, 2))
pr_sero <- grab("pr_sero") |>
  mutate(sero_group = if_else(serotype %in% vt, "VT", "NVT"),
         estimate = fmt_ci(median, q2.5, q97.5, 2))
pd_sero <- grab("pd_sero") |>
  mutate(estimate = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 1))
dens_mult <- grab("dens_mult") |>
  mutate(estimate = case_when(
    variable == "meanD" ~ fmt_ci(log10(median), log10(q2.5), log10(q97.5), 2),
    variable == "P_m3"  ~ fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 1),
    TRUE                ~ fmt_ci(median, q2.5, q97.5, 2)
  ))

# Wide tables: one row per method x run type x class/serotype
widen <- function(df, ...) df |>
  select(group, method, ..., year, estimate) |>
  pivot_wider(names_from = year, values_from = estimate)

# Density and multiplicity must agree across approaches (Monte Carlo error)
shared_check <- bind_rows(lapply(RUN_TYPES, function(rt) {
  idx <- which(runs$run_type == rt)
  if (length(idx) < 2) return(NULL)
  f <- fits[idx]
  ref <- if (4 %in% runs$a[idx]) "density" else APPROACHES[runs$a[idx][1]]
  compare_approaches(f, reference_approach = ref)$shared |> mutate(run_type = rt, .before = 1)
}))
if (nrow(shared_check) > 0) {
  message("Shared quantities: max |z| = ", round(max(abs(shared_check$z), na.rm = TRUE), 2),
          "; flagged = ", sum(shared_check$flag, na.rm = TRUE))
}

allocation_totals <- bind_rows(lapply(1:4, function(a) {
  w <- allocate_weights(dat, APPROACHES[a], population = "multiple")
  k <- dat$m_obs >= 2
  as_tibble(rowsum(w, dat$survey[dat$pos_idx][k]), rownames = "year") |>
    mutate(method = paste0("Method ", a), .before = 1)
}))

xlsx_file <- file.path(out_dir, paste0("Multi_M1to4_results_re0_", STAMP, ".xlsx"))
write.xlsx(
  list(Run_metadata         = grab("meta"),
       Share_class          = share_class,
       Share_sero           = share_sero,
       PR_class             = pr_class,
       PR_sero              = pr_sero,
       PD_sero              = pd_sero,
       Share_class_table    = widen(share_class, class),
       Share_sero_table     = widen(share_sero, sero_group, serotype),
       PR_class_table       = widen(pr_class, class),
       PR_sero_table        = widen(pr_sero, sero_group, serotype),
       Density_multiplicity = dens_mult,
       Crude_check          = grab("crude"),
       Shared_check         = shared_check,
       Allocation_totals    = allocation_totals),
  file = xlsx_file, overwrite = TRUE
)
message("Written: ", xlsx_file)
