############################################################
# Sensitivity analysis: serotype 6A reclassified as vaccine-associated
# Approaches i-iv, unified framework (Fig S7)
#
# 6A is added to the vaccine-type group (cross-protection from the 6B
# component of PCV10): VT = PCV10 + 6A, NVT = the remaining non-PCV10
# serotypes. 6A keeps its own serotype category, so the data,
# likelihood and posterior are those of the primary analysis; only the
# VT/NVT sums change. They are therefore computed here, draw by draw,
# from the primary fits of 02_approaches_1_4 (no new fits):
#   P_VT  = sum of P_sero over PCV10 + 6A;   P_NVT = the other serotypes
#   share_VT, share_NVT likewise;  ratios versus Oct. 2016
#   15B/C + 15A + 23A as a share of the remaining NVT carriage
#
# Reads:  results/02_approaches_1_4/fits/M1to4_fit_M<a>_<run>_<STAMP>.rds
# Writes: results/06_sensitivity/6A/Sensi_6A_All_M<a>_<run>_re0_<STAMP>.xlsx
#         (layout read by 07_plots/Sensi_6A_Preva_PR_plot_<STAMP>.R)
#
# Minutes on the Mac. Run from the project directory, after the
# primary fits:
#   Rscript "06_sensitivity/6A/Sensi_6A_20261015.R"
############################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(posterior)
library(openxlsx)

source("00_setup/settings_20261015.R")
source("00_setup/helpers_20261015.R")

VT_6A <- c(VT_SEROS, "6A")
out_dir <- dir_create("results/06_sensitivity/6A")
fit_file <- function(a, run_type) {
  file.path("results/02_approaches_1_4/fits",
            paste0("M1to4_fit_M", a, "_", run_type, "_", STAMP, ".rds"))
}

sensi_6A <- function(fit, run_type) {
  lev  <- fit$data$survey_levels
  b    <- fit$data$reference
  cats <- fit$labels$categories
  stopifnot("6A" %in% cats, !"6A" %in% fit$data$vaccine_types)
  vt  <- cats %in% VT_6A
  nch <- fit$sampler$chains
  method <- method_of(fit)

  # sum of an rvar [survey, category] over a set of categories
  arr_sum <- function(rv, idx) {
    rvar(apply(draws_of(rv)[, , idx, drop = FALSE], c(1, 2), sum), nchains = nch)
  }
  P_sero <- stan_rvar(fit, "P_sero")
  share  <- stan_rvar(fit, "share")
  P <- list(overall = stan_rvar(fit, "P_overall"),
            VT = arr_sum(P_sero, vt), NVT = arr_sum(P_sero, !vt))
  S <- list(VT = arr_sum(share, vt), NVT = arr_sum(share, !vt))

  # draw-level check: VT + NVT = overall carriage
  stopifnot(max(abs(draws_of(P$VT + P$NVT - P$overall))) < 1e-6)

  by_class <- function(x, stem, scale = 1) {
    bind_rows(lapply(names(x), function(cl) {
      summarise_rvar(x[[cl]] * scale, paste0(stem, cl), lev) |> mutate(class = cl)
    }))
  }

  preva_class <- by_class(P, "P_") |>
    mutate(prevalence = fmt_ci(100 * median, 100 * q2.5, 100 * q97.5, 2), method = method)
  pr_class <- by_class(lapply(P, function(p) p / p[b]), "PR_") |>
    filter(year != BASELINE) |>
    mutate(PR = fmt_ci(median, q2.5, q97.5, 2), method = method)
  share_class <- by_class(S, "share_") |>
    mutate(Prop_among_carriers = fmt_ci(median, q2.5, q97.5, 3), method = method)

  target <- cats %in% TARGET_GROUP[[1]]
  stopifnot(sum(target) == length(TARGET_GROUP[[1]]))
  prop_3 <- summarise_rvar(100 * arr_sum(P_sero, target) / P$NVT,
                           "prop_15BC_15A_23A_among_NVT", lev) |>
    mutate(serotype_group = names(TARGET_GROUP),
           proportion = fmt_ci(median, q2.5, q97.5, 1), method = method)

  list(
    Run_metadata     = fit_metadata(fit, run_type,
                                    analysis = "Sensitivity: 6A reclassified as vaccine-associated",
                                    vt_definition = paste(VT_6A, collapse = ", ")),
    Preva_Class_long = preva_class,
    PR_Class_long    = pr_class,
    Share_Class_long = share_class,
    Prop_3_among_NVT = prop_3
  )
}

for (run_type in RUN_TYPES) {
  for (a in 1:4) {
    f <- fit_file(a, run_type)
    if (!file.exists(f)) { message("Missing: ", f); next }
    fit <- readRDS(f)
    xlsx_file <- file.path(out_dir, paste0("Sensi_6A_All_M", a, "_", run_type, "_re0_", STAMP, ".xlsx"))
    write.xlsx(sensi_6A(fit, run_type), xlsx_file, overwrite = TRUE)
    message("Written: ", xlsx_file)
  }
}
