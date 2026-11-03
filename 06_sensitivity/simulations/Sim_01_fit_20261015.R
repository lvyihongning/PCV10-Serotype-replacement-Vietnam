############################################################
# 06_sensitivity/simulations/Sim_01_fit_20261015.R
#
# High-carriage simulations (Figs S9, S10) - step 1: fit.
#
# Each simulated population (Sim_00_simulate_data) is analysed exactly
# as the observed data, with pneumoalloc:
#   approaches 1-4  fit_carriage()  (unified framework; same allocation
#                   rules, serotype categories, priors and design)
#   presence-based  fit_presence()
# each unadjusted and adjusted. Baseline Oct. 2016, comparison Oct. 2019.
#
# 20 fits: sim (2) x run type (2) x approach (1-4, presence-based).
# Each fit's summaries are saved as it finishes, then combined.
#
# Mode (first command-line argument; default "all"):
#   all       fit all 20 runs in turn, then combine (local use)
#   compile   compile the two Stan models only (cluster, before the array)
#   <i>       fit run i (1-20, row of `runs`) and save it (cluster array;
#             cluster/run_Sim_fits.sbatch)
#   combine   combine the 20 saved fits and write the outputs
#
# Output (results/06_sensitivity/simulations/):
#   fits/Sim_fit_<sim>_<run_type>_a<approach>_<STAMP>.rds   one per fit
#   Sim_results_<STAMP>.xlsx     long tables for both simulations
#   Sim_results_<STAMP>.rds      the same tables, for Sim_02_figures
# Next: Sim_02_figures_20261015.R
#
# Populations of 2,000 children: a few minutes per fit on the Mac.
############################################################

rm(list = ls())
# Core tidyverse packages loaded one by one: the cluster has no tidyverse
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)
library(tibble)
library(openxlsx)

source("00_setup/settings_20261015.R")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "all" else args[1]


############################################################
#### 0) Settings ###########################################
############################################################

date_tag <- STAMP
out_dir <- "results/06_sensitivity/simulations"
data_file <- file.path(out_dir, paste0("Sim_data_", date_tag, ".rds"))
fit_dir <- file.path(out_dir, "fits")

sims_to_run <- c("sim1", "sim2")
approaches <- 1:4
run_types <- RUN_TYPES
fit_presence_based <- TRUE
reuse_saved_fits <- FALSE   # mode "all": skip runs whose fit file exists for the same data

runs <- expand_grid(
  sim = sims_to_run,
  run_type = run_types,
  approach = c(approaches, if (fit_presence_based) 5L)
)
fit_file <- function(i) {
  file.path(fit_dir, sprintf("Sim_fit_%s_%s_a%d_%s.rds",
                             runs$sim[i], runs$run_type[i], runs$approach[i], date_tag))
}

approach_info <- tibble(
  approach = 1:5,
  name = c(APPROACHES, "presence"),
  method = c("Dominant serotype", "Equal share", "Relative abundance",
             "Serotype-specific density", "Presence-based")
)

baseline_level <- BASELINE


############################################################
#### 1) Functions: analysis data ###########################
############################################################

# Covariate coding as in 00_setup/01_prepare_data, standardised within
# the simulated population
prepare_covariates <- function(sim) {
  sim |>
    mutate(
      survey_yearmon = factor(simulated_period, levels = c("Oct. 2016", "Oct. 2019")),
      famnum_cat = fct_relevel(factor(famnum_cat), "<=3"),
      moedu1 = fct_relevel(factor(moedu1), "1"),
      trial_arm = fct_drop(factor(trial_arm)) |> fct_relevel("0p+1"),
      ill1 = as.integer(ill1 == 1),   # raw coding 1 = yes, 2 = no -> 1 = ill, 0 = not ill
      agemonth_z = as.numeric(scale(agemonth)),
      income_z = as.numeric(scale(log(pmin(income, 100))))
    )
}

# pneumoalloc analysis data of one simulated population
build_analysis_data <- function(sim) {
  prepare_carriage_data(
    prepare_covariates(sim),
    survey = "survey_yearmon",
    serotypes = c("sero1", "sero2", "sero3"),
    abundance = c("ra_sero1", "ra_sero2", "ra_sero3"),
    density = "dens_tot",
    vaccine_types = VT_SEROS,
    categories = CATEGORIES,
    other_label = OTHER_NVT,
    covariates = COVARIATE_VARS,
    reference = baseline_level,
    dominant = "abundance"
  )
}


############################################################
#### 2) Functions: fit and summarise #######################
############################################################

fmt_ci <- function(med, lo, hi, digits = 2) {
  paste0(formatC(med, digits = digits, format = "f"), " (",
         formatC(lo, digits = digits, format = "f"), "-",
         formatC(hi, digits = digits, format = "f"), ")")
}

summary_cols <- c("mean", "median", "sd", "q2.5", "q97.5", "rhat", "ess_bulk", "ess_tail")
class_of <- c(Overall = "overall", VT = "VT", NVT = "NVT", `Multiple carriage` = "multiple")

# Common output shape: one row per year x class (or serotype)
summarise_fit <- function(fit) {
  mp <- marginal_prevalence(fit) |> filter(measure == "prevalence")
  pr <- prevalence_ratios(fit) |> filter(measure == "prevalence")
  tidy <- function(tab, key) {
    tab |>
      mutate(year = as.character(survey), "{key}" := outcome) |>
      select(year, all_of(key), all_of(summary_cols))
  }
  classes <- names(class_of)
  list(
    prev_class = mp |> filter(outcome %in% classes) |> tidy("class") |>
      mutate(class = unname(class_of[class])),
    pr_class = pr |> filter(outcome %in% classes) |> tidy("class") |>
      mutate(class = unname(class_of[class])),
    prev_sero = mp |> filter(group == "serotype") |> tidy("serotype"),
    pr_sero = pr |> filter(group == "serotype") |> tidy("serotype")
  )
}

# Crude allocation-weighted prevalence (density-weighted for approach 4)
# or crude presence; the unadjusted model should reproduce it up to
# prior shrinkage
crude_class <- function(dat, name) {
  cr <- describe_prevalence(dat, name) |>
    filter(group %in% c("overall", "VT", "NVT")) |>
    transmute(year = as.character(survey), class = group, crude = value)
  if (name == "presence") {
    cr <- bind_rows(cr, tibble(
      year = dat$survey_levels,
      class = "multiple",
      crude = as.vector(tapply(dat$n_sero >= 2, dat$survey, mean))
    ))
  }
  cr
}

sampler_diagnostics <- function(fit) {
  d <- fit$fit$diagnostic_summary(quiet = TRUE)
  tibble(divergences = sum(d$num_divergent), max_treedepth = sum(d$num_max_treedepth),
         min_ebfmi = min(d$ebfmi))
}

fit_run <- function(dat, approach, run_type) {
  name <- approach_info$name[approach]
  args <- c(list(dat, covariates = covariates_for(run_type), priors = PRIORS), SAMPLER)
  fit <- if (name == "presence") {
    do.call(fit_presence, args)
  } else {
    do.call(fit_carriage, c(args, list(approach = name)))
  }
  out <- summarise_fit(fit)
  out$crude_class <- crude_class(dat, name)
  out$diagnostics <- sampler_diagnostics(fit)
  out
}


############################################################
#### 3) Fit ################################################
############################################################

stopifnot(mode %in% c("all", "compile", "combine") || grepl("^[0-9]+$", mode))

if (mode == "compile") {
  carriage_model("unified")
  carriage_model("presence")
  message("Models compiled.")
  quit(save = "no")
}

sim_input <- readRDS(data_file)
# Fingerprint of each simulated population, so that changing one
# simulation leaves the fits of the other valid. Built from the exact
# stored values (numbers as IEEE doubles, text as UTF-8), so that it is
# the same on the Mac and the cluster.
data_fingerprint <- function(df) {
  tf <- tempfile()
  on.exit(unlink(tf))
  con <- file(tf, "wb")
  for (nm in names(df)) {
    x <- df[[nm]]
    writeBin(enc2utf8(nm), con)
    if (is.numeric(x) || is.logical(x)) {
      writeBin(as.double(x), con, endian = "little")
    } else {
      writeBin(enc2utf8(ifelse(is.na(x), "\001NA", as.character(x))), con)
    }
  }
  close(con)
  unname(tools::md5sum(tf))
}
sim_md5 <- map_chr(set_names(sims_to_run), \(s) data_fingerprint(sim_input$sim_data[[s]]))
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data <- map(set_names(sims_to_run), \(s) build_analysis_data(sim_input$sim_data[[s]]))

# Fit run i and save its summaries, tagged with the simulated data used
fit_one <- function(i) {
  sim <- runs$sim[i]
  run_type <- runs$run_type[i]
  approach <- runs$approach[i]
  message("\n================ ", i, ": ", sim, " | ", run_type, " | ",
          approach_info$method[approach], " ================")
  out <- fit_run(analysis_data[[sim]], approach, run_type)
  out$data_md5 <- sim_md5[[sim]]
  out$pneumoalloc_version <- pneumoalloc_version()
  saveRDS(out, fit_file(i))
  message("Written: ", fit_file(i))
}

fit_is_current <- function(i) {
  file.exists(fit_file(i)) &&
    identical(readRDS(fit_file(i))$data_md5, sim_md5[[runs$sim[i]]])
}

if (mode == "all") {
  todo <- seq_len(nrow(runs))
  if (reuse_saved_fits) todo <- discard(todo, fit_is_current)
  walk(todo, fit_one)
} else if (mode != "combine") {
  i <- as.integer(mode)
  stopifnot(i >= 1, i <= nrow(runs))
  fit_one(i)
  quit(save = "no")
}

# Every fit must exist and come from the current simulated data
stale <- discard(seq_len(nrow(runs)), fit_is_current)
if (length(stale) > 0) {
  stop("Missing fits, or fits of other simulated data: ",
       paste(basename(fit_file(stale)), collapse = ", "))
}
fits <- map(seq_len(nrow(runs)), \(i) readRDS(fit_file(i)))


############################################################
#### 4) Combine and export #################################
############################################################

collect <- function(field) {
  map2(fits, seq_len(nrow(runs)), \(f, i) {
    f[[field]] |>
      mutate(simulation = runs$sim[i],
             scenario = sim_input$scenarios[[runs$sim[i]]]$label,
             run_type = runs$run_type[i],
             approach = runs$approach[i],
             method = approach_info$method[runs$approach[i]],
             .before = 1)
  }) |>
    list_rbind()
}

results <- list(
  prev_class = collect("prev_class") |>
    mutate(prevalence = fmt_ci(median * 100, q2.5 * 100, q97.5 * 100, 1)),
  pr_class = collect("pr_class") |>
    filter(year != baseline_level) |>
    mutate(PR = fmt_ci(median, q2.5, q97.5)),
  prev_sero = collect("prev_sero") |>
    mutate(sero_group = if_else(serotype %in% VT_SEROS, "VT", "NVT"),
           prevalence = fmt_ci(median * 100, q2.5 * 100, q97.5 * 100, 1)),
  pr_sero = collect("pr_sero") |>
    filter(year != baseline_level) |>
    mutate(sero_group = if_else(serotype %in% VT_SEROS, "VT", "NVT"),
           PR = fmt_ci(median, q2.5, q97.5)),
  crude_check = collect("prev_class") |>
    select(simulation, run_type, method, year, class, median) |>
    left_join(collect("crude_class") |> select(simulation, run_type, method, year, class, crude),
              by = c("simulation", "run_type", "method", "year", "class")) |>
    mutate(diff_pct_points = (median - crude) * 100),
  diagnostics = collect("diagnostics")
)

message("\nSampler diagnostics:")
print(results$diagnostics, n = Inf)
message("Max |model - crude| (percentage points), unadjusted: ",
        round(max(abs(filter(results$crude_check, run_type == "unadj")$diff_pct_points),
                  na.rm = TRUE), 3))

# Wide summary of the key estimates: unadjusted prevalence, adjusted PR
key_table <- bind_rows(
  results$prev_class |> filter(run_type == "unadj", class %in% c("VT", "NVT")) |>
    transmute(simulation, class, method, estimate = paste0("Prevalence ", year), value = prevalence),
  results$pr_class |> filter(run_type == "adj", class %in% c("VT", "NVT")) |>
    transmute(simulation, class, method, estimate = "aPR Oct. 2019", value = PR)
) |>
  pivot_wider(names_from = estimate, values_from = value) |>
  arrange(simulation, desc(class), match(method, approach_info$method))

saveRDS(c(results, list(key_table = key_table,
                        sero_levels = map(analysis_data, "sero_levels"))),
        file.path(out_dir, paste0("Sim_results_", date_tag, ".rds")))

write.xlsx(
  list(Key_estimates = key_table, Preva_Class_long = results$prev_class,
       PR_Class_long = results$pr_class, Preva_Sero_long = results$prev_sero,
       PR_Sero_long = results$pr_sero, Crude_check = results$crude_check,
       Sampler_diagnostics = results$diagnostics),
  file.path(out_dir, paste0("Sim_results_", date_tag, ".xlsx")),
  overwrite = TRUE
)

message("\nDone. Written to ", out_dir, "/.")
