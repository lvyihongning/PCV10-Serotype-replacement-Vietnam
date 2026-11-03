############################################################
# Helpers shared by the analysis scripts
#
# The figure scripts in 07_plots read result workbooks with fixed
# sheet and column names (year, class / serotype / category, median,
# q2.5, q97.5, method "Method <a>", ...). These helpers put the
# pneumoalloc results into that layout, so the figure code is
# unchanged, and summarise the few Stan quantities that the package
# does not tabulate itself (e.g. prevalence differences).
#
# Sourced after settings:
#   source("00_setup/settings_20261015.R")
#   source("00_setup/helpers_20261015.R")
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(posterior)
})

SUMMARY_COLS <- c("mean", "median", "sd", "q2.5", "q97.5",
                  "rhat", "ess_bulk", "ess_tail", "mcse_mean")

# "median (lower, upper)"
fmt_ci <- function(med, lo, hi, digits = 2, sep = ", ") {
  f <- function(v) formatC(v, digits = digits, format = "f")
  paste0(f(med), " (", f(lo), sep, f(hi), ")")
}

# Approach label of the earlier result files: "Method 1" ... "Method 4"
# (1 dominant serotype, 2 equal share, 3 relative abundance,
# 4 serotype-specific density)
method_of <- function(fit) paste0("Method ", fit$approach$approach)

# pneumoalloc rows -> earlier layout: survey -> year (text), the outcome
# in a column named `key` (class, serotype or category) and the Stan
# variable name in `variable`; `keep` names further columns to keep
legacy <- function(tab, key, variable, keep = NULL) {
  tab |>
    mutate(variable = variable, year = as.character(survey), "{key}" := outcome) |>
    select(variable, year, all_of(key), any_of(SUMMARY_COLS), any_of(keep))
}

# Posterior summaries of an rvar indexed [survey] or [survey, k]
summarise_rvar <- function(rv, variable, survey_levels, labels = NULL, key = "serotype") {
  s <- summarise_draws(
    rv, "mean", "median", "sd",
    q2.5  = ~ quantile2(.x, probs = 0.025),
    q97.5 = ~ quantile2(.x, probs = 0.975),
    "rhat", "ess_bulk", "ess_tail", "mcse_mean"
  )
  idx <- str_match(s$variable, "\\[(\\d+)(?:,(\\d+))?\\]")
  s$year <- survey_levels[as.integer(idx[, 2])]
  if (!is.null(labels)) s[[key]] <- labels[as.integer(idx[, 3])]
  s$variable <- variable
  s |> select(variable, year, any_of(key), any_of(SUMMARY_COLS))
}

# The same for a generated quantity of a pneumoalloc fit
summarise_stan <- function(fit, variable, labels = NULL, key = "serotype") {
  rv <- as_draws_rvars(fit$fit$draws(variable))[[variable]]
  summarise_rvar(rv, variable, fit$data$survey_levels, labels, key)
}

# Draws of a generated quantity as an rvar ([survey] or [survey, k])
stan_rvar <- function(fit, variable) {
  as_draws_rvars(fit$fit$draws(variable))[[variable]]
}

# Run metadata of one fit
fit_metadata <- function(fit, run_type, ...) {
  diag <- fit$fit$diagnostic_summary(quiet = TRUE)
  chk <- check_fit(fit)$summary
  extra <- list(...)
  tibble(
    item = c("analysis_stamp", "pneumoalloc_version", "model", "method",
             "run_type", "formula", "n", "chains", "iter_warmup", "iter_sampling",
             "seed", "num_divergent", "num_max_treedepth", "min_ebfmi",
             "max_rhat", "min_ess_bulk", "min_ess_tail",
             paste0("prior_", names(fit$priors)), names(extra)),
    value = as.character(c(
      STAMP, pneumoalloc_version(), fit$type,
      if (is.null(fit$approach)) NA else fit$approach$label,
      run_type,
      paste(c("~ 0 + survey", attr(stats::terms(fit$covariates %||% ~ 1), "term.labels")),
            collapse = " + "),
      length(fit$rows), fit$sampler$chains, fit$sampler$iter_warmup,
      fit$sampler$iter_sampling, fit$sampler$seed,
      sum(diag$num_divergent), sum(diag$num_max_treedepth), round(min(diag$ebfmi), 3),
      round(chk$value[chk$check == "max R-hat"], 4),
      round(chk$value[chk$check == "min bulk ESS"]),
      round(chk$value[chk$check == "min tail ESS"]),
      unlist(fit$priors), unlist(extra)
    ))
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
