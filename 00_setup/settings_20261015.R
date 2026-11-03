############################################################
# Shared settings for every analysis script
#
# Sourced at the top of each script, from the project directory:
#   source("00_setup/settings_20261015.R")
#
# One place for the stamp, paths, survey rounds, serotype sets,
# covariates, priors (Table S2) and sampler settings, so that all
# analyses use the same definitions. The models are those of the
# pneumoalloc package (github.com/lvyihongning/pneumoalloc); the
# installed version is recorded in every result workbook.
############################################################

suppressPackageStartupMessages(library(pneumoalloc))

# ==========================================================
# Stamp and paths
# ==========================================================
STAMP     <- "20261015"
DATA_FILE <- "data/carr_dat_clean_20260917.qs"

# Built by 00_setup/01_prepare_data_<STAMP>.R (base-R .rds, so the
# cluster needs no qs package)
ANALYSIS_DATA <- file.path("results/00_setup", paste0("analysis_data_", STAMP, ".rds"))

# ==========================================================
# Survey rounds and serotype sets
# ==========================================================
SURVEY_LEVELS <- c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018",
                   "Oct. 2019", "Oct. 2020", "July 2022")
BASELINE <- "Oct. 2016"

# PCV10 serotypes, in the order of the earlier result files
VT_SEROS  <- c("1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F")
stopifnot(setequal(VT_SEROS, pcv10_types()))

# Non-PCV10 serotypes with their own category: the 9 most common NVTs
# (union of the per-survey top 7). All other NVTs are pooled into
# "Other NVT", the reference category of the composition model.
NVT_SEROS  <- c("6A", "NT2", "15A", "19A", "15B/C", "11A", "35B", "23A", "34")
OTHER_NVT  <- "Other NVT"
CATEGORIES <- c(VT_SEROS, NVT_SEROS)

# Serotype group reported as a share of non-PCV10 carriage
TARGET_GROUP <- list(`15B/C + 15A + 23A` = c("15B/C", "15A", "23A"))

# ==========================================================
# Covariates
# ==========================================================
# Raw variables that define the complete-case analysis population
ADJ_VARS <- c("agemonth", "compchild", "trial_arm", "income",
              "famnum_cat", "moedu1", "ill1")

# Model variables (built in 00_setup/01_prepare_data) and the adjusted
# formula; survey round is always in the model (~ 0 + survey + ...)
COVARIATE_VARS <- c("agemonth_z", "compchild", "trial_arm", "income_z",
                    "famnum_cat", "moedu1", "ill1")
COVARIATES <- ~ agemonth_z + compchild + trial_arm + income_z + famnum_cat + moedu1 + ill1

RUN_TYPES <- c("unadj", "adj")
RUN_LABEL <- c(unadj = "unadjusted", adj = "adjusted")
covariates_for <- function(run_type) if (run_type == "adj") COVARIATES else NULL

# ==========================================================
# Priors (Table S2) and sampler
# ==========================================================
# Package defaults: carriage 2.5, composition 2, multiplicity 2,
# density survey 3 (centre log 1e6) and covariates 1, multiplicity
# density effects 1, sigma ~ exponential(0.5); carriage states 2.5,
# types of multiple carriage 2, presence-based 2.5.
PRIORS <- carriage_priors()

# 4 chains x (1000 warm-up + 2000 sampling) = 8000 draws per fit
SAMPLER <- list(chains = 4, iter_warmup = 1000, iter_sampling = 2000,
                seed = 123, max_treedepth = 12, refresh = 500)

APPROACHES <- c("dominant", "equal", "abundance", "density")   # approaches 1-4

# ==========================================================
# Small shared helpers
# ==========================================================
load_analysis_data <- function() {
  if (!file.exists(ANALYSIS_DATA)) {
    stop("Run 00_setup/01_prepare_data_", STAMP, ".R first: ", ANALYSIS_DATA, " not found.")
  }
  readRDS(ANALYSIS_DATA)
}

dir_create <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Version of pneumoalloc used, for the run metadata
pneumoalloc_version <- function() {
  d <- utils::packageDescription("pneumoalloc")
  paste0(d$Version, if (!is.null(d$RemoteSha)) paste0(" (", substr(d$RemoteSha, 1, 7), ")"))
}
