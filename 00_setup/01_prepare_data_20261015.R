############################################################
# Step 0: analysis data for every model
#
# Builds, once, the pneumoalloc analysis data used by all analyses:
#   - analysis population: complete cases for the adjusted covariates
#     (and commune, number of serotypes), so unadjusted and adjusted
#     models of every analysis use the same children
#   - covariate coding (as in the earlier scripts)
#   - carriers, density, number of detected serotypes, dominant
#     serotype, type of multiple carriage
#   - serotype categories: PCV10 serotypes, the 9 retained non-PCV10
#     serotypes and "Other NVT" (reference category, last)
# and checks them against the raw fields of the data set.
#
# Input:  data/carr_dat_clean_20260917.qs   (needs qs; run on the Mac)
# Output: results/00_setup/analysis_data_<STAMP>.rds
#         results/00_setup/analysis_data_checks_<STAMP>.xlsx
# Next:   any analysis script; copy the .rds to the cluster for the
#         cluster runs (see README).
#
# Run from the project directory.
############################################################

rm(list = ls())

library(qs)
library(dplyr)
library(tidyr)
library(forcats)
library(openxlsx)

source("00_setup/settings_20261015.R")

out_dir <- dir_create("results/00_setup")

# Columns kept in the analysis data (models, crude checks and figures);
# columns absent from the data set are skipped
KEEP_COLS <- c(
  "id", "survey_yearmon", "year", "commune_ID",
  "agemonth", "compchild", "trial_arm", "income", "famnum_cat", "moedu1", "ill1",
  "pcv", "tot", "pos", "pneu_carr", "pclytqpcr",
  "num_sero", "sero1", "sero2", "sero3", "ra_sero1", "ra_sero2", "ra_sero3",
  "dens_sero1", "dens_sero2", "dens_sero3", "dens_tot", "log_dens_tot",
  "domin_sero", "domin_subdomin_sero_1"
)

# ==========================================================
# 1) LOAD AND CODE COVARIATES
# ==========================================================
raw <- qread(DATA_FILE) |> ungroup()
View(raw)
table(raw$pneu_carr)

n_income_capped <- sum(raw$income > 100, na.rm = TRUE)
message("Income values capped at 100: ", n_income_capped)

# Survey labels: "Oct. 2017" is the November 2017 survey (as relabelled
# in the earlier figure scripts)
survey_raw <- as.character(raw$survey_yearmon)
if ("Oct. 2017" %in% survey_raw) {
  stopifnot(!"Nov. 2017" %in% survey_raw)
  survey_raw[survey_raw == "Oct. 2017"] <- "Nov. 2017"
  message("Survey label 'Oct. 2017' relabelled 'Nov. 2017'.")
}
unknown <- setdiff(unique(survey_raw), SURVEY_LEVELS)
if (length(unknown) > 0) stop("Unexpected survey labels: ", paste(unknown, collapse = ", "))

# Standardisation on the full data set, before the complete-case
# restriction, as in the earlier scripts
df <- raw |>
  select(any_of(KEEP_COLS)) |>
  mutate(
    survey_yearmon = factor(survey_raw, levels = intersect(SURVEY_LEVELS, survey_raw)),
    famnum_cat     = fct_relevel(factor(famnum_cat), "<=3"),
    moedu1         = fct_relevel(factor(moedu1), "1"),
    trial_arm      = fct_drop(trial_arm) |> fct_relevel("0p+1"),  # drop empty "Unvaccinated" level
    ill1           = as.integer(ill1 == 1),   # raw coding 1 = yes, 2 = no -> 1 = ill, 0 = not ill
    agemonth_z     = as.numeric(scale(agemonth)),
    income_z       = as.numeric(scale(log(pmin(income, 100)))),
    across(c(sero1, sero2, sero3), ~ na_if(trimws(as.character(.x)), ""))
  )
View(df)
# ==========================================================
# 2) ANALYSIS POPULATION
# ==========================================================
cc <- complete.cases(df[, c(ADJ_VARS, "survey_yearmon", "commune_ID", "num_sero")])
n_dropped <- sum(!cc)
message("Children excluded for missing covariate, commune or serotype data: ",
        n_dropped, " of ", nrow(df), " (", round(100 * n_dropped / nrow(df), 2), "%)")
df <- df[cc, , drop = FALSE]

# ==========================================================
# 3) PNEUMOALLOC ANALYSIS DATA
# ==========================================================
dat <- prepare_carriage_data(
  df,
  survey        = "survey_yearmon",
  serotypes     = c("sero1", "sero2", "sero3"),
  abundance     = c("ra_sero1", "ra_sero2", "ra_sero3"),
  density       = "dens_tot",
  vaccine_types = VT_SEROS,
  categories    = CATEGORIES,
  other_label   = OTHER_NVT,
  covariates    = COVARIATE_VARS,
  reference     = BASELINE,
  dominant      = "abundance"     # highest relative abundance; ties by slot order
)
print(dat)
stopifnot(dat$n_dropped == 0)      # complete cases were taken above

# ==========================================================
# 4) CHECKS AGAINST THE RAW FIELDS
# ==========================================================
checks <- list()

# Number of detected serotypes = num_sero
stopifnot(all(dat$n_sero == dat$data$num_sero))

# Serotype categories: PCV10, retained NVTs, "Other NVT" last
sero_long <- dat$sero_long
expected_levels <- c(VT_SEROS[VT_SEROS %in% sero_long$category],
                     NVT_SEROS[NVT_SEROS %in% sero_long$category], OTHER_NVT)
stopifnot(identical(dat$sero_levels, expected_levels))
checks$Categories <- tibble(serotype = dat$sero_levels, vaccine_type = dat$vt_mask == 1)

# Dominant serotype (highest relative abundance) = the raw sero1 field
dom <- sero_long[sero_long$dominant, ]
dom <- dom[order(dom$k), ]
sero1_cat <- if_else(dat$data$sero1[dat$pos_idx] %in% CATEGORIES,
                     dat$data$sero1[dat$pos_idx], OTHER_NVT)
n_dom_mismatch <- sum(dom$category != sero1_cat)
message("Carriers whose dominant category differs from sero1: ", n_dom_mismatch)
stopifnot(n_dom_mismatch == 0)

# Type of multiple carriage versus the raw dominant/subdominant pattern
if ("domin_subdomin_sero_1" %in% names(dat$data)) {
  CLASS_MAP <- c(
    "VT+VT"  = "VT-only",  "VT+VT+VT"    = "VT-only",
    "NVT+NVT" = "NVT-only", "NVT+NVT+NVT" = "NVT-only",
    "VT+NVT" = "Mixed, VT-dominant", "VT+VT+NVT" = "Mixed, VT-dominant",
    "VT+NVT+NVT" = "Mixed, VT-dominant", "VT+NVT+VT" = "Mixed, VT-dominant",
    "NVT+VT" = "Mixed, NVT-dominant", "NVT+NVT+VT" = "Mixed, NVT-dominant",
    "NVT+VT+VT" = "Mixed, NVT-dominant", "NVT+VT+NVT" = "Mixed, NVT-dominant"
  )
  multi <- dat$m_obs >= 2
  type_raw <- unname(CLASS_MAP[dat$data$domin_subdomin_sero_1[dat$pos_idx][multi]])
  checks$Multiple_types <- as.data.frame(table(
    package = dat$multiple_type[multi], raw_pattern = type_raw, useNA = "ifany"
  ))
  n_type_mismatch <- sum(dat$multiple_type[multi] != type_raw, na.rm = TRUE)
  message("Multiple carriers whose type differs from domin_subdomin_sero_1: ",
          n_type_mismatch, " of ", sum(multi))
}

# Carrier definition (>= 1 serotype) versus the lytA qPCR result
if ("pclytqpcr" %in% names(dat$data)) {
  checks$Carrier_vs_lytA <- as.data.frame(table(
    survey = dat$survey, lytA = dat$data$pclytqpcr, carrier = dat$carrier, useNA = "ifany"
  ))
}

# Crude counts and allocation totals per survey, category and approach
checks$States <- describe_carriage_states(dat)
checks$Allocation_totals <- bind_rows(lapply(APPROACHES, function(a) {
  w <- allocate_weights(dat, a)
  as_tibble(rowsum(w, dat$survey[dat$pos_idx]), rownames = "survey") |>
    mutate(approach = a, .before = 1)
}))
checks$Population <- tibble(
  item  = c("children", "carriers", "multiple carriers", "excluded (missing data)",
            "income capped at 100", "max detected serotypes"),
  value = c(dat$N, dat$N_pos, sum(dat$m_obs >= 2), n_dropped, n_income_capped, dat$M_observed)
)

# ==========================================================
# 5) SAVE
# ==========================================================
saveRDS(
  list(stamp = STAMP, data_file = DATA_FILE, dat = dat,
       n_dropped = n_dropped, n_income_capped = n_income_capped,
       pneumoalloc_version = pneumoalloc_version()),
  ANALYSIS_DATA
)
write.xlsx(checks, file.path(out_dir, paste0("analysis_data_checks_", STAMP, ".xlsx")),
           overwrite = TRUE)
message("Written: ", ANALYSIS_DATA)
