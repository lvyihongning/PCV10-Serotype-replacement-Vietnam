############################################################
# 06_sensitivity/arm/Arm_comparison_20261015.R
#
# Sensitivity analysis: can the four trial arms be pooled? (Table S10)
#
# Within each trial arm (0p+1, 1p+1, 2p+1, 3p+0) of the analysis
# population:
#   1. unadjusted multinomial model for the number of detected
#      serotypes (pneumoalloc::fit_carriage_states()), giving overall
#      carriage prevalence, multiple serotype carriage prevalence and
#      prevalence ratios versus October 2016
#   2. crude prevalences, as a check on the model
#   3. observed density among carriers (median and IQR, log10 GE/mL)
#
# Output (results/06_sensitivity/arm/):
#   Sensi_Arm_comparison_<STAMP>.xlsx       all tables, including Table_S10
#   figures/NotInManuscript_Prevalence and density by arm_<STAMP>.{tiff,pptx}
#
# 4 small fits; minutes on the Mac. Run from the project directory.
############################################################

rm(list = ls())
library(tidyverse)
library(openxlsx)
library(patchwork)
library(eoffice)

source("00_setup/settings_20261015.R")


############################################################
#### 0) Settings ###########################################
############################################################

arms <- c("0p+1", "1p+1", "2p+1", "3p+0")
baseline_level <- BASELINE
out_dir <- dir_create("results/06_sensitivity/arm")
dir_create(out_dir, "figures")

carriage_label <- "Carriage (>=1 serotype)"
multiple_label <- "Multiple carriage (>=2 serotypes)"
# pneumoalloc outcome labels -> labels of the earlier tables
outcome_label <- c("No carriage" = "No carriage", "1 serotype" = "One serotype",
                   "2 serotypes" = "Two serotypes", "3 serotypes" = "Three serotypes",
                   "3+ serotypes" = "Three serotypes",
                   "Carriage" = carriage_label, "Multiple carriage" = multiple_label)


############################################################
#### 1) Functions ##########################################
############################################################

fmt_ci <- function(med, lo, hi, digits = 1) {
  paste0(
    formatC(med, digits = digits, format = "f"), " (",
    formatC(lo, digits = digits, format = "f"), "-",
    formatC(hi, digits = digits, format = "f"), ")"
  )
}

# Fit the multinomial model within one arm and summarise it
fit_arm <- function(arm) {

  message("\n=== Trial arm ", arm, " ===")

  # pneumoalloc analysis data of the children of this arm
  dat_arm <- prepare_carriage_data(
    dplyr::filter(base_dat, trial_arm == .env$arm),   # the data also has a column `arm`
    survey = "survey_yearmon", serotypes = c("sero1", "sero2", "sero3"),
    density = "dens_tot", vaccine_types = VT_SEROS, categories = CATEGORIES,
    other_label = OTHER_NVT, reference = baseline_level
  )

  fit <- do.call(fit_carriage_states, c(list(dat_arm, priors = PRIORS), SAMPLER))
  chk <- check_fit(fit)

  diagnostics <- chk$summary |>
    dplyr::select(check, value) |>
    tidyr::pivot_wider(names_from = check, values_from = value) |>
    dplyr::mutate(trial_arm = arm, .before = 1)

  prevalence <- marginal_prevalence(fit) |>
    dplyr::transmute(
      year = survey, category = unname(outcome_label[outcome]),
      mean, median, sd, q2.5, q97.5, rhat, ess_bulk, ess_tail,
      trial_arm = arm,
      prevalence = fmt_ci(median * 100, q2.5 * 100, q97.5 * 100, digits = 1)
    )

  pr <- prevalence_ratios(fit) |>
    dplyr::transmute(
      year = survey, category = unname(outcome_label[outcome]),
      mean, median, sd, q2.5, q97.5, rhat, ess_bulk, ess_tail,
      trial_arm = arm,
      prevalence_ratio = fmt_ci(median, q2.5, q97.5, digits = 2)
    )

  # Crude prevalences, compared with the model
  crude_vs_model <- prevalence |>
    dplyr::select(trial_arm, year, category, median, prevalence) |>
    dplyr::left_join(
      describe_carriage_states(dat_arm) |>
        dplyr::filter(measure == "prevalence") |>
        dplyr::transmute(year = survey, category = unname(outcome_label[outcome]),
                         crude_prev = value),
      by = c("year", "category")
    ) |>
    dplyr::mutate(diff_pct_points = (median - crude_prev) * 100)

  list(
    prevalence = prevalence,
    pr = pr,
    crude_vs_model = crude_vs_model,
    diagnostics = diagnostics
  )
}

# Observed density among carriers, median (IQR) on the log10 scale
summarise_density <- function(dat, ...) {
  dat |>
    dplyr::group_by(survey_yearmon, trial_arm, ...) |>
    dplyr::summarise(
      n = dplyr::n(),
      dens_median = median(dens_tot, na.rm = TRUE),
      dens_q1 = quantile(dens_tot, 0.25, na.rm = TRUE),
      dens_q3 = quantile(dens_tot, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      density = fmt_ci(log10(dens_median), log10(dens_q1), log10(dens_q3), digits = 2)
    )
}


############################################################
#### 2) Data ###############################################
############################################################

# Analysis population (00_setup), with its covariate coding
base_dat <- load_analysis_data()$dat$data |>
  dplyr::mutate(trial_arm = forcats::fct_drop(factor(trial_arm)))

table(base_dat$trial_arm, base_dat$survey_yearmon)


############################################################
#### 3) Prevalence by arm ##################################
############################################################

arm_results <- purrr::map(arms, fit_arm)

collect <- function(field) {
  purrr::map(arm_results, \(r) r[[field]]) |>
    purrr::list_rbind() |>
    dplyr::mutate(trial_arm = factor(trial_arm, levels = arms))
}

prevalence <- collect("prevalence")
pr <- collect("pr")
crude_vs_model <- collect("crude_vs_model")
sampler_diagnostics <- collect("diagnostics")
print(sampler_diagnostics)

prevalence_table <- prevalence |>
  dplyr::select(category, trial_arm, year, prevalence) |>
  tidyr::pivot_wider(names_from = year, values_from = prevalence) |>
  dplyr::arrange(category, trial_arm)

pr_table <- pr |>
  dplyr::select(category, trial_arm, year, prevalence_ratio) |>
  tidyr::pivot_wider(names_from = year, values_from = prevalence_ratio) |>
  dplyr::arrange(category, trial_arm)


############################################################
#### 4) Density by arm #####################################
############################################################

data_dens <- base_dat |>
  dplyr::filter(num_sero > 0)

density_overall <- summarise_density(data_dens)
density_by_m <- summarise_density(data_dens, num_sero)

density_overall_table <- density_overall |>
  dplyr::select(trial_arm, survey_yearmon, density) |>
  tidyr::pivot_wider(names_from = survey_yearmon, values_from = density)

density_by_m_table <- density_by_m |>
  dplyr::select(trial_arm, num_sero, survey_yearmon, density) |>
  tidyr::pivot_wider(names_from = survey_yearmon, values_from = density) |>
  dplyr::arrange(num_sero, trial_arm)


############################################################
#### 5) Table S10 and export ###############################
############################################################

s10_measures <- c(
  "Overall pneumococcal carriage prevalence",
  "Multiple serotype carriage prevalence",
  "Overall pneumococcal density (median, IQR, log10 GE per mL)"
)

table_s10 <- dplyr::bind_rows(
  prevalence |>
    dplyr::filter(category %in% c(carriage_label, multiple_label)) |>
    dplyr::transmute(
      measure = dplyr::if_else(category == carriage_label, s10_measures[1], s10_measures[2]),
      trial_arm, year, value = prevalence
    ),
  density_overall |>
    dplyr::transmute(measure = s10_measures[3], trial_arm, year = survey_yearmon, value = density)
) |>
  dplyr::mutate(measure = factor(measure, levels = s10_measures)) |>
  dplyr::arrange(measure, trial_arm, year) |>
  tidyr::pivot_wider(names_from = year, values_from = value)

openxlsx::write.xlsx(
  list(
    Table_S10 = table_s10,
    Prevalence_table = prevalence_table,
    PR_table = pr_table,
    Density_overall_table = density_overall_table,
    Density_by_m_table = density_by_m_table,
    Prevalence = prevalence,
    PR = pr,
    Crude_vs_model = crude_vs_model,
    Density_overall = density_overall,
    Density_by_m = density_by_m,
    Sampler_diagnostics = sampler_diagnostics
  ),
  file = file.path(out_dir, paste0("Sensi_Arm_comparison_", STAMP, ".xlsx")),
  overwrite = TRUE
)


############################################################
#### 6) Figure: prevalence and density by arm ##############
############################################################
# Not in the current manuscript (Table S10 only); kept for reference.

theme_1 <- theme(
  legend.position = "right",
  panel.border = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0, size = 14, vjust = 1),
  plot.title.position = "plot",
  axis.line = element_line(colour = "black", linewidth = 0.5),
  axis.text.x = element_text(colour = "black", size = 10),
  axis.text.y = element_text(colour = "black", size = 10),
  axis.title = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black", linewidth = 0.5)
)

arm_colours <- c("#ccebc5", "#b3cde3", "#decbe4", "#f37974")

plot_prevalence <- function(cat, title, y_max, y_step, y_lab) {
  ggplot(
    prevalence |> dplyr::filter(category == cat),
    aes(x = year, y = median * 100, fill = trial_arm)
  ) +
    annotate("rect", xmin = 5.4, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "grey80", alpha = 0.25) +
    geom_col(position = position_dodge(width = 0.5), width = 0.5,
             color = "black", linewidth = 0.3, alpha = 0.8) +
    geom_errorbar(aes(ymin = q2.5 * 100, ymax = q97.5 * 100),
                  position = position_dodge(width = 0.5),
                  width = 0.2, linewidth = 0.3) +
    geom_vline(xintercept = 4.7 / 3, linetype = "dashed", linewidth = 0.4) +
    labs(title = title, x = "", y = y_lab) +
    scale_fill_manual(values = arm_colours, name = NULL) +
    scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = y_step),
                       expand = c(0, 0)) +
    theme_1
}

p_overall <- plot_prevalence(carriage_label, "A. Overall carriage prevalence",
                             41, 10, "Prevalence (%)")
p_multiple <- plot_prevalence(multiple_label, "B. Multiple serotype carriage prevalence",
                              10, 2, "")

p_density <- ggplot(
  data_dens,
  aes(x = survey_yearmon, y = log10(dens_tot), fill = trial_arm)
) +
  annotate("rect", xmin = 5.5, xmax = Inf, ymin = 0.02, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_violin(position = position_dodge(0.8), alpha = 0.8, color = NA) +
  geom_boxplot(
    aes(group = interaction(survey_yearmon, trial_arm)),
    position = position_dodge(0.8), width = 0.25, outlier.shape = NA,
    fill = NA, color = "grey30"
  ) +
  geom_vline(xintercept = 4.7 / 3, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = arm_colours, name = NULL) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    title = "C. Density by trial arm",
    x = "Survey time",
    y = expression("Density (log"[10]*" GE per mL)")
  ) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), expand = c(0, 0)) +
  theme_1

p_com <- (p_overall + p_multiple + plot_layout(guides = "collect")) / p_density
p_com

fig_stem <- file.path(out_dir, "figures", paste0("NotInManuscript_Prevalence and density by arm_", STAMP))
ggsave(paste0(fig_stem, ".tiff"),
       p_com, width = 10, height = 7, units = "in", dpi = 600, compression = "lzw")

topptx(p_com, paste0(fig_stem, ".pptx"), width = 10, height = 7)

message("\nDone. Written to ", out_dir, "/.")
