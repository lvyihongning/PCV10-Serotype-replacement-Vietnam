############################################################
# Dominant versus subdominant serotypes (SF comment, Results:
# serotype composition before PCV10 introduction)
#
# Dominant serotype carriage:    the dominant serotype of every carrier,
#                                single and multiple (highest relative
#                                abundance = sero1, checked in
#                                00_setup/01_prepare_data); weight 1.
# Subdominant serotype carriage: the other serotypes of multiple
#                                carriers, weighted as in Figure S3: 1
#                                for the second serotype of a two-
#                                serotype carrier, 0.5 each for the second
#                                and third of a three-serotype carrier,
#                                so every multiple carrier counts once.
#
# For each, the three most prevalent non-PCV10 serotypes (ranked within
# that carriage) and their combined share of it; also the PCV10 share
# and the share of every serotype.
#
# Estimates are Bayesian: binomial likelihood (weighted counts for
# subdominant carriage) with a flat Beta(1, 1) prior, i.e. the
# unadjusted intercept-only model, whose posterior is
# Beta(1 + n, 1 + N - n). Reported as posterior median (95% CrI).
#
# Input:  results/00_setup/analysis_data_20261015.rds
# Output: results/04_relative_abundance/Dominant_vs_subdominant_20261027.xlsx
#
# Run from the project directory (seconds; no MCMC).
############################################################
rm(list = ls())
library(tidyverse)
library(openxlsx)

source("00_setup/settings_20261015.R")

out_file <- "results/04_relative_abundance/Dominant_vs_subdominant_20261027.xlsx"

data <- load_analysis_data()$dat$data


############################################################
# 1) Long data: one row per carrier x detected serotype, with weights
############################################################
long <- data |>
  mutate(row = row_number()) |>
  filter(num_sero >= 1) |>
  select(row, survey_yearmon, num_sero, sero1, sero2, sero3) |>
  pivot_longer(c(sero1, sero2, sero3), names_to = "slot", values_to = "serotype") |>
  filter(!is.na(serotype)) |>
  mutate(
    role = if_else(slot == "sero1", "Dominant", "Subdominant"),
    vt   = if_else(serotype %in% VT_SEROS, "PCV10", "Non-PCV10"),
    w    = if_else(role == "Dominant", 1, 1 / (num_sero - 1))
  )

# every carrier: num_sero serotypes, one dominant (weight 1), and the
# subdominant weights of a multiple carrier sum to 1
chk <- long |>
  group_by(row) |>
  summarise(n = n(), num_sero = first(num_sero),
            w_dom = sum(w[role == "Dominant"]), w_sub = sum(w[role == "Subdominant"]))
stopifnot(all(chk$n == chk$num_sero), all(chk$w_dom == 1),
          all(abs(chk$w_sub - (chk$num_sero >= 2)) < 1e-9))


############################################################
# 2) Weighted counts and Bayesian shares
############################################################
post_text <- function(n, N) {
  q <- qbeta(c(0.5, 0.025, 0.975), 1 + n, 1 + N - n) * 100
  sprintf("%.1f%% (95%% CrI %.1f-%.1f)", q[1], q[2], q[3])
}

# N = carriers (dominant) or multiple carriers (subdominant)
count_by <- function(long, weight) {
  long |>
    group_by(survey_yearmon, role, vt, serotype) |>
    summarise(n = sum({{ weight }}), .groups = "drop") |>
    group_by(survey_yearmon, role) |>
    mutate(N = sum(n)) |>
    ungroup()
}

# the three most prevalent non-PCV10 serotypes within each carriage.
# Ties for third place are all listed (tie = TRUE); the share uses the
# three largest counts, so it does not depend on which one is chosen.
top3_of <- function(counts) {
  counts |>
    filter(vt == "Non-PCV10") |>
    group_by(survey_yearmon, role) |>
    arrange(desc(n), .by_group = TRUE) |>
    mutate(rank = min_rank(desc(n))) |>
    summarise(
      top3 = paste0(serotype[rank <= 3], " (", n[rank <= 3], ")", collapse = ", "),
      tie  = sum(rank <= 3) > 3,
      n    = sum(head(n, 3)),
      N    = first(N),
      .groups = "drop"
    ) |>
    mutate(estimate = map2_chr(n, N, post_text))
}

counts <- count_by(long, w)

# a) the three most prevalent non-PCV10 serotypes
top3_share <- top3_of(counts)

# a2) check: every subdominant serotype weighted 1 (per detection, so
#     N = subdominant detections and three-serotype carriers count twice)
top3_share_detection <- top3_of(count_by(mutate(long, w_detection = 1), w_detection))

n_three <- data |>
  filter(num_sero >= 2) |>
  count(survey_yearmon, three_serotypes = num_sero == 3)

# b) PCV10 serotypes
pcv10_share <- counts |>
  group_by(survey_yearmon, role) |>
  summarise(n = sum(n[vt == "PCV10"]), N = first(N), .groups = "drop") |>
  mutate(estimate = map2_chr(n, N, post_text))

# c) every serotype
sero_share <- counts |>
  mutate(estimate = map2_chr(n, N, post_text)) |>
  arrange(survey_yearmon, role, desc(n))


############################################################
# 3) Print the baseline results and write the workbook
############################################################
message("\n==== ", BASELINE, ": three most prevalent non-PCV10 serotypes ====")
print(top3_share |> filter(survey_yearmon == BASELINE) |>
        select(role, top3, tie, n, N, estimate), width = Inf)

message("Check, subdominant serotypes weighted 1 each (per detection):")
print(top3_share_detection |> filter(survey_yearmon == BASELINE, role == "Subdominant") |>
        select(role, top3, tie, n, N, estimate), width = Inf)
message("Multiple carriers with two / three serotypes:")
print(n_three |> filter(survey_yearmon == BASELINE))

message("\n==== ", BASELINE, ": PCV10 serotypes ====")
print(pcv10_share |> filter(survey_yearmon == BASELINE) |>
        select(role, n, N, estimate), width = Inf)

message("\n==== ", BASELINE, ": every serotype ====")
print(sero_share |> filter(survey_yearmon == BASELINE) |>
        select(vt, serotype, role, n, estimate) |>
        pivot_wider(names_from = role, values_from = c(n, estimate)) |>
        arrange(vt, desc(n_Dominant)),
      n = 40, width = Inf)

write.xlsx(
  list(Top3_NVT_share = top3_share, Top3_NVT_share_per_detection = top3_share_detection,
       PCV10_share = pcv10_share, By_serotype = sero_share, Multiple_carriers = n_three),
  out_file, overwrite = TRUE
)
message("\nWritten to ", out_file)
