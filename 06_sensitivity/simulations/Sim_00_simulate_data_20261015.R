############################################################
# 06_sensitivity/simulations/Sim_00_simulate_data_20261015.R
#
# High-carriage, high multiple-serotype-carriage (MSC) simulations
# (Supplement; Figs S9 and S10) - step 0: simulated populations.
#
# For each simulation and each of Oct. 2016 and Oct. 2019, a population
# of 1,000 children is built with fixed structure:
#   300 non-carriers, 350 single, 245 two-serotype, 105 three-serotype.
# Non-carriers and single-serotype carriers are resampled with
# replacement from the same survey in both simulations.
#
# Scenarios
#   sim1  "empirical resampling"
#         Two-serotype carriers are resampled intact from the same
#         survey (serotypes, relative abundance, density, covariates),
#         so serotype composition is the observed one. Three-serotype
#         carriers extend a resampled two-serotype carrier by a third,
#         lowest-abundance serotype drawn (frequency-weighted) from the
#         minority serotypes of multiple carriers in that survey;
#         abundance profile from observed three-serotype carriers;
#         density = two-serotype density x the adjusted three- versus
#         two-serotype density ratio of the primary density model
#         (05_density, all pneumococci).
#   sim2  "structured minority replacement"
#         MSC children are built on observed backbones (covariates,
#         density, relative abundance) and their serotypes ASSIGNED to
#         meet targets: VT-dominant share fixed at 90% in both periods,
#         so replacement happens only in minority positions; NVT
#         minority share 10% -> 95%, always filling the lowest-abundance
#         positions first. Serotypes are drawn from the Oct. 2016
#         pools in both periods, except minority non-PCV10 serotypes,
#         which come from the same survey (the observed replacement
#         serotypes). Three-serotype carriers: covariates from observed
#         MSC children, abundance AND density from observed
#         three-serotype carriers.
#
# sim1 was revised on 12 October 2026 to match the supplement
# (20261008): it no longer reassigns serotypes, and the 2.85/2.32
# density factor was replaced by the model estimate. sim2 was revised
# on 14 October 2026: the VT-dominant share is held at 90% in Oct. 2019
# (previously 90% -> 45%), so that non-PCV10 replacement is confined to
# minority positions; and its PCV10 and dominant serotypes are drawn
# from the Oct. 2016 pools in both periods. The Oct. 2019 pools hold only
# 4 PCV10-dominant multiple carriers, which had inflated serotype 14.
#
# Input: data/carr_dat_clean_20260917.qs
#        results/05_density/Density_summaries_adj_20261015.rds
#        (run 05_density/Density_20261015.R first; pneumoalloc density model)
#
# Output (results/06_sensitivity/simulations/):
#   Sim_data_20261015.rds               list(sim1 = ..., sim2 = ...)
#   Sim_data_checks_20261015.xlsx       scenario checks
# Next: Sim_01_fit_20261015.R
############################################################

rm(list = ls())
library(qs)
library(tidyverse)
library(openxlsx)


############################################################
#### 0) Settings ###########################################
############################################################

date_tag <- "20261015"
out_dir <- "results/06_sensitivity/simulations"
density_file <- "results/05_density/Density_summaries_adj_20261015.rds"
periods <- c("Oct. 2016", "Oct. 2019")
n_check_reps <- 0L   # > 0: also simulate this many replicates for a stability check

PCV10_SEROTYPES <- c("1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F")

# Columns kept in the simulated data. Other columns of the source rows
# (e.g. pcv10_carr, sero1_3) would be stale for re-assigned MSC children.
keep_cols <- c(
  "id", "survey_yearmon", "commune_ID",
  "agemonth", "compchild", "trial_arm", "income", "famnum_cat", "moedu1", "ill1",
  "num_sero", "dens_tot", "sero1", "sero2", "sero3",
  "ra_sero1", "ra_sero2", "ra_sero3", "dens_sero1", "dens_sero2", "dens_sero3"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
#### 1) Data ###############################################
############################################################

data <- qread("data/carr_dat_clean_20260917.qs") |>
  ungroup() |>
  mutate(
    sero1 = as.character(sero1),
    sero2 = as.character(sero2),
    sero3 = as.character(sero3),
    .source_row = row_number()
  )

# Adjusted density ratio, three versus two detected serotypes, all
# pneumococci (primary density model; common to all survey rounds)
DENS_RATIO_3_VS_2 <- readRDS(density_file)$multiplicity_effects |>
  filter(group == "Overall", estimand == "Density ratio",
         contrast == "Three versus two serotypes", year == "Oct. 2016") |>
  pull(median)
stopifnot(length(DENS_RATIO_3_VS_2) == 1, is.finite(DENS_RATIO_3_VS_2), DENS_RATIO_3_VS_2 > 0)

message("Three- vs two-serotype density ratio (adjusted, posterior median): ",
        round(DENS_RATIO_3_VS_2, 3))


############################################################
#### 2) Scenarios ##########################################
############################################################

# composition: "observed" keeps the serotypes of resampled carriers;
#   "assigned" re-assigns MSC serotypes to the targets below
# minority_rule (assigned only): how minority positions are ordered
# before the first n_NVT of them receive an NVT
#   "random"             random order
#   "vt_dominant_low_ra" VT-dominant children first, then low relative
#                        abundance, mixed with noise by ra_strength
#   "lowest_ra"          strictly lowest relative abundance first
# triple_source: how three-serotype carriers are built
#   "extended_double"    observed two-serotype carrier + one minority
#                        serotype; density x density_ratio_3_vs_2
#   "triple_profile"     covariates from an observed MSC child; abundance
#                        and density from an observed three-serotype carrier
# fixed_pools (assigned only): serotype pools taken from one survey in
#   both periods instead of the same survey (NULL = all same-survey)
scenarios <- list(
  sim1 = list(
    label = "Empirical resampling",
    seed = 20260923,
    composition = "observed",
    vt_dominant = c("Oct. 2016" = NA_real_, "Oct. 2019" = NA_real_),
    nvt_minority = c("Oct. 2016" = NA_real_, "Oct. 2019" = NA_real_),
    minority_rule = c("Oct. 2016" = NA_character_, "Oct. 2019" = NA_character_),
    ra_strength = c("Oct. 2016" = NA_real_, "Oct. 2019" = NA_real_),
    triple_source = "extended_double",
    density_ratio_3_vs_2 = DENS_RATIO_3_VS_2,
    fixed_pools = NULL
  ),
  sim2 = list(
    label = "Structured minority replacement",
    seed = 20260909,
    composition = "assigned",
    vt_dominant = c("Oct. 2016" = 0.90, "Oct. 2019" = 0.90),
    nvt_minority = c("Oct. 2016" = 0.10, "Oct. 2019" = 0.95),
    minority_rule = c("Oct. 2016" = "lowest_ra", "Oct. 2019" = "lowest_ra"),
    ra_strength = c("Oct. 2016" = NA, "Oct. 2019" = NA),
    triple_source = "triple_profile",
    density_ratio_3_vs_2 = NA,
    fixed_pools = list(period = "Oct. 2016",
                       pools = c("dominant_VT", "dominant_NVT", "minority_VT", "all_VT"))
  )
)

population <- list(N = 1000, carriage_prev = 0.70, msc_among_carriers = 0.50,
                   prop_two_among_msc = 0.70)

target_counts <- function(p = population) {
  n_carrier <- round(p$N * p$carriage_prev)
  n_msc <- round(n_carrier * p$msc_among_carriers)
  n_two <- round(n_msc * p$prop_two_among_msc)
  list(noncarrier = p$N - n_carrier, single = n_carrier - n_msc,
       two = n_two, three = n_msc - n_two)
}


############################################################
#### 3) Functions ##########################################
############################################################

# Observed three-serotype relative-abundance profiles (sorted), with
# density for "triple_profile". Pooled over Oct. 2016 - Oct. 2019.
make_triple_pool <- function(dat, require_density) {
  x <- dat |>
    filter(
      as.character(survey_yearmon) %in%
        c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019"),
      num_sero == 3
    )
  if (require_density) x <- filter(x, is.finite(dens_tot), dens_tot > 0)

  ra <- as.matrix(x[c("ra_sero1", "ra_sero2", "ra_sero3")])
  valid <- apply(ra, 1, \(z) all(is.finite(z)) && all(z > 0))
  ra <- t(apply(ra[valid, , drop = FALSE], 1, \(z) sort(z / sum(z), decreasing = TRUE)))
  if (nrow(ra) == 0) stop("No valid observed three-serotype profiles.")

  tibble(ra_rank1 = ra[, 1], ra_rank2 = ra[, 2], ra_rank3 = ra[, 3],
         dens_tot = x$dens_tot[valid])
}

triple_pools <- list(
  extended_double = make_triple_pool(data, require_density = FALSE),
  triple_profile = make_triple_pool(data, require_density = TRUE)
)

serotypes_long <- function(dat) {
  bind_rows(
    transmute(dat, .source_row, num_sero, slot = 1L, serotype = sero1, ra = as.numeric(ra_sero1)),
    transmute(dat, .source_row, num_sero, slot = 2L, serotype = sero2, ra = as.numeric(ra_sero2)),
    transmute(dat, .source_row, num_sero, slot = 3L, serotype = sero3, ra = as.numeric(ra_sero3))
  ) |>
    filter(!is.na(serotype), serotype != "", is.finite(ra))
}

# Same-survey pools of dominant and minority serotypes, VT and NVT.
# Vectors keep repeats, so sampling is frequency-weighted.
make_serotype_pools <- function(dat_period) {
  all_long <- serotypes_long(dat_period)
  ranked <- all_long |>
    filter(num_sero >= 2) |>
    group_by(.source_row) |>
    arrange(desc(ra), slot, .by_group = TRUE) |>
    mutate(rank = row_number()) |>
    ungroup()

  pick <- function(d, vt) d$serotype[(d$serotype %in% PCV10_SEROTYPES) == vt]
  all_vt <- pick(all_long, TRUE)
  all_nvt <- pick(all_long, FALSE)
  or_all <- function(x, fallback) if (length(x) == 0) fallback else x

  list(
    dominant_VT = or_all(pick(filter(ranked, rank == 1), TRUE), all_vt),
    dominant_NVT = or_all(pick(filter(ranked, rank == 1), FALSE), all_nvt),
    minority_VT = or_all(pick(filter(ranked, rank > 1), TRUE), all_vt),
    minority_NVT = or_all(pick(filter(ranked, rank > 1), FALSE), all_nvt),
    minority_all = or_all(filter(ranked, rank > 1)$serotype, all_long$serotype),
    all_VT = all_vt,
    all_NVT = all_nvt
  )
}

sample_unique_serotype <- function(pool, fallback, existing = character()) {
  x <- pool[!(pool %in% existing)]
  if (length(x) == 0) x <- fallback[!(fallback %in% existing)]
  if (length(x) == 0) stop("No unique serotype candidate available.")
  sample(x, 1)
}

# Sort the two relative abundances of a two-serotype backbone, largest first
sort_double_ra <- function(dat) {
  r <- cbind(dat$ra_sero1, dat$ra_sero2)
  ok <- is.finite(r[, 1]) & is.finite(r[, 2]) & rowSums(r) > 0
  r[ok, ] <- t(apply(r[ok, , drop = FALSE], 1, \(z) sort(z / sum(z), decreasing = TRUE)))
  dat$ra_sero1 <- r[, 1]
  dat$ra_sero2 <- r[, 2]
  dat$ra_sero3 <- NA_real_
  dat
}

double_pool <- function(dat_period) {
  pool <- dat_period |>
    filter(num_sero == 2, is.finite(dens_tot), dens_tot > 0,
           is.finite(ra_sero1), is.finite(ra_sero2))
  if (nrow(pool) == 0) stop("No valid observed two-serotype carriers.")
  pool
}

generate_doubles <- function(dat_period, n_two) {
  double_pool(dat_period) |>
    slice_sample(n = n_two, replace = TRUE) |>
    mutate(num_sero = 2L, sero3 = NA_character_, .synthetic_triple = FALSE) |>
    sort_double_ra()
}

generate_triples <- function(dat_period, n_three, sc, pools) {
  pool <- triple_pools[[sc$triple_source]]

  if (sc$triple_source == "extended_double") {
    # Observed two-serotype carrier (covariates, serotypes, density),
    # extended by a third, lowest-abundance serotype drawn from the same
    # survey's minority serotypes; density scaled by the model ratio
    out <- double_pool(dat_period) |> slice_sample(n = n_three, replace = TRUE)
    p <- pool[sample(seq_len(nrow(pool)), n_three, replace = TRUE), ]
    swap <- out$ra_sero2 > out$ra_sero1
    out[swap, c("sero1", "sero2")] <- out[swap, c("sero2", "sero1")]
    out$sero3 <- map2_chr(out$sero1, out$sero2, \(s1, s2) {
      sample_unique_serotype(pools$minority_all, c(pools$all_VT, pools$all_NVT), c(s1, s2))
    })
    out$dens_tot <- out$dens_tot * sc$density_ratio_3_vs_2
  } else {
    # Covariates from an observed MSC child; abundance and density from
    # an observed three-serotype profile; serotypes assigned later
    msc_pool <- filter(dat_period, num_sero >= 2)
    if (nrow(msc_pool) == 0) stop("No observed MSC carriers.")
    out <- slice_sample(msc_pool, n = n_three, replace = TRUE)
    p <- pool[sample(seq_len(nrow(pool)), n_three, replace = TRUE), ]
    out$dens_tot <- p$dens_tot
    out$sero1 <- NA_character_
    out$sero2 <- NA_character_
    out$sero3 <- NA_character_
  }

  out$num_sero <- 3L
  out$ra_sero1 <- p$ra_rank1
  out$ra_sero2 <- p$ra_rank2
  out$ra_sero3 <- p$ra_rank3
  out$.synthetic_triple <- TRUE
  out
}

# Dominant serotype (slot 1): VT for the target share, in each of the
# two- and three-serotype groups
assign_dominant <- function(msc, pools, target_vt) {
  for (k in c(2, 3)) {
    idx <- which(msc$num_sero == k)
    if (length(idx) == 0) next
    n_vt <- round(length(idx) * target_vt)
    type <- sample(c(rep("VT", n_vt), rep("NVT", length(idx) - n_vt)))
    for (j in seq_along(idx)) {
      msc$sero1[idx[j]] <- if (type[j] == "VT") {
        sample_unique_serotype(pools$dominant_VT, pools$all_VT)
      } else {
        sample_unique_serotype(pools$dominant_NVT, pools$all_NVT)
      }
    }
  }
  msc
}

# Minority serotypes (slots 2 and 3): the first n_NVT positions in the
# order set by `rule` receive an NVT, the rest a VT
assign_minority <- function(msc, pools, target_nvt, rule, ra_strength) {
  slots <- bind_rows(
    msc |> mutate(row_id = row_number()) |> filter(num_sero >= 2) |>
      transmute(row_id, slot = 2L, rank = 2L, ra = ra_sero2),
    msc |> mutate(row_id = row_number()) |> filter(num_sero == 3) |>
      transmute(row_id, slot = 3L, rank = 3L, ra = ra_sero3)
  ) |>
    mutate(dominant_is_vt = msc$sero1[row_id] %in% PCV10_SEROTYPES,
           random_tie = runif(n()))

  slots <- switch(
    rule,
    random = arrange(slots, random_tie),
    vt_dominant_low_ra = slots |>
      mutate(ra_priority = ra_strength * (1 - ra) + (1 - ra_strength) * runif(n())) |>
      arrange(desc(dominant_is_vt), desc(ra_priority), random_tie),
    lowest_ra = arrange(slots, ra, desc(rank), random_tie),
    stop("Unknown minority_rule: ", rule)
  )

  n_nvt <- round(nrow(slots) * target_nvt)
  slots$target_type <- if_else(seq_len(nrow(slots)) <= n_nvt, "NVT", "VT")

  for (row_i in unique(slots$row_id)) {
    info <- slots |> filter(row_id == row_i) |> arrange(slot)
    existing <- msc$sero1[row_i]
    for (j in seq_len(nrow(info))) {
      replacement <- if (info$target_type[j] == "NVT") {
        sample_unique_serotype(pools$minority_NVT, pools$all_NVT, existing)
      } else {
        sample_unique_serotype(pools$minority_VT, pools$all_VT, existing)
      }
      msc[[paste0("sero", info$slot[j])]][row_i] <- replacement
      existing <- c(existing, replacement)
    }
  }
  msc
}

simulate_survey <- function(sc, period) {
  dat_period <- filter(data, as.character(survey_yearmon) == period)
  n <- target_counts()

  sample_group <- function(k, size) {
    pool <- filter(dat_period, num_sero == k)
    if (nrow(pool) == 0) stop("No observations with num_sero = ", k)
    pool |>
      slice_sample(n = size, replace = TRUE) |>
      mutate(.synthetic_triple = FALSE, .synthetic_composition = FALSE)
  }

  noncarrier <- sample_group(0, n$noncarrier)
  single <- sample_group(1, n$single)

  pools <- make_serotype_pools(dat_period)
  if (!is.null(sc$fixed_pools)) {
    fixed <- make_serotype_pools(filter(data, as.character(survey_yearmon) == sc$fixed_pools$period))
    pools[sc$fixed_pools$pools] <- fixed[sc$fixed_pools$pools]
  }
  msc <- if (sc$composition == "observed") {
    bind_rows(sample_group(2, n$two),
              generate_triples(dat_period, n$three, sc, pools) |>
                mutate(.synthetic_composition = TRUE))
  } else {
    bind_rows(generate_doubles(dat_period, n$two),
              generate_triples(dat_period, n$three, sc, pools)) |>
      slice_sample(prop = 1) |>
      assign_dominant(pools, sc$vt_dominant[[period]]) |>
      assign_minority(pools, sc$nvt_minority[[period]],
                      sc$minority_rule[[period]], sc$ra_strength[[period]]) |>
      mutate(.synthetic_composition = TRUE)
  }

  bind_rows(noncarrier, single, msc) |>
    slice_sample(prop = 1) |>
    mutate(
      dens_sero1 = if_else(!is.na(sero1) & sero1 != "", dens_tot * ra_sero1, NA_real_),
      dens_sero2 = if_else(!is.na(sero2) & sero2 != "", dens_tot * ra_sero2, NA_real_),
      dens_sero3 = if_else(!is.na(sero3) & sero3 != "", dens_tot * ra_sero3, NA_real_),
      simulated_id = row_number(),
      simulated_period = period
    ) |>
    select(simulated_period, simulated_id, all_of(keep_cols),
           .synthetic_triple, .synthetic_composition)
}

simulate_pair <- function(sc, replicate_id = 0L) {
  set.seed(sc$seed + replicate_id)
  map(periods, \(p) simulate_survey(sc, p)) |>
    list_rbind() |>
    mutate(simulation_rep = replicate_id)
}


# ---------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------

# One row per detected serotype in MSC children, ranked by abundance
msc_rank_table <- function(sim) {
  sim |>
    mutate(row_id = row_number()) |>
    filter(num_sero >= 2) |>
    select(simulated_period, row_id, sero1:sero3, ra_sero1:ra_sero3) |>
    pivot_longer(c(sero1:sero3, ra_sero1:ra_sero3),
                 names_to = c(".value", "slot"), names_pattern = "(sero|ra_sero)(\\d)") |>
    filter(!is.na(sero), sero != "") |>
    group_by(simulated_period, row_id) |>
    arrange(desc(ra_sero), .by_group = TRUE) |>
    mutate(rank = row_number(), is_nvt = !(sero %in% PCV10_SEROTYPES)) |>
    ungroup()
}

check_simulation <- function(sim) {
  structure <- sim |>
    group_by(simulated_period) |>
    summarise(
      N = n(), carriage = mean(num_sero > 0), MSC = mean(num_sero >= 2),
      MSC_among_carriers = sum(num_sero >= 2) / sum(num_sero > 0),
      two_among_MSC = sum(num_sero == 2) / sum(num_sero >= 2),
      .groups = "drop"
    )

  ranks <- msc_rank_table(sim)
  scenario <- ranks |>
    group_by(simulated_period) |>
    summarise(
      n_MSC = n_distinct(row_id),
      prop_VT_dominant_MSC = mean(!is_nvt[rank == 1]),
      n_minority_positions = sum(rank > 1),
      prop_NVT_minority = mean(is_nvt[rank > 1]),
      median_RA_NVT_minority = median(ra_sero[rank > 1 & is_nvt]),
      median_RA_VT_minority = median(ra_sero[rank > 1 & !is_nvt]),
      .groups = "drop"
    )
  by_rank <- ranks |>
    group_by(simulated_period, rank) |>
    summarise(n = n(), prop_NVT = mean(is_nvt), median_RA = median(ra_sero), .groups = "drop")

  # Crude NVT prevalence under approaches 1-3 and presence, per child
  crude <- sim |>
    mutate(row_id = row_number()) |>
    left_join(
      ranks |>
        group_by(row_id) |>
        mutate(w = ra_sero / sum(ra_sero)) |>
        summarise(dom_nvt = is_nvt[rank == 1], prop_nvt = mean(is_nvt),
                  ra_nvt = sum(w[is_nvt]), any_nvt = any(is_nvt), .groups = "drop"),
      by = "row_id"
    ) |>
    mutate(
      single_nvt = num_sero == 1 & !(sero1 %in% PCV10_SEROTYPES),
      `Dominant serotype` = if_else(num_sero == 1, single_nvt, dom_nvt),
      `Equal share` = if_else(num_sero == 1, as.numeric(single_nvt), prop_nvt),
      `Relative abundance` = if_else(num_sero == 1, as.numeric(single_nvt), ra_nvt),
      `Presence-based` = if_else(num_sero == 1, single_nvt, any_nvt)
    ) |>
    mutate(across(c(`Dominant serotype`:`Presence-based`), \(x) as.numeric(coalesce(x, 0)))) |>
    group_by(simulated_period) |>
    summarise(across(c(`Dominant serotype`:`Presence-based`), mean), .groups = "drop") |>
    pivot_longer(-simulated_period, names_to = "method", values_to = "crude_NVT_prev") |>
    pivot_wider(names_from = simulated_period, values_from = crude_NVT_prev) |>
    mutate(crude_NVT_PR = `Oct. 2019` / `Oct. 2016`)

  integrity <- tibble(
    check = c("duplicate serotypes within a child", "RA does not sum to 1",
              "serotype densities do not sum to total"),
    n_children = c(
      sum(apply(sim[c("sero1", "sero2", "sero3")], 1, \(x) anyDuplicated(x[!is.na(x) & x != ""]) > 0)),
      sum(sim$num_sero > 0 &
            abs(rowSums(sim[c("ra_sero1", "ra_sero2", "ra_sero3")], na.rm = TRUE) - 1) > 1e-6),
      sum(sim$num_sero > 0 &
            abs(rowSums(sim[c("dens_sero1", "dens_sero2", "dens_sero3")], na.rm = TRUE) -
                  sim$dens_tot) > 1e-6 * sim$dens_tot)
    )
  )
  stopifnot(all(integrity$n_children == 0))

  list(structure = structure, scenario = scenario, by_rank = by_rank,
       crude_NVT = crude, integrity = integrity)
}


############################################################
#### 4) Simulate ###########################################
############################################################

sim_data <- map(scenarios, simulate_pair)
checks <- map(sim_data, check_simulation)

# Observed MSC composition of the source surveys, for comparison
observed_msc <- data |>
  filter(as.character(survey_yearmon) %in% periods) |>
  select(all_of(keep_cols)) |>
  mutate(simulated_period = as.character(survey_yearmon)) |>
  msc_rank_table() |>
  group_by(simulated_period) |>
  summarise(
    n_MSC = n_distinct(row_id),
    prop_VT_dominant_MSC = mean(!is_nvt[rank == 1]),
    n_minority_positions = sum(rank > 1),
    prop_NVT_minority = mean(is_nvt[rank > 1]),
    .groups = "drop"
  )

message("\n==== Observed MSC composition ====")
print(observed_msc)

iwalk(checks, \(ch, nm) {
  message("\n==== ", nm, ": ", scenarios[[nm]]$label, " ====")
  print(ch$scenario)
  print(ch$crude_NVT)
})

saveRDS(
  list(sim_data = sim_data, scenarios = scenarios, periods = periods,
       pcv10_serotypes = PCV10_SEROTYPES),
  file.path(out_dir, paste0("Sim_data_", date_tag, ".rds"))
)

check_sheets <- imap(checks, \(ch, nm) {
  set_names(ch, paste0(nm, "_", names(ch)))
}) |>
  list_flatten(name_spec = "{inner}")

scenario_table <- imap(scenarios, \(sc, nm) {
  tibble(
    simulation = nm, label = sc$label, seed = sc$seed, period = periods,
    vt_dominant = unname(sc$vt_dominant[periods]),
    nvt_minority = unname(sc$nvt_minority[periods]),
    minority_rule = unname(sc$minority_rule[periods]),
    ra_strength = unname(sc$ra_strength[periods]),
    composition = sc$composition,
    triple_source = sc$triple_source,
    density_ratio_3_vs_2 = sc$density_ratio_3_vs_2,
    fixed_pools = if (is.null(sc$fixed_pools)) NA_character_ else
      paste0(sc$fixed_pools$period, ": ", paste(sc$fixed_pools$pools, collapse = ", "))
  )
}) |>
  list_rbind()

write.xlsx(c(list(Scenarios = scenario_table, Observed_MSC = observed_msc), check_sheets),
           file.path(out_dir, paste0("Sim_data_checks_", date_tag, ".xlsx")),
           overwrite = TRUE)


############################################################
#### 5) Optional: stability across replicates ##############
############################################################
# Crude NVT PRs over n_check_reps replicate populations, to show that
# the single analysed replicate (replicate 0) is typical.

if (n_check_reps > 0) {
  stability <- imap(scenarios, \(sc, nm) {
    map(seq_len(n_check_reps), \(r) {
      check_simulation(simulate_pair(sc, r))$crude_NVT |> mutate(rep = r)
    }) |>
      list_rbind() |>
      group_by(method) |>
      summarise(median_PR = median(crude_NVT_PR),
                q2.5 = quantile(crude_NVT_PR, 0.025),
                q97.5 = quantile(crude_NVT_PR, 0.975), .groups = "drop") |>
      mutate(simulation = nm, .before = 1)
  }) |>
    list_rbind()
  print(stability)
  write.xlsx(stability, file.path(out_dir, paste0("Sim_stability_", date_tag, ".xlsx")),
             overwrite = TRUE)
}

message("\nDone. Written to ", out_dir, "/.")
