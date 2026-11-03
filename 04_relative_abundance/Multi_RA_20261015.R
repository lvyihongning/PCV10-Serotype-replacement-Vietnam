############################################################
# Relative abundance of serotypes in multiple carriage (Figs S2, S4)
#
# A. (Fig S4) Relative abundance of each common serotype among children with
#    multiple carriage, by survey: beta regression per serotype
#    (04_relative_abundance/0_RA_20261015.stan; this model is not part
#    of pneumoalloc).
# B. (Fig S2) Dominant and subdominant serotypes among multiple carriers, Oct.
#    2016 versus Oct. 2019 (crude).
#
# Uses the analysis population and serotype data of
# 00_setup/01_prepare_data (pneumoalloc): detected serotypes, relative
# abundance (normalised to sum to one within each child) and serotype
# categories, the same as in all other analyses.
#
# Run from the project directory (a few minutes on the Mac).
############################################################
rm(list = ls())
library(tidyverse)
library(forcats)
library(openxlsx)
library(stringr)
library(patchwork)
library(ggrepel)
library(cmdstanr)
library(posterior)

source("00_setup/settings_20261015.R")
dir.create("results/04_relative_abundance/figures", recursive = TRUE, showWarnings = FALSE)

#### 2) Load & prepare core factors ####
dat  <- load_analysis_data()$dat
data <- dat$data          # analysis population, one row per child

vt_seros  <- VT_SEROS
nvt_seros <- NVT_SEROS
keep_seros <- c(vt_seros, nvt_seros)

survey_levels <- dat$survey_levels


############################################################
# A) Long data: one row per multiple carrier x detected serotype
############################################################
sl <- dat$sero_long
sl <- sl[dat$m_obs[sl$k] >= 2, ]

data_ra <- tibble(
  id             = data$id[sl$child],
  survey_yearmon = factor(as.character(dat$survey[sl$child]), levels = survey_levels),
  num_sero       = dat$m_obs[sl$k],
  sero_all       = if_else(sl$vaccine_type, "VT", "NVT"),
  sero           = factor(sl$category),       # retained serotype or "Other NVT"
  w              = sl$abundance,              # relative abundance
  rank           = sl$slot                    # slot 1 = dominant
)

data_ra

table(data_ra$survey_yearmon, data_ra$sero)

mod_beta <- cmdstan_model("04_relative_abundance/0_RA_20261015.stan")

common_serotypes <- c("6A", "NT2", "15A", "15B/C", "19F", "23F", "23A", "6B", "19A")

prepare_serotype_data <- function(dat, serotype_name, eps = 1e-6) {
  
  d <- dat %>%
    filter(as.character(sero) == serotype_name) %>%
    filter(!is.na(w)) %>%
    mutate(
      survey_yearmon = trimws(as.character(survey_yearmon))
    )
  
  # Check for unexpected survey labels
  bad_surveys <- setdiff(unique(d$survey_yearmon), survey_levels)
  if (length(bad_surveys) > 0) {
    stop(
      "Unexpected survey labels for serotype ", serotype_name, ": ",
      paste(bad_surveys, collapse = ", ")
    )
  }
  
  # Re-factor explicitly, then drop empty levels for this serotype
  d <- d %>%
    mutate(
      survey_yearmon = factor(survey_yearmon, levels = survey_levels)
    ) %>%
    filter(!is.na(survey_yearmon)) %>%
    droplevels() %>%
    mutate(
      y = pmin(pmax(w, eps), 1 - eps),
      survey_id = as.integer(survey_yearmon)
    ) %>%
    arrange(survey_id)
  
  if (nrow(d) == 0) {
    stop("No data available for serotype ", serotype_name)
  }
  
  if (any(is.na(d$survey_id))) {
    stop("survey_id still contains NA for serotype ", serotype_name)
  }
  
  stan_data <- list(
    N = nrow(d),
    T = nlevels(d$survey_yearmon),
    survey = d$survey_id,
    y = d$y
  )
  
  list(
    data = d,
    stan_data = stan_data,
    survey_levels = levels(d$survey_yearmon)
  )
}

#-----------------------------
# 6) Helper: fit one serotype
#-----------------------------
fit_one_serotype <- function(dat, serotype_name,
                             iter_warmup = 1000,
                             iter_sampling = 2000,
                             chains = 4,
                             parallel_chains = 4,
                             seed = 1234) {
  
  prep <- prepare_serotype_data(dat, serotype_name)
  
  fit <- mod_beta$sample(
    data = prep$stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 500,
    adapt_delta = 0.95
  )
  
  list(
    serotype = serotype_name,
    fit = fit,
    prep = prep
  )
}

#-----------------------------
# 7) Helper: extract posterior summary
#-----------------------------
extract_mu_summary <- function(fit_obj) {
  
  mu_draws <- fit_obj$fit$draws("mu_out") %>%
    as_draws_df()
  
  survey_levels <- fit_obj$prep$survey_levels
  serotype_name <- fit_obj$serotype
  sero_all <- fit_obj$prep$sero_all
  
  out <- map_dfr(seq_along(survey_levels), function(j) {
    x <- mu_draws[[paste0("mu_out[", j, "]")]]
    
    tibble(
      sero = serotype_name,
      ssero_all = sero_all,
      survey_yearmon = survey_levels[j],
      mean_ra = median(x),
      l95 = quantile(x, 0.025),
      u95 = quantile(x, 0.975)
    )
  })
  
  counts <- fit_obj$prep$data %>%
    count(survey_yearmon, name = "n")
  
  out %>%
    left_join(counts, by = "survey_yearmon")
}

#-----------------------------
# 8) Fit all common serotypes
#-----------------------------
fit_list <- vector("list", length(common_serotypes))
names(fit_list) <- common_serotypes

for (s in common_serotypes) {
  cat("Fitting:", s, "\n")
  fit_list[[s]] <- fit_one_serotype(data_ra, s)
}

fit_list[["15A"]]$fit$summary()
fit_list[["15A"]]$fit$diagnostic_summary()
fit_list[["15A"]]$fit$cmdstan_diagnose()

diag_tab <- purrr::map_dfr(fit_list, function(x) {
  s <- x$fit$summary()
  d <- x$fit$diagnostic_summary()
  
  tibble::tibble(
    sero = x$serotype,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    total_divergent = sum(d$num_divergent),
    total_treedepth = sum(d$num_max_treedepth),
    min_ebfmi = min(d$ebfmi, na.rm = TRUE)
  )
})

diag_tab


#-----------------------------
# 9) Combine summaries
#-----------------------------
ra_summary <- map_dfr(fit_list, extract_mu_summary) %>%
  mutate(
    survey_yearmon = factor(survey_yearmon, levels = survey_levels),
    sero = factor(sero, levels = common_serotypes),
    sero_all = sero %in% vt_seros %>% if_else("VT", "NVT"),
    RA = paste0(round(mean_ra, 2), " (", round(l95, 2), "-", round(u95 , 2), ")"),
  ) %>%
  arrange(sero_all, sero, survey_yearmon)
ra_summary

ra_summary %>% 
  select(sero_all, sero,  survey_yearmon, RA) %>%
  pivot_wider(
    names_from = survey_yearmon,
    values_from = RA, 
    values_fill = "0"
  ) %>% 
  rename(Oct_2016 = `Oct. 2016`, June_2017 = `June 2017`, Nov_2017 = `Nov. 2017`, Oct_2018 = `Oct. 2018`,
         Oct_2019 = `Oct. 2019`,  July_2022 = "July 2022") -> ra_summary_wide

# save summary table
write.xlsx(list(
  ra_summary = ra_summary,
  ra_summary_wide = ra_summary_wide), 
           "results/04_relative_abundance/FigS4_Relative abundance summary_20261015.xlsx", overwrite = TRUE)

#-----------------------------
# 10) Plot
#-----------------------------
theme_1 <- theme(
  legend.position = "right",
  panel.border = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 14),
  axis.line = element_line(colour = "black", linewidth = 0.5),
  axis.line.y = element_line(colour = "black", linewidth = 0.5),
  axis.text.x = element_text(colour = "black", size = 10),
  axis.text.y = element_text(colour = "black", size = 10),
  axis.title = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black", linewidth = 0.5)
)

ra_summary$sero_all <- factor(ra_summary$sero_all, levels = c("VT", "NVT"))
ra_summary$sero <- factor(ra_summary$sero, levels = c("19F", "6B", "23F",
                                                      "NT2", "6A", "15A", "15B/C", "23A", "19A"))

ra_summary %>% 
  filter(survey_yearmon %in% c("Oct. 2016", "Oct. 2019")) %>%
  select(sero_all, sero, survey_yearmon, mean_ra, l95, u95) %>%
  pivot_wider(
    id_cols = c(sero_all, sero),
    names_from = survey_yearmon,
    values_from = c(mean_ra, l95, u95)) %>% 
  rename(x_mid = `mean_ra_Oct. 2016`, x_min = `l95_Oct. 2016`, x_max = `u95_Oct. 2016`,
         y_mid = `mean_ra_Oct. 2019`, y_min = `l95_Oct. 2019`, y_max = `u95_Oct. 2019`) -> ra_plot

table(data_ra$sero)
table(ra_summary$sero)

data_ra %>% 
  filter(sero %in% unique(ra_summary$sero)) -> data_ra1

p_ra <- ggplot(ra_summary, aes(x = survey_yearmon, y = mean_ra, group = sero, color = sero_all)) +
  annotate("rect",
           xmin = 5.4, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey35") +
  geom_jitter(
    data = data_ra1, aes(
      x = survey_yearmon,
      y = w,
      color = sero_all ),
    inherit.aes = FALSE,
    width = 0.12,
    height = 0,
    shape = 16,
    size = 1.5,
    alpha = 0.25,
    show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = l95, ymax = u95),
    width = 0.08,
    linewidth = 0.55,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  # geom_line(
  #   linewidth = 0.65,
  #   alpha = 0.95
  # ) +
  geom_point(
    size = 2.4,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c( "VT"  = "#CC6677",  # muted rose
                "NVT" = "#4477AA" ), 
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")
  ) +
  scale_y_continuous(limits = c(0, 1), 
                     labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE))+
  facet_wrap(~ sero) +
  labs(
    x = "Survey time",
    y = "Relative abundance among multiple serotype carriage",
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  theme_1+
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(size = 14), # Change 14 to your desired size
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank()
  )
p_ra


ggsave(
  "results/04_relative_abundance/figures/FigS4_Relative abundance in multiple carriage_20261015.tiff",
  p_ra,
  width = 12, height = 8, units = "in", dpi = 600, compression = "lzw")

library(eoffice)
topptx(p_ra, "results/04_relative_abundance/figures/FigS4_Relative abundance in multiple carriage_20261015.pptx",
       width = 12, height = 8)


#### domianant and subdominant serotypes ####
vt_seros  <- c("1","4","5","6B","7F","9V","14","18C","19F","23F")
nvt_seros <- c("6A", "NT2", "15A", "19A", "15B/C", "11A", "35B", "23A", "34")
keep_seros <- c(vt_seros, nvt_seros)

data %>%
  filter(num_sero >1 ) %>% 
  mutate(sero = if_else(sero1 %in% keep_seros, sero1, "Other NVT")) %>% 
  group_by(year, survey_yearmon) %>%
  mutate(tot_co = n()) %>%
  group_by(year, survey_yearmon, domin_sero, sero) %>%
  summarise(n_dom = n(), tot_co = first(tot_co)) %>%
  mutate(prop_dom = n_dom / tot_co)   -> data_co_domin


##subdominant serotype: sample-weighted pooling for three serotypes
data %>% 
  select(year, survey_yearmon, sero1, sero2, sero3, num_sero) %>%
  filter(num_sero >1 ) %>% 
  group_by(year, survey_yearmon) %>%
  mutate(tot_co = n()) %>%
  filter(num_sero == 3) %>%
  pivot_longer(cols = c(sero2, sero3), names_to = "sero_slot", values_to = "sero") %>%
  mutate(sero = if_else(sero %in% keep_seros, sero, "Other NVT"),
         n = 0.5) -> data_co_subdomin

data %>% 
  select(year, survey_yearmon, sero1, sero2, sero3, num_sero) %>%
  filter(num_sero >1 ) %>% 
  group_by(year, survey_yearmon) %>%
  mutate(tot_co = n()) %>%
  filter(num_sero == 2) %>%
  mutate(sero = if_else(sero2 %in% keep_seros, sero2, "Other NVT"),
         n = 1) %>% 
  bind_rows(data_co_subdomin) %>% 
  group_by(year, survey_yearmon,  sero) %>%
  summarise(n_subdom = sum(n), tot_co = first(tot_co)) %>%
  mutate(prop_subdom = n_subdom / tot_co)  -> data_co_subdomin1


data_co_domin %>%
  left_join(data_co_subdomin1, by = c("year", "survey_yearmon", "sero", "tot_co")) %>% 
  mutate(
    prop_dom = if_else(is.na(prop_dom), 0, prop_dom),
    prop_subdom = if_else(is.na(prop_subdom), 0, prop_subdom)
  ) -> data_co_comb

data_co_comb$domin_sero <- factor(data_co_comb$domin_sero, levels = c("VT", "NVT"))
data_co_comb$sero <- factor(data_co_comb$sero, levels = c(keep_seros, "Other NVT"))


p_co_domin <- ggplot(
  data = data_co_comb %>% 
    filter(survey_yearmon %in% c("Oct. 2016", "Oct. 2019")) %>%
    mutate(
      prop_subdom = ifelse(is.na(prop_subdom), 0, prop_subdom),
      prop_total = prop_dom + prop_subdom)) +
  # Total bar: dominant + subdominant
  geom_col(aes(
    x = sero,
    y = prop_total * 100,
    group = survey_yearmon),
    fill = "grey95",
    color = "black",
    position = position_dodge(width = 0.75),
    width = 0.55,
    linewidth = 0.35 ) +
  # Dominant component
  geom_col(aes(
    x = sero,
    y = prop_dom * 100,
    fill = survey_yearmon,
    group = survey_yearmon ),
    color = "black",
    position = position_dodge(width = 0.75),
    width = 0.55,
    linewidth = 0.35) +
  scale_fill_manual(
    values = c(
      "Oct. 2016" = "#A6CEE3",
      "Oct. 2019" =  "#F4A3A8"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 65),
    breaks = seq(0, 60, by = 20),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Serotype",
    y = "Proportion among multiple serotype carriage (%)",
    fill = NULL) +
  theme_1 +
  theme(
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 10),
    axis.line = element_line(color = "black"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(face = "bold", size = 11)
  )

p_co_domin

ggsave(
  "results/04_relative_abundance/figures/FigS2_Dominant and minority serotypes_20261015.tiff",
  p_co_domin,
  width = 11, height = 6, units = "in", dpi = 600, compression = "lzw")

library(eoffice)
topptx(p_co_domin, "results/04_relative_abundance/figures/FigS2_Dominant and minority serotypes_20261015.pptx",
       width = 11, height = 6)

