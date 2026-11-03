############################################################
# Figure 3: pneumococcal density by survey round and multiplicity
#
# Observed density: the analysis population of 00_setup (pneumoalloc
# analysis data). Model-derived density by multiplicity:
# results/05_density/Density_PR_adj_<STAMP>.xlsx.
# Run from the project directory.
############################################################

#### 0) Load packages ####
rm(list = ls())
library(tidyverse)
library(forcats)
library(openxlsx)
library(cmdstanr)
library(posterior)
library(stringr)   # for str_match
library(patchwork) # for combining ggplots

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

#### 1) Load & prepare core factors ####
source("00_setup/settings_20261015.R")
data <- load_analysis_data()$dat$data %>%     # analysis population
  select(id, survey_yearmon, agemonth,  compchild, trial_arm, income,  famnum_cat, 
         pcv, year, tot, pos, pneu_carr,
         num_sero, sero1, sero2, sero3, num_sero,
         ra_sero1, ra_sero2, ra_sero3,
         dens_sero1, dens_sero2, dens_sero3, dens_tot, log_dens_tot) %>%
  filter(num_sero >= 1) %>% 
  mutate(
    survey_yearmon = as.character(survey_yearmon),
    survey_yearmon = ifelse(survey_yearmon == "Oct. 2017", "Nov. 2017", survey_yearmon) 
  )

data$survey_yearmon <- factor(data$survey_yearmon, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

#### Density: absolute abundance####
data %>%
  group_by(survey_yearmon) %>%
  mutate(dens_total = sum(dens_tot)) %>%
  ungroup() %>%
  filter(pneu_carr == 1) %>%
  # dplyr::select(
  #   id, agemonth, agemonth_z, compchild, trial_arm, income, income_z, famnum_cat, moedu1, pcv,
  #   year, survey_yearmon, tot, pos, num_sero, sero1, sero2, sero3, num_sero1,
  #   ra_sero1, ra_sero2, ra_sero3, dens_sero1, dens_sero2, dens_sero3, dens_tot, dens_total
  # ) %>%
  pivot_longer(
    cols = c(sero1, sero2, sero3),
    names_to = "serotype", values_to = "sero"
  ) %>%
  mutate(sero_all = case_when(
    is.na(sero) ~ "NO",
    sero %in% c("1", "4", "5", "6B", "7F", "9V", "14", "18C", "19F", "23F") ~ "VT",
    TRUE ~ "NVT"
  )) %>%
  filter(!(ra_sero2 == 0 & serotype == "sero2")) %>%
  filter(!(ra_sero3 == 0 & serotype == "sero3")) %>%
  mutate(
    sero = ifelse(is.na(sero), "NO", sero),
    w = case_when(
      serotype == "sero1" ~ ra_sero1,
      serotype == "sero2" ~ ra_sero2,
      serotype == "sero3" ~ ra_sero3
    ),
    w = ifelse(sero_all == "NO", 1, w),
    dens = case_when(
      serotype == "sero1" ~ dens_sero1,
      serotype == "sero2" ~ dens_sero2,
      serotype == "sero3" ~ dens_sero3
    )) -> data_dens

#data_dens$survey_yearmon <- factor(data_dens$survey_yearmon, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))
#### 2) Define VT/NVT sets ####
vt_seros  <- c("1","4","5","6B","7F","9V","14","18C","19F","23F")
nvt_seros <- c("6A", "NT2", "15A", "19A", "15B/C", "11A", "35B", "23A", "34")
keep_seros <- c(vt_seros, nvt_seros)

data_dens %>%
  mutate(
    sero = if_else(sero %in% keep_seros, sero, "Other NVT")
  ) -> data_dens

data_dens$sero <- factor(data_dens$sero,
                         levels = c(
                           "19F", "6B", "23F", "14", "9V",  "1",
                           "6A", "NT2",  "15A", "19A", "15B/C", "11A", "35B", "23A", "34", "13", "3", "6C", "Other NVT"))


#### 3) Density and prevalence by class ####
###mean and median density by all, vt and nvt
data_dens %>% 
  group_by(id, survey_yearmon) %>%
  summarise(dens = sum(dens)) %>%
  group_by(survey_yearmon) %>%
  summarise(mean_dens = mean(dens),
            median_dens = median(dens),
            q25 = quantile(dens, 0.25),
            q75 = quantile(dens, 0.75),
            mean_dens_log = log10(mean(dens)),
            median_dens_log = log10(median(dens)),
            q25_log = log10(q25),
            q75_log = log10(q75),
            mean_log = mean(log10(dens)),
            median_log = median(log10(dens)),
            
  ) %>% 
  mutate(sero_all = "overall")-> dens_all

data_dens %>%
  group_by(id, survey_yearmon, sero_all) %>%
  summarise(dens = sum(dens)) %>%
  group_by(survey_yearmon, sero_all) %>%
  summarise(mean_dens = mean(dens),
            median_dens = median(dens),
            q25 = quantile(dens, 0.25),
            q75 = quantile(dens, 0.75),
            mean_dens_log = log10(mean(dens)),
            median_dens_log = log10(median(dens)),
            q25_log = log10(q25),
            q75_log = log10(q75),
            mean_log = mean(log10(dens)),
            median_log = median(log10(dens))
  ) %>% 
  bind_rows(dens_all) %>% 
  arrange(sero_all)-> dens_class


p_all <- ggplot(data = data, aes(x = survey_yearmon, y = log_dens_tot)) +
  annotate("rect",
           xmin = 5.5, xmax = Inf, ymin = 0.02, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_vline(aes(xintercept = 4.7 / 3),
             linetype = "dashed", linewidth = 0.4)+
  geom_violin(fill = "#F4A3A8", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.25, outlier.shape = NA, 
               fill = "white", color = "grey30") +
  # geom_jitter(color="red", size=0.4, alpha= 0.5,
  #             position = position_jitter(width = 0.1, height = 0)) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(
    limits = c(0, 10), 
    expand = c(0,0),
    name = expression("Density (log"[10]*" GE per mL)"),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  labs(
    x = "", 
    title = "A. All pneumococci"
  ) +
  theme_1+
  theme(
    legend.title = element_blank(),
    #axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(
      size = 14,
      hjust = 0, # Horizontal: 0 = left, 0.5 = center, 1 = right
      vjust = 1 # Vertical: 0 = bottom, 1 = top
    )
  )
p_all

data_dens %>% 
  group_by(id,survey_yearmon, sero_all) %>%
  summarise(dens = sum(dens)) -> data_dens_class

p_vt <- ggplot(data = data_dens_class %>% filter(sero_all == "VT"), aes(x = survey_yearmon, y = log10(dens))) +
  annotate("rect",
           xmin = 5.5, xmax = Inf, ymin = 0.02, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_vline(aes(xintercept = 4.7 / 3),
             linetype = "dashed", linewidth = 0.4)+
  geom_violin(fill = "#F4A3A8", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", color = "grey30") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(
    limits = c(0, 10), 
    expand = c(0,0),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  labs(
    x = "", 
    y = "",
    title = "B. PCV10 serotypes"
  ) +
  theme_1+
  theme(
    legend.title = element_blank(),
    #axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(
      size = 14,
      hjust = 0, # Horizontal: 0 = left, 0.5 = center, 1 = right
      vjust = 1 # Vertical: 0 = bottom, 1 = top
    )
  )
p_vt

p_nvt <- ggplot(data = data_dens_class %>% filter(sero_all == "NVT"), aes(x = survey_yearmon, y = log10(dens))) +
  annotate("rect",
           xmin = 5.5, xmax = Inf, ymin = 0.02, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_vline(aes(xintercept = 4.7 / 3),
             linetype = "dashed", linewidth = 0.4)+
  geom_violin(fill = "#F4A3A8", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", color = "grey30") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(
    limits = c(0, 10), 
    expand = c(0,0),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  labs(
    x = "", 
    y = "",
    title = "C. Non-PCV10 serotypes"
  ) +
  theme_1+
  theme(
    legend.title = element_blank(),
    #axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(
      size = 14,
      hjust = 0, # Horizontal: 0 = left, 0.5 = center, 1 = right
      vjust = 1 # Vertical: 0 = bottom, 1 = top
    )
  )
p_nvt

###density by multiple serotypes###
data %>% 
ungroup() %>%
  group_by(survey_yearmon, num_sero) %>%
  summarise(n = n(),
            mean_dens = mean(dens_tot),
            median_dens = median(dens_tot),
            q25 = quantile(dens_tot, 0.25),
            q75 = quantile(dens_tot, 0.75),
            mean_dens_log = log10(mean(dens_tot)),
            median_dens_log = log10(median(dens_tot)),
            q25_log = log10(q25),
            q75_log = log10(q75),
            mean_log = mean(log10(dens_tot)),
            median_log = median(log10(dens_tot))
  ) |> 
  mutate(
    sero_all = case_when(
      num_sero == 1 ~ "One serotype",
      num_sero == 2 ~ "Two serotypes",
      num_sero == 3 ~ "Three serotypes"
    )
  )-> dens_multi

data %>% 
  filter(num_sero > 1) %>%
  group_by(survey_yearmon) %>%
  summarise(n = n(),
            mean_dens = mean(dens_tot),
            median_dens = median(dens_tot),
            q25 = quantile(dens_tot, 0.25),
            q75 = quantile(dens_tot, 0.75),
            mean_dens_log = log10(mean(dens_tot)),
            median_dens_log = log10(median(dens_tot)),
            q25_log = log10(q25),
            q75_log = log10(q75),
            mean_log = mean(log10(dens_tot)),
            median_log = median(log10(dens_tot))
  ) |> 
  mutate(sero_all = "Multiple serotypes")-> dens_multi1

dens_class |> 
  bind_rows(dens_multi) |> 
  bind_rows(dens_multi1) |> 
  mutate(
    density = paste0(
      formatC(median_log, format = "f", digits = 2), " (",
      formatC(q25_log, format = "f", digits = 2), "-",
      formatC(q75_log, format = "f", digits = 2), ")"
    ),
  )-> dens_class_multi

write.xlsx(dens_class_multi, "results/07_plots/All_density_summary_20261015.xlsx", 
           sheetName = "Density_by_multiplicity", overwrite = TRUE)

table(data$survey_yearmon, data$num_sero)

plot_data <- data %>%
  filter(num_sero > 0) %>%                     # carriers only
  mutate(num_sero_f = factor(pmin(num_sero, 3),
                             levels = 1:3,
                             labels = c("1 serotype", "2 serotypes", "3 serotypes")))

plot_data_f <- plot_data %>%
  group_by(survey_yearmon, num_sero_f) %>%
  filter(n() >= 2) %>%
  ungroup()

plot_data_f <- data %>%
  filter(num_sero > 0) %>%
  mutate(num_sero_f = factor(pmin(as.integer(num_sero), 3),
                             levels = 1:3,
                             labels = c("1 serotype", "2 serotypes", "3 serotypes"))) %>%
  group_by(survey_yearmon, num_sero_f) %>%
  filter(n() >= 2) %>%
  ungroup()

count(plot_data_f, survey_yearmon, num_sero_f)

dodge <- position_dodge(width = 0.8, preserve = "single")

p_multi_observed <- ggplot(plot_data_f, aes(x = survey_yearmon, y = log10(dens_tot), fill = num_sero_f)) +
  annotate("rect", xmin = 5.5, xmax = Inf, ymin = 0.02, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_vline(xintercept = 4.7 / 3, linetype = "dashed", linewidth = 0.4) +
  geom_violin(alpha = 0.8, colour = NA, width = 0.75,
              position = dodge, trim = F) +
  # geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", colour = "grey30",
  #              position = dodge, show.legend = FALSE) +
  scale_fill_manual(values = c("#9EC9E2", "#F4A3A8", "#C7A6D6"), drop = FALSE) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x), drop = FALSE) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0),
                     labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  labs(x = "", y = "", title = "C. Observed density by number of detected serotypes") +
  theme_1 +
  theme(legend.title = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = 14, hjust = 0, vjust = 1))
p_multi_observed

# ggsave("results/07_plots/figures/Fig3C_Observed density by multiplicity_20261015.tiff", p_multi_observed,
#        width = 8, height = 4, units = "in", dpi = 600, compression = "lzw")


### model derived density by class####
dens_m <- read.xlsx("results/05_density/Density_PR_adj_20261015.xlsx", sheet = "Density_by_multiplicity") 

table(dens_m$m)
dens_m$m <- factor(dens_m$m, 
                   levels = c("One serotype", "Two serotypes", "Three serotypes"),
                   labels = c("One serotype", "Two serotypes", "Three serotypes"))
table(dens_m$year)
view(dens_m)

dens_m$year <- factor(dens_m$year, levels = c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"),
                      labels = c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))


p_model <- ggplot(
  dens_m |> filter(group == "Overall"), aes(x = year, fill = m)) +
  annotate("rect",
           xmin = 5.5, xmax = Inf, ymin = 0.02, ymax = Inf,
           fill = "grey80", alpha = 0.25) +
  geom_boxplot(
    aes(
      ymin = log10(median_q025D_by_m),
      lower = q25_log10,
      middle = median_log10,
      upper = q75_log10,
      ymax = log10(median_q975D_by_m)
    ),
    stat = "identity",
    position = position_dodge(width = 0.7),
    alpha = 0.75,
    width = 0.5
  ) +
  geom_vline(aes(xintercept = 4.7 / 3),
             linetype = "dashed", linewidth = 0.4
  )+
  scale_fill_manual(
    values = c( "#9ECAE1", "#CAB2D6", "#FDBB84"),
    labels = c("1" = "Carrying one serotype",
               "2" = "Carrying two serotypes",
               "3" = "Carrying three serotypes"),
    name = ""
  ) +
  labs(
    title = "D. Model-derived density by number of detected serotypes",
    x = "Survey time",
    y = expression("Density (log"[10]*" GE per mL)")  ) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2), expand = c(0,0)) +
  theme_1 +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),             # 2. Controls the text size in the legend
    #legend.key.size = unit(0.4, "cm"),                  # 2. Controls the size of the colored boxes
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    #axis.text.x = element_text(angle = 30, hjust = 1)
  ) 
p_model

library(patchwork)
p_dens1 <- (p_all + p_vt + p_nvt)/ p_model
p_dens1

ggsave("results/07_plots/figures/Fig3_Density_20261015.tiff", p_dens1,
       width = 11, height = 7, units = "in", dpi = 600, compression = "lzw")

library(eoffice)
topptx(p_dens1, "results/07_plots/figures/Fig3_Density_20261015.pptx", width = 11, height = 7)

