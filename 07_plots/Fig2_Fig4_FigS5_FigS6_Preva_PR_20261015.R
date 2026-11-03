rm(list = ls())
library(tidyverse)
library(openxlsx)

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


# Load data
#### 1. Unadjusted Prevalence ####
#### 1.1 Prevalence of VT and NVT ####
##unadjusted VT, NVT and Overall##
data_m1 <- read.xlsx("results/02_approaches_1_4/All_M1_unadj_re0_20261015.xlsx", sheet = "Preva_Class_long")
data_m2 <- read.xlsx("results/02_approaches_1_4/All_M2_unadj_re0_20261015.xlsx", sheet = "Preva_Class_long")
data_m3 <- read.xlsx("results/02_approaches_1_4/All_M3_unadj_re0_20261015.xlsx", sheet = "Preva_Class_long")
data_m4 <- read.xlsx("results/02_approaches_1_4/All_M4_unadj_re0_20261015.xlsx", sheet = "Preva_Class_long")

data_m1_adj <- read.xlsx("results/02_approaches_1_4/All_M1_adj_re0_20261015.xlsx", sheet = "Preva_Class_long")
data_m2_adj <- read.xlsx("results/02_approaches_1_4/All_M2_adj_re0_20261015.xlsx", sheet = "Preva_Class_long")
data_m3_adj <- read.xlsx("results/02_approaches_1_4/All_M3_adj_re0_20261015.xlsx", sheet = "Preva_Class_long")
data_m4_adj <- read.xlsx("results/02_approaches_1_4/All_M4_adj_re0_20261015.xlsx", sheet = "Preva_Class_long")

data_m1 %>% 
  bind_rows(data_m2) %>%
  bind_rows(data_m3) %>%
  bind_rows(data_m4) %>% 
  mutate(group = "unadjusted")-> data_prev_unadj

data_m1_adj %>%
  bind_rows(data_m2_adj) %>%
  bind_rows(data_m3_adj) %>%
  bind_rows(data_m4_adj) %>% 
  mutate(group = "adjusted") %>% 
  bind_rows(data_prev_unadj) %>% 
  mutate(median_1 = round(median * 100, 1),
         q2.5_1   = round(q2.5 * 100, 1),
         q97.5_1  = round(q97.5 * 100, 1),
         approach = case_when(
           method == "Method 4" ~ "Approach 1",
           method == "Method 1" ~ "Approach 2",
           method == "Method 2" ~ "Approach 3",
           method == "Method 3" ~ "Approach 4"),
         prevalence = paste0(format(median_1, nsmall = 1), " (", format(q2.5_1, nsmall = 1), "-", format(q97.5_1, nsmall = 1), ")"),
         year = ifelse(year == "Oct. 2017", "Nov. 2017", year)
  ) -> data_prev

data_prev$method <- factor(data_prev$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))

#data_prev$approach <- factor(data_prev$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
data_prev$year <- factor(data_prev$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

data_prev %>% 
  select(year, class, method,  approach, group, prevalence) %>%
  pivot_wider(names_from = year,
              values_from = prevalence) %>% 
  arrange(group, class, approach)-> data_prev_table

view(data_prev_table)

write.xlsx(list(
  Prevalence_VT_NVT = data_prev,
  Prevalence_VT_NVT_Table = data_prev_table
), "results/07_plots/Fig2_Prevalence of VT and NVT_table_20261015.xlsx")

##### Plot unadjusted prevalence VT and NVT #####
## vt
data_prev %>% 
  filter(group == "unadjusted") %>%
  filter(class == "VT") -> data_prev_vt

p_vt <- ggplot(data_prev_vt, aes(x = year, y = median_1, group = method, fill = method)) +
  annotate("rect",
           xmin = 5.4, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.25
  ) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), 
           width = 0.5, color = "black", linewidth = 0.3, alpha = 0.8) +
  geom_errorbar(aes(ymin = q2.5_1, ymax = q97.5_1),
                position = position_dodge(width = 0.5),
                width = 0.2, linewidth = 0.3) +
  labs(
    x = "", 
    y = "Prevalence (%)",
    title = "A. Prevalence of PCV10 serotypes"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0), breaks = seq(0, 40, 10)) +
  scale_fill_manual(
    values = c(
      "#B2DFDB",   "#9ECAE1","#CAB2D6", "#EF8A91",
      "#f37974", "#F4A3A8", "#FDBB84", 
      "#cbd5e8", "#f4cae4", "#fdcdac", "#f37974",
      "#9CD5E6",  "#D0C5DE","#fdcdac","#f37974", "#fbb4b9" , "#fc8d62","#fb9a99",  "#EA8731",
      "#ccebc5", "#b3cde3", "#decbe4", "#f37974","#FBDAC7",
      "#e0f4e7",   "#d6edf5", "#e7dce9", "#f9a19a",
      "#ccebc5", "#b3cde3", "#decbe4","#fbb4ae"
    ), name = NULL,
    labels = c(
      "Method 1" = "Dominant serotype",
      "Method 2" = "Equal share",
      "Method 3" = "Relative abundance",
      "Method 4" = "Serotype-specific density"
    )
  ) +
  theme_1 +
  theme(
    # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    #axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title.position = "plot",     
    plot.title = element_text(
      size = 14,
      hjust = 0, # Horizontal: 0 = left, 0.5 = center, 1 = right
      vjust = 1 # Vertical: 0 = bottom, 1 = top
    ),
    legend.text = element_text(size = 12),             # 2. Controls the text size in the legend
    #legend.key.size = unit(0.4, "cm")                  # 2. Controls the size of the colored boxes
  ) +
  geom_vline(aes(color = "PCV10 introduction", xintercept = 4.7 / 3),
             linetype = "dashed",
             linewidth = 0.4
  ) +
  scale_color_manual(
    name = "",
    values = c("PCV10 introduction" = "black")
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  )
p_vt

table(data_prev_vt$approach)

## nvt
data_prev %>%
  filter(group == "unadjusted") %>%
  filter(class == "NVT") -> data_prev_nvt

p_nvt <- ggplot(data_prev_nvt, aes(x = year, y = median_1, group = method, fill = method)) +
  annotate("rect",
           xmin = 5.4, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.25
  ) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), 
           width = 0.5, color = "black", linewidth = 0.3, alpha = 0.8) +
  geom_errorbar(aes(ymin = q2.5_1, ymax = q97.5_1),
                position = position_dodge(width = 0.5),
                width = 0.2, linewidth = 0.3
  ) +
  labs(
    x = "", y = "",
    title = "B. Prevalence of non-PCV10 serotypes"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0), breaks = seq(0, 40, 10)) +
  scale_fill_manual(
    values = c( 
      "#B2DFDB",   "#9ECAE1","#CAB2D6", "#EF8A91",
      "#ccebc5", "#b3cde3", "#decbe4", "#f37974",
      "#e0f4e7",   "#d6edf5", "#e7dce9", "#f9a19a",
      "#ccebc5", "#b3cde3", "#decbe4","#fbb4ae"
    ), name = NULL,
    labels = c(
      "Method 1" = "Dominant serotype",
      "Method 2" = "Equal share",
      "Method 3" = "Relative abundance",
      "Method 4" = "Serotype-specific density"
    )
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",     
    #axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(
      size = 14,
      hjust = 0, # Horizontal: 0 = left, 0.5 = center, 1 = right
      vjust = 1 # Vertical: 0 = bottom, 1 = top
    ),
    legend.text = element_text(size = 12),             # 2. Controls the text size in the legend
    #legend.key.size = unit(0.4, "cm")                  # 2. Controls the size of the colored boxes
  ) +
  geom_vline(aes(color = "PCV10 introduction", xintercept = 4.7 / 3),
             linetype = "dashed",
             linewidth = 0.4
  ) +
  scale_color_manual(
    name = "",
    values = c("PCV10 introduction" = "black")
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  )
p_nvt


#### 2. Adjsuted PR ####
#### 2.1 Adjusted PR of VT and NVT ####
##unadjusted PR of VT, NVT and Overall##
data_pr_m1 <- read.xlsx("results/02_approaches_1_4/All_M1_unadj_re0_20261015.xlsx", sheet = "PR_Class_long")
data_pr_m2 <- read.xlsx("results/02_approaches_1_4/All_M2_unadj_re0_20261015.xlsx", sheet = "PR_Class_long")
data_pr_m3 <- read.xlsx("results/02_approaches_1_4/All_M3_unadj_re0_20261015.xlsx", sheet = "PR_Class_long")
data_pr_m4 <- read.xlsx("results/02_approaches_1_4/All_M4_unadj_re0_20261015.xlsx", sheet = "PR_Class_long")

##adjusted PR of VT, NVT and Overall##
data_pr_m1_adj <- read.xlsx("results/02_approaches_1_4/All_M1_adj_re0_20261015.xlsx", sheet = "PR_Class_long")
data_pr_m2_adj <- read.xlsx("results/02_approaches_1_4/All_M2_adj_re0_20261015.xlsx", sheet = "PR_Class_long")
data_pr_m3_adj <- read.xlsx("results/02_approaches_1_4/All_M3_adj_re0_20261015.xlsx", sheet = "PR_Class_long")
data_pr_m4_adj <- read.xlsx("results/02_approaches_1_4/All_M4_adj_re0_20261015.xlsx", sheet = "PR_Class_long")

data_pr_m1 %>% 
  bind_rows(data_pr_m2) %>%
  bind_rows(data_pr_m3) %>%
  bind_rows(data_pr_m4) %>% 
  mutate(group = "unadjusted")-> data_pr_unadj

data_pr_m1_adj %>%
  bind_rows(data_pr_m2_adj) %>%
  bind_rows(data_pr_m3_adj) %>%
  bind_rows(data_pr_m4_adj) %>%
  mutate(group = "adjusted") %>%
  bind_rows(data_pr_unadj) %>% 
  mutate(
    median_1 = round(median, 2),
    q2.5_1   = round(q2.5, 2),
    q97.5_1  = round(q97.5, 2),
    PR = paste0(median_1, " (", q2.5_1, "-", q97.5_1, ")"),
    #PR = paste0(format(median_1, nsmall = 2), " (", format(q2.5_1, nsmall = 2), "-", format(q97.5_1, nsmall = 2), ")"),
    approach = case_when(
      method == "Method 4" ~ "Approach 1",
      method == "Method 1" ~ "Approach 2",
      method == "Method 2" ~ "Approach 3",
      method == "Method 3" ~ "Approach 4"
    ),
    year = ifelse(year == "Oct. 2017", "Nov. 2017", year)
  )-> data_pr


data_pr$method <- factor(data_pr$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))

#data_pr$approach <- factor(data_pr$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
data_pr$year <- factor(data_pr$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

data_pr %>% 
  select(year, class, approach,method,  group, PR) %>%
  pivot_wider(names_from = year,
              values_from = PR) %>% 
  arrange(group, class, method)-> data_pr_table
view(data_pr_table)

write.xlsx(list(
  PR_VT_NVT = data_pr,
  PR_VT_NVT_Table = data_pr_table
), "results/07_plots/Fig2_PR of VT and NVT_table_20261015.xlsx")

##### Plot adjusted PR VT and NVT ####
dodge <- position_dodge(width = 0.5)

pr_vt <- ggplot(
  data = data_pr %>% filter(class == "VT" & group == "adjusted"),
  aes(
    x = year, y = median_1,
    ymin = q2.5_1, ymax = q97.5_1, group = method, color = method
  )) +
  annotate("rect",
           xmin = 4.4, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey80", alpha = 0.25
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  geom_point(position = dodge, size = 1) +
  geom_linerange(position = dodge, linewidth = 0.5) +
  labs(
    title = "C. Adjusted PRs of PCV10 serotypes",
    x = "Survey time",
    y = "Prevalence ratio"
  ) +
  scale_color_manual(
    values = c( 
      "#5AB4AC",  "#4393C3", "#9970AB", "#D6604D",
      "#B2DFDB",   "#9ECAE1","#CAB2D6", "#EF8A91",
      "#99CD4E",  "#5ea5f4", "#B1AFFF","#f15d58", "black", "#5cdbdf", "#FFD1D1", "#789DBC", "#B1AFFF"),
    labels = c(
      "Method 1" = "Dominant serotype",
      "Method 2" = "Equal share",
      "Method 3" = "Relative abundance",
      "Method 4" = "Serotype-specific density"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_log10(
    limits = c(0.05, 1.1),
    breaks = c(0.05, 0.1, 0.5,  1),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE))+
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    #axis.text.x = element_text(angle = 30, hjust = 1)
    legend.text = element_text(size = 12),             # 2. Controls the text size in the legend
    #legend.key.size = unit(0.4, "cm")                  # 2. Controls the size of the colored boxes
  )
pr_vt

pr_nvt <- ggplot(
  data = data_pr %>% filter(class == "NVT" & group == "adjusted"),
  aes(
    x = year, y = median_1,
    ymin = q2.5_1, ymax = q97.5_1, group = method, color = method
  )
) +
  annotate("rect",
           xmin = 4.4, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey80", alpha = 0.25
  ) +
  geom_point(position = dodge, size = 1) +
  geom_linerange(position = dodge, linewidth = 0.5) +
  labs(
    title = "D. Adjusted PRs of non-PCV10 serotypes",
    x = "Survey time",
    y = "" ) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
  scale_color_manual(
    values = c(
      "#5AB4AC",  "#4393C3", "#9970AB", "#D6604D",
      "#99CD4E",  "#5ea5f4", "#B1AFFF","#f15d58", "black", "#5cdbdf", "#FFD1D1", "#789DBC", "#B1AFFF"),
    labels = c(
      "Method 1" = "Dominant serotype",
      "Method 2" = "Equal share",
      "Method 3" = "Relative abundance",
      "Method 4" = "Serotype-specific density"
    ),
    name = NULL
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  # scale_y_continuous(
  #   #expand = c(0, 0),
  #   limits = c(1, 2), breaks = seq(0, 10, 0.2)) +
  scale_y_log10(
    limits = c(0.9, 2.1),
    breaks = c(0.8, 1, 1.5, 2),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    #axis.text.x = element_text(angle = 30, hjust = 1)
    legend.text = element_text(size = 12),             # 2. Controls the text size in the legend
    #legend.key.size = unit(0.4, "cm")                  # 2. Controls the size of the colored boxes
  )
pr_nvt

library(patchwork)

p_all <- (p_vt+p_nvt)/(pr_vt+pr_nvt) + plot_layout(guides = "collect") + plot_layout(ncol = 1)
p_all

ggsave("results/07_plots/figures/Fig2_Prevalence and PR of VT and NVT_20261015.tiff", p_all,
       width = 12, height = 7, units = "in", dpi = 600, compression = "lzw")

library(eoffice)
topptx(p_all, "results/07_plots/figures/Fig2_Prevalence and PR of VT and NVT_20261015.pptx",
       width = 12, height = 7)


### 1.2 Serotype prevalence ####
##unadjusted serotype prevalence##
data_sero_m1 <- read.xlsx("results/02_approaches_1_4/All_M1_unadj_re0_20261015.xlsx", sheet = "Preva_Sero_long")
data_sero_m2 <- read.xlsx("results/02_approaches_1_4/All_M2_unadj_re0_20261015.xlsx", sheet = "Preva_Sero_long")
data_sero_m3 <- read.xlsx("results/02_approaches_1_4/All_M3_unadj_re0_20261015.xlsx", sheet = "Preva_Sero_long")
data_sero_m4 <- read.xlsx("results/02_approaches_1_4/All_M4_unadj_re0_20261015.xlsx", sheet = "Preva_Sero_long")

data_sero_m1_adj <- read.xlsx("results/02_approaches_1_4/All_M1_adj_re0_20261015.xlsx", sheet = "Preva_Sero_long")
data_sero_m2_adj <- read.xlsx("results/02_approaches_1_4/All_M2_adj_re0_20261015.xlsx", sheet = "Preva_Sero_long")
data_sero_m3_adj <- read.xlsx("results/02_approaches_1_4/All_M3_adj_re0_20261015.xlsx", sheet = "Preva_Sero_long")
data_sero_m4_adj <- read.xlsx("results/02_approaches_1_4/All_M4_adj_re0_20261015.xlsx", sheet = "Preva_Sero_long")

data_sero_m1 %>%
  bind_rows(data_sero_m2) %>%
  bind_rows(data_sero_m3) %>%
  bind_rows(data_sero_m4) %>%
  mutate(group = "unadjusted")-> data_sero_prev_unadj

data_sero_m1_adj %>%
  bind_rows(data_sero_m2_adj) %>%
  bind_rows(data_sero_m3_adj) %>%
  bind_rows(data_sero_m4_adj) %>%
  mutate(group = "adjusted") %>%
  bind_rows(data_sero_prev_unadj) %>%
  rename(sero = sero_group) %>%
  mutate(median_1 = round(median * 100, 1),
         q2.5_1   = round(q2.5 * 100, 1),
         q97.5_1  = round(q97.5 * 100, 1),
         approach = case_when(
           method == "Method 4" ~ "Approach 1",
           method == "Method 1" ~ "Approach 2",
           method == "Method 2" ~ "Approach 3",
           method == "Method 3" ~ "Approach 4"),
         prevalence = paste0(format(median_1, nsmall = 1), " (", format(q2.5_1, nsmall = 1), "-", format(q97.5_1, nsmall = 1), ")"),
         year = ifelse(year == "Oct. 2017", "Nov. 2017", year))-> data_sero_prev

#data_sero_prev$approach <- factor(data_sero_prev$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
data_sero_prev$year <- factor(data_sero_prev$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

str(data_sero_prev)

data_sero_prev %>%
  select(year, serotype, sero,method,  approach, group, prevalence) %>%
  pivot_wider(names_from = year,
              values_from = prevalence) %>%
  arrange(group, sero, serotype, approach) -> data_sero_prev_table
view(data_sero_prev_table)

write.xlsx(list(
  Serotype_Prevalence = data_sero_prev,
  Serotype_Prevalence_Table = data_sero_prev_table
), "results/07_plots/Fig4_Serotype prevalence_table_20261015.xlsx")

####pre-vaccination carriage in October 2016####
data_prev %>% 
  select(class, method,  prevalence, group, year) %>%
  mutate(sero = class) %>% 
  bind_rows(data_sero_prev %>% 
              select(sero, serotype, method,  prevalence, group, year) %>% 
              rename(class = serotype)) %>%
  filter(group == "unadjusted") %>%
  filter(year == "Oct. 2016") %>% 
  filter(class != "NO") %>% 
  pivot_wider(
    names_from = method,
    values_from = prevalence
  ) %>% 
  arrange(desc(sero))-> data_prev_sero_2016

write.xlsx(data_prev_sero_2016, "results/07_plots/TableS4_Prevalence in October 2016_20261015.xlsx")


###Prevalence in method 4 and method 1 in October 2016 and October 2019###
data_sero_prev %>% 
  filter(method %in% c( "Method 4", "Method 1")) %>% 
  filter(group == "unadjusted") %>%
  filter(year %in% c("Oct. 2016",  "Oct. 2019")) %>% 
  select(sero, serotype, year, method, approach, group, prevalence, median_1, q2.5_1, q97.5_1) -> data_sero_prev_m4_1


# 2. Reshape data: Include BOTH Method 1 and Method 4
data_plot_wide <- data_sero_prev_m4_1 %>%
  filter(serotype != "NO") %>% 
  # Remove the filter step so we keep both methods
  pivot_wider(
    id_cols = c(serotype, sero, method), # Add method to id_cols
    names_from = year,
    values_from = c(median_1, q2.5_1, q97.5_1)
  ) %>%
  rename(
    x_mid = `median_1_Oct. 2016`, x_min = `q2.5_1_Oct. 2016`, x_max = `q97.5_1_Oct. 2016`,
    y_mid = `median_1_Oct. 2019`, y_min = `q2.5_1_Oct. 2019`, y_max = `q97.5_1_Oct. 2019`
  )

data_plot_wide$sero <- factor(data_plot_wide$sero, levels = c("VT", "NVT"))

# 3. Create the plot
library(ggrepel)
plot1 <- ggplot(data_plot_wide, aes(x = x_mid, y = y_mid, color = sero, shape = method)) +
  # Reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_line(aes(group = serotype), color = "gray65", alpha = 0.7, linewidth = 0.35, show.legend = FALSE, linetype = "dashed") +
  # Error bars: Filtered to ONLY show for Method 4
  geom_errorbarh(
    data = subset(data_plot_wide, method == "Method 4"), 
    aes(xmin = x_min, xmax = x_max), 
    alpha = 0.6,
    linewidth = 0.5,
    width = 0,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = subset(data_plot_wide, method == "Method 4"), 
    aes(ymin = y_min, ymax = y_max), 
    alpha = 0.6,
    linewidth = 0.5,
    width = 0,
    show.legend = FALSE
  )+
  # Points: Draws both Method 1 (hollow) and Method 4 (solid)
  geom_point(size = 2.6, stroke = 0.65) + 
  # Labels: Filtered to only Method 4 to prevent duplicate labels
  geom_text_repel(
    data = subset(data_plot_wide, method == "Method 4"), 
    aes(label = serotype), 
    size = 3.5, 
    show.legend = FALSE,
    box.padding = 0.3,
    max.overlaps = Inf ) +
  # Aesthetics for Color (VT vs NVT)
  scale_color_manual(
    values = c( "VT"  = "#CC6677",  # muted rose
                "NVT" = "#4477AA" ), 
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")
  ) +
  # 2. Exact match for shape aesthetic from plot1 (Removed name = NULL)
  scale_shape_manual(
    values = c("Method 4" = 16, "Method 1" = 1),
    labels = c("Method 4" = "Serotype-specific density", "Method 1" = "Dominant serotype")
  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 3, 6, 15), # Adjust breaks to log-friendly numbers
    limits = c(0.05, 15),         # Log can't start at 0; use a small floor
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  scale_y_log10(
    breaks = c(0.1, 1, 3, 6, 15), # Adjust breaks to log-friendly numbers
    limits = c(0.05, 15),         # Log can't start at 0; use a small floor
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  labs(
    title = "A. Prevalence",
    x = "Prevalence in Oct. 2016 (%)",
    y = "Prevalence in Oct. 2019 (%)"
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical", 
    axis.text = element_text(color = "black") )+
  coord_cartesian(
    xlim = c(0.04, 18),
    ylim = c(0.04, 18),
    clip = "off"
  ) 

plot1

###Prevalence in method 2 and method 1 in October 2016 and October 2019###
data_sero_prev %>% 
  filter(method %in% c( "Method 2", "Method 1")) %>% 
  filter(group == "unadjusted") %>%
  filter(year %in% c("Oct. 2016",  "Oct. 2019")) %>% 
  select(sero, serotype, year, method, approach, group, prevalence, median_1, q2.5_1, q97.5_1) -> data_sero_prev_m2_1


# 2. Reshape data: Include BOTH Method 1 and Method 4
data_plot_wide_m2 <- data_sero_prev_m2_1 %>%
  filter(serotype != "NO") %>% 
  # Remove the filter step so we keep both methods
  pivot_wider(
    id_cols = c(serotype, sero, method), # Add method to id_cols
    names_from = year,
    values_from = c(median_1, q2.5_1, q97.5_1)
  ) %>%
  rename(
    x_mid = `median_1_Oct. 2016`, x_min = `q2.5_1_Oct. 2016`, x_max = `q97.5_1_Oct. 2016`,
    y_mid = `median_1_Oct. 2019`, y_min = `q2.5_1_Oct. 2019`, y_max = `q97.5_1_Oct. 2019`
  )

data_plot_wide_m2$sero <- factor(data_plot_wide_m2$sero, levels = c("VT", "NVT"))

# 3. Create the plot
library(ggrepel)
plot1_m2 <- ggplot(data_plot_wide_m2, aes(x = x_mid, y = y_mid, color = sero, shape = method)) +
  # Reference line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_line(aes(group = serotype), color = "gray65", alpha = 0.7, linewidth = 0.35, show.legend = FALSE, linetype = "dashed") +
  # Error bars: Filtered to ONLY show for Method 4
  geom_errorbarh(
    data = subset(data_plot_wide_m2, method == "Method 2"), 
    aes(xmin = x_min, xmax = x_max), 
    alpha = 0.6,
    linewidth = 0.5,
    width = 0,
    show.legend = FALSE ) +
  geom_errorbar(
    data = subset(data_plot_wide_m2, method == "Method 2"), 
    aes(ymin = y_min, ymax = y_max), 
    alpha = 0.6,
    linewidth = 0.5,
    width = 0,
    show.legend = FALSE )+
  geom_point(size = 2.6, stroke = 0.65) + 
  geom_text_repel(
    data = subset(data_plot_wide_m2, method == "Method 2"), 
    aes(label = serotype), 
    size = 3.5, 
    show.legend = FALSE,
    box.padding = 0.3,
    max.overlaps = Inf ) +
  scale_color_manual(
    values = c( "VT"  = "#CC6677",  # muted rose
                "NVT" = "#4477AA" ), 
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")) +
  scale_shape_manual(
    values = c("Method 2" = 16, "Method 1" = 1),
    labels = c("Method 2" = "Equal share", "Method 1" = "Dominant serotype")  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 3, 6, 15), # Adjust breaks to log-friendly numbers
    limits = c(0.05, 15),         # Log can't start at 0; use a small floor
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  scale_y_log10(
    breaks = c(0.1, 1, 3, 6, 15), # Adjust breaks to log-friendly numbers
    limits = c(0.05, 15),         # Log can't start at 0; use a small floor
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  labs(
    title = "A. Prevalence",
    x = "Prevalence in Oct. 2016 (%)",
    y = "Prevalence in Oct. 2019 (%)" ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical", 
    axis.text = element_text(color = "black") )+
  coord_cartesian(
    xlim = c(0.04, 18),
    ylim = c(0.04, 18),
    clip = "off"
  ) 

plot1_m2


###Prevalence in method 3 and method 1 in October 2016 and October 2019###
data_sero_prev %>%
  filter(method %in% c( "Method 3", "Method 1")) %>% 
  filter(group == "unadjusted") %>%
  filter(year %in% c("Oct. 2016",  "Oct. 2019")) %>% 
  select(sero, serotype, year, method, approach, group, prevalence, median_1, q2.5_1, q97.5_1) -> data_sero_prev_m3_1

# 2. Reshape data: Include BOTH Method 1 and Method 4
data_plot_wide_m3 <- data_sero_prev_m3_1 %>%
  filter(serotype != "NO") %>%
  # Remove the filter step so we keep both methods
  pivot_wider(
    id_cols = c(serotype, sero, method), # Add method to id_cols
    names_from = year,
    values_from = c(median_1, q2.5_1, q97.5_1)
  ) %>%
  rename(
    x_mid = `median_1_Oct. 2016`, x_min = `q2.5_1_Oct. 2016`, x_max = `q97.5_1_Oct. 2016`,
    y_mid = `median_1_Oct. 2019`, y_min = `q2.5_1_Oct. 2019`, y_max = `q97.5_1_Oct. 2019`
  )

data_plot_wide_m3$sero <- factor(data_plot_wide_m3$sero, levels = c("VT", "NVT"))

# 3. Create the plot
plot1_m3 <- ggplot(data_plot_wide_m3, aes(x = x_mid, y = y_mid, color = sero, shape = method)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey20", linewidth = 0.4) +
  geom_line(aes(group = serotype), color = "gray65", alpha = 0.7, linewidth = 0.35, show.legend = FALSE, linetype = "dashed") +
  # Error bars: Filtered to ONLY show for Method 4
  geom_errorbarh(
    data = subset(data_plot_wide_m3, method == "Method 3"), 
    aes(xmin = x_min, xmax = x_max), 
    alpha = 0.6,
    linewidth = 0.5,
    width = 0,
    show.legend = FALSE ) +
  geom_errorbar(
    data = subset(data_plot_wide_m3, method == "Method 3"), 
    aes(ymin = y_min, ymax = y_max), 
    alpha = 0.6,
    linewidth = 0.5,
    width = 0,
    show.legend = FALSE )+
  geom_point(size = 2.6, stroke = 0.65) + 
  geom_text_repel(
    data = subset(data_plot_wide_m3, method == "Method 3"), 
    aes(label = serotype), 
    size = 3.5, 
    show.legend = FALSE,
    box.padding = 0.3,
    max.overlaps = Inf ) +
  scale_color_manual(
    values = c( "VT"  = "#CC6677",  # muted rose
                "NVT" = "#4477AA" ), 
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")) +
  scale_shape_manual(
    values = c("Method 3" = 16, "Method 1" = 1),
    labels = c("Method 3" = "Relative abundance", "Method 1" = "Dominant serotype")  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 3, 6, 15), # Adjust breaks to log-friendly numbers
    limits = c(0.05, 15),         # Log can't start at 0; use a small floor
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  scale_y_log10(
    breaks = c(0.1, 1, 3, 6, 15), # Adjust breaks to log-friendly numbers
    limits = c(0.05, 15),         # Log can't start at 0; use a small floor
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  labs(
    title = "A. Prevalence",
    x = "Prevalence in Oct. 2016 (%)",
    y = "Prevalence in Oct. 2019 (%)" ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical", 
    axis.text = element_text(color = "black") )+
  coord_cartesian(
    xlim = c(0.04, 18),
    ylim = c(0.04, 18),
    clip = "off"
  ) 
plot1_m3

#### 2. Adjsuted PR ####
### 2.2 Adjusted PR of individual serotypes ####
##unadjusted PR of individual serotypes##
pr_sero_m1 <- read.xlsx("results/02_approaches_1_4/All_M1_unadj_re0_20261015.xlsx", sheet = "PR_Sero_long")
pr_sero_m2 <- read.xlsx("results/02_approaches_1_4/All_M2_unadj_re0_20261015.xlsx", sheet = "PR_Sero_long")
pr_sero_m3 <- read.xlsx("results/02_approaches_1_4/All_M3_unadj_re0_20261015.xlsx", sheet = "PR_Sero_long")
pr_sero_m4 <- read.xlsx("results/02_approaches_1_4/All_M4_unadj_re0_20261015.xlsx", sheet = "PR_Sero_long")

##adjusted PR of individual serotypes##
pr_sero_m1_adj <- read.xlsx("results/02_approaches_1_4/All_M1_adj_re0_20261015.xlsx", sheet = "PR_Sero_long")
pr_sero_m2_adj <- read.xlsx("results/02_approaches_1_4/All_M2_adj_re0_20261015.xlsx", sheet = "PR_Sero_long")
pr_sero_m3_adj <- read.xlsx("results/02_approaches_1_4/All_M3_adj_re0_20261015.xlsx", sheet = "PR_Sero_long")
pr_sero_m4_adj <- read.xlsx("results/02_approaches_1_4/All_M4_adj_re0_20261015.xlsx", sheet = "PR_Sero_long")

pr_sero_m1 %>%
  bind_rows(pr_sero_m2) %>%
  bind_rows(pr_sero_m3) %>%
  bind_rows(pr_sero_m4) %>%
  mutate(group = "unadjusted")-> pr_sero_unadj

pr_sero_m1_adj %>%
  bind_rows(pr_sero_m2_adj) %>%
  bind_rows(pr_sero_m3_adj) %>%
  bind_rows(pr_sero_m4_adj) %>%
  mutate(group = "adjusted") %>%
  bind_rows(pr_sero_unadj) %>%
  rename(sero = sero_group) %>%
  mutate(
    median_1 = round(median, 2),
    q2.5_1   = round(q2.5, 2),
    q97.5_1  = round(q97.5, 2),
    PR = paste0(median_1, " (", q2.5_1, "-", q97.5_1, ")"),
    #PR = paste0(format(median_1, nsmall = 2), " (", format(q2.5_1, nsmall = 2), "-", format(q97.5_1, nsmall = 2), ")"),
    approach = case_when(
      method == "Method 4" ~ "Approach 1",
      method == "Method 1" ~ "Approach 2",
      method == "Method 2" ~ "Approach 3",
      method == "Method 3" ~ "Approach 4"
    )
  )-> pr_sero

pr_sero$approach <- factor(pr_sero$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
#pr_sero$method <- factor(pr_sero$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))
pr_sero$year <- factor(pr_sero$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

pr_sero %>%
  select(year, serotype, sero, method, approach, group, PR) %>%
  pivot_wider(names_from = year,
              values_from = PR) %>%
  arrange(group, sero, serotype,  approach) -> pr_sero_table

write.xlsx(list(
  Serotype_PR       = pr_sero,
  Serotype_PR_Table = pr_sero_table
), "results/07_plots/Fig4_Serotype PR_table_20261015.xlsx")

##### Plot adjusted PR of individual serotypes ####
pr_sero %>%
  filter(sero != "NO") %>%
  filter(group == "adjusted") %>%
  #filter(method %in% c("Method 4", "Method 1")) %>%
  filter(year == "Oct. 2019") %>%
  select(sero, serotype, year, method, approach, group, PR, median_1, q2.5_1, q97.5_1) -> pr_sero_m4_adj1
#filter(serotype %in% c("19F", "6B", "23F", "15B/C", "6A", "15A", "23A", "NT2"))
#filter(!(year %in% c("Oct. 2020", "July 2022")))

pr_sero_m4_adj1$sero <- factor(pr_sero_m4_adj1$sero, levels =  c("VT", "NVT"))

pr_sero_m4_adj1$serotype <- factor(pr_sero_m4_adj1$serotype,
                                   levels = c(
                                     "19F", "6B", "23F", "14", "9V",  "1",
                                     "15B/C", "6A",  "15A", "19A", "11A", "35B", "23A","NT2",  "34", "13", "3", "6C", "Other NVT"))

########################################################################
### Figure 1: Method 4 vs Method 1 in October 2019 #####################
########################################################################
lower_cap   <- 0.01
upper_cap   <- 150
dodge_width <- 0.55
n_mod       <- 2

pr_sero_plot_data <- pr_sero_m4_adj1 |>
  filter(method %in% c("Method 1", "Method 4")) |>
  mutate(
    serotype = droplevels(serotype),
    method   = factor(method, levels = c("Method 1", "Method 4")),
    # manual dodge: Method 1 on the left, Method 4 on the right
    x_pos    = as.integer(serotype) +
      (as.integer(method) - (n_mod + 1) / 2) * dodge_width / n_mod,
    below_cap  = !is.na(q2.5_1)  & q2.5_1  < lower_cap,
    over_cap   = !is.na(q97.5_1) & q97.5_1 > upper_cap,
    lower_plot = pmin(pmax(q2.5_1,  lower_cap), upper_cap),
    upper_plot = pmin(pmax(q97.5_1, lower_cap), upper_cap),
    median_plot = ifelse(median_1 < lower_cap | median_1 > upper_cap, NA, median_1),
    arrow_end = case_when(
      below_cap & over_cap ~ "both",
      below_cap            ~ "first",   # arrow pointing down
      over_cap             ~ "last",    # arrow pointing up
      TRUE                 ~ "none"
    )
  )

sero_levels <- levels(pr_sero_plot_data$serotype)

pr_sero_plot <- ggplot(pr_sero_plot_data, aes(x = x_pos, color = sero)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             linewidth = 0.4, color = "grey30") +
  # CIs fully inside the axis
  geom_segment(
    data = filter(pr_sero_plot_data, arrow_end == "none"),
    aes(xend = x_pos, y = lower_plot, yend = upper_plot),
    linewidth = 0.45, alpha = 0.75, show.legend = FALSE, na.rm = TRUE
  ) +
  # CIs cut off at the bottom, top, or both
  lapply(c("first", "last", "both"), function(e) {
    geom_segment(
      data = filter(pr_sero_plot_data, arrow_end == e),
      aes(xend = x_pos, y = lower_plot, yend = upper_plot),
      linewidth = 0.45, alpha = 0.75, show.legend = FALSE,
      arrow = arrow(length = unit(1.2, "mm"), type = "closed", ends = e)
    )
  }) +
  geom_point(aes(y = median_plot, shape = method),
             size = 2.6, stroke = 0.65, na.rm = TRUE) +
  scale_x_continuous(
    breaks = seq_along(sero_levels),
    labels = gsub(" ", "\n", sero_levels),   # "Other NVT" on two lines
    expand = expansion(add = 0.5)
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 150),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  coord_cartesian(ylim = c(lower_cap, upper_cap), clip = "off") +
  labs(
    x = "Serotype",
    y = "Prevalence ratio",
    title = "B. Adjusted PRs, Oct. 2019 vs Oct. 2016"
  ) +
  scale_color_manual(
    values = c("VT" = "#CC6677", "NVT" = "#4477AA"),
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")
  ) +
  scale_shape_manual(
    values = c("Method 4" = 16, "Method 1" = 1),
    labels = c("Method 4" = "Serotype-specific density",
               "Method 1" = "Dominant serotype")
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    legend.title = element_blank()
  )

pr_sero_plot

plot2 <- plot1 / pr_sero_plot + plot_layout(guides = "collect")
plot2

ggsave("results/07_plots/figures/Fig4_Serotype prevalence and PR_dominant vs density_20261015.tiff", plot2,
       width = 9.5, height = 9, units = "in", dpi = 600, compression = "lzw")

topptx(plot2, "results/07_plots/figures/Fig4_Serotype prevalence and PR_dominant vs density_20261015.pptx",
       width = 9.5, height = 9)


########################################################################
### Figure 2: Method 2 vs Method 1 in October 2019 #####################
########################################################################

lower_cap   <- 0.01
upper_cap   <- 100
dodge_width <- 0.55
n_mod       <- 2

pr_sero_plot_data_m2 <- pr_sero_m4_adj1 |>
  filter(method %in% c("Method 1", "Method 2")) |>
  mutate(
    serotype = droplevels(serotype),
    method   = factor(method, levels = c("Method 1", "Method 2")),
    # manual dodge: Method 1 on the left, Method 2 on the right
    x_pos    = as.integer(serotype) +
      (as.integer(method) - (n_mod + 1) / 2) * dodge_width / n_mod,
    below_cap  = !is.na(q2.5_1)  & q2.5_1  < lower_cap,
    over_cap   = !is.na(q97.5_1) & q97.5_1 > upper_cap,
    lower_plot = pmin(pmax(q2.5_1,  lower_cap), upper_cap),
    upper_plot = pmin(pmax(q97.5_1, lower_cap), upper_cap),
    median_plot = ifelse(median_1 < lower_cap | median_1 > upper_cap, NA, median_1),
    arrow_end = case_when(
      below_cap & over_cap ~ "both",
      below_cap            ~ "first",
      over_cap             ~ "last",
      TRUE                 ~ "none"
    )
  )

sero_levels_m2 <- levels(pr_sero_plot_data_m2$serotype)

pr_sero_plot_m2 <- ggplot(pr_sero_plot_data_m2, aes(x = x_pos, color = sero)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             linewidth = 0.4, color = "grey30") +
  geom_segment(
    data = filter(pr_sero_plot_data_m2, arrow_end == "none"),
    aes(xend = x_pos, y = lower_plot, yend = upper_plot),
    linewidth = 0.45, alpha = 0.75, show.legend = FALSE, na.rm = TRUE
  ) +
  lapply(c("first", "last", "both"), function(e) {
    geom_segment(
      data = filter(pr_sero_plot_data_m2, arrow_end == e),
      aes(xend = x_pos, y = lower_plot, yend = upper_plot),
      linewidth = 0.45, alpha = 0.75, show.legend = FALSE,
      arrow = arrow(length = unit(1.2, "mm"), type = "closed", ends = e)
    )
  }) +
  geom_point(aes(y = median_plot, shape = method),
             size = 2.6, stroke = 0.65, na.rm = TRUE) +
  scale_x_continuous(
    breaks = seq_along(sero_levels_m2),
    labels = gsub(" ", "\n", sero_levels_m2),
    expand = expansion(add = 0.5)
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  coord_cartesian(ylim = c(lower_cap, upper_cap), clip = "off") +
  labs(
    x = "Serotype",
    y = "Prevalence ratio",
    title = "B. Adjusted PRs, Oct. 2019 vs Oct. 2016"
  ) +
  scale_color_manual(
    values = c("VT" = "#CC6677", "NVT" = "#4477AA"),
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")
  ) +
  scale_shape_manual(
    values = c("Method 2" = 16, "Method 1" = 1),
    labels = c("Method 2" = "Equal share",
               "Method 1" = "Dominant serotype")
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    legend.title = element_blank()
  )

pr_sero_plot_m2

plot2_m2 <- plot1_m2 / pr_sero_plot_m2 + plot_layout(guides = "collect")
plot2_m2

ggsave("results/07_plots/figures/FigS5_Serotype prevalence and PR_dominant vs equal share_20261015.tiff", plot2_m2,
       width = 9.5, height = 9, units = "in", dpi = 600, compression = "lzw")

topptx(plot2_m2, "results/07_plots/figures/FigS5_Serotype prevalence and PR_dominant vs equal share_20261015.pptx",
       width = 9.5, height = 9)


########################################################################
### Figure 3: Method 3 vs Method 1 in October 2019 #####################
########################################################################

lower_cap   <- 0.01
upper_cap   <- 100
dodge_width <- 0.55
n_mod       <- 2

pr_sero_plot_data_m3 <- pr_sero_m4_adj1 |>
  filter(method %in% c("Method 1", "Method 3")) |>
  mutate(
    serotype = droplevels(serotype),
    method   = factor(method, levels = c("Method 1", "Method 3")),
    # manual dodge: Method 1 on the left, Method 3 on the right
    x_pos    = as.integer(serotype) +
      (as.integer(method) - (n_mod + 1) / 2) * dodge_width / n_mod,
    below_cap  = !is.na(q2.5_1)  & q2.5_1  < lower_cap,
    over_cap   = !is.na(q97.5_1) & q97.5_1 > upper_cap,
    lower_plot = pmin(pmax(q2.5_1,  lower_cap), upper_cap),
    upper_plot = pmin(pmax(q97.5_1, lower_cap), upper_cap),
    median_plot = ifelse(median_1 < lower_cap | median_1 > upper_cap, NA, median_1),
    arrow_end = case_when(
      below_cap & over_cap ~ "both",
      below_cap            ~ "first",
      over_cap             ~ "last",
      TRUE                 ~ "none"
    )
  )

sero_levels_m3 <- levels(pr_sero_plot_data_m3$serotype)

pr_sero_plot_m3 <- ggplot(pr_sero_plot_data_m3, aes(x = x_pos, color = sero)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             linewidth = 0.4, color = "grey30") +
  geom_segment(
    data = filter(pr_sero_plot_data_m3, arrow_end == "none"),
    aes(xend = x_pos, y = lower_plot, yend = upper_plot),
    linewidth = 0.45, alpha = 0.75, show.legend = FALSE, na.rm = TRUE
  ) +
  lapply(c("first", "last", "both"), function(e) {
    geom_segment(
      data = filter(pr_sero_plot_data_m3, arrow_end == e),
      aes(xend = x_pos, y = lower_plot, yend = upper_plot),
      linewidth = 0.45, alpha = 0.75, show.legend = FALSE,
      arrow = arrow(length = unit(1.8, "mm"), type = "closed", ends = e)
    )
  }) +
  geom_point(aes(y = median_plot, shape = method),
             size = 2.6, stroke = 0.65, na.rm = TRUE) +
  scale_x_continuous(
    breaks = seq_along(sero_levels_m3),
    labels = gsub(" ", "\n", sero_levels_m3),
    expand = expansion(add = 0.5)
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  coord_cartesian(ylim = c(lower_cap, upper_cap), clip = "off") +
  labs(
    x = "Serotype",
    y = "Prevalence ratio",
    title = "B. Adjusted PRs, Oct. 2019 vs Oct. 2016"
  ) +
  scale_color_manual(
    values = c("VT" = "#CC6677", "NVT" = "#4477AA"),
    labels = c("VT" = "PCV10 serotypes", "NVT" = "Non-PCV10 serotypes")
  ) +
  scale_shape_manual(
    values = c("Method 3" = 16, "Method 1" = 1),
    labels = c("Method 3" = "Relative abundance",
               "Method 1" = "Dominant serotype")
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    legend.title = element_blank()
  )

pr_sero_plot_m3

plot2_m3 <- plot1_m3 / pr_sero_plot_m3 + plot_layout(guides = "collect")
plot2_m3

ggsave("results/07_plots/figures/FigS6_Serotype prevalence and PR_dominant vs relative abundance_20261015.tiff", plot2_m3,
       width = 9.5, height = 9, units = "in", dpi = 600, compression = "lzw")

topptx(plot2_m3, "results/07_plots/figures/FigS6_Serotype prevalence and PR_dominant vs relative abundance_20261015.pptx",
       width = 9.5, height = 9)


#
