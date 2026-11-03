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
# One results file per approach (Method 1-4) and adjustment, from 06_sensitivity/6A
dir_6A <- "results/06_sensitivity/6A"
read_6A <- function(method, adj, sheet) {
  read.xlsx(file.path(dir_6A, paste0("Sensi_6A_All_M", method, "_", adj, "_re0_20261015.xlsx")),
            sheet = sheet)
}

data_m1 <- read_6A(1, "unadj", "Preva_Class_long")
data_m2 <- read_6A(2, "unadj", "Preva_Class_long")
data_m3 <- read_6A(3, "unadj", "Preva_Class_long")
data_m4 <- read_6A(4, "unadj", "Preva_Class_long")

data_m1_adj <- read_6A(1, "adj", "Preva_Class_long")
data_m2_adj <- read_6A(2, "adj", "Preva_Class_long")
data_m3_adj <- read_6A(3, "adj", "Preva_Class_long")
data_m4_adj <- read_6A(4, "adj", "Preva_Class_long")

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
  select(year, class, approach,method, group, prevalence) %>%
  pivot_wider(names_from = year,
              values_from = prevalence) %>% 
  arrange(group, class, method)-> data_prev_table
view(data_prev_table)

write.xlsx(list(
  Prevalence_VT_NVT = data_prev,
  Prevalence_VT_NVT_Table = data_prev_table
), "results/07_plots/FigS7_6A prevalence of VT and NVT_table_20261015.xlsx")

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
    title = "A. Prevalence of PCV10+6A-type"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0), breaks = seq(0, 40, 10)) +
  scale_fill_manual(
    values = c(
      "#B2DFDB",   "#9ECAE1","#CAB2D6", "#EF8A91"
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
    title = "B. Prevalence of remaining non-PCV10 serotypes"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0), breaks = seq(0, 40, 10)) +
  scale_fill_manual(
    values = c( 
      "#B2DFDB",   "#9ECAE1","#CAB2D6", "#EF8A91"
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
data_pr_m1 <- read_6A(1, "unadj", "PR_Class_long")
data_pr_m2 <- read_6A(2, "unadj", "PR_Class_long")
data_pr_m3 <- read_6A(3, "unadj", "PR_Class_long")
data_pr_m4 <- read_6A(4, "unadj", "PR_Class_long")

##adjusted PR of VT, NVT and Overall##
data_pr_m1_adj <- read_6A(1, "adj", "PR_Class_long")
data_pr_m2_adj <- read_6A(2, "adj", "PR_Class_long")
data_pr_m3_adj <- read_6A(3, "adj", "PR_Class_long")
data_pr_m4_adj <- read_6A(4, "adj", "PR_Class_long")

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
  select(year, class, approach,method, group, PR) %>%
  pivot_wider(names_from = year,
              values_from = PR) %>% 
  arrange(group, class, method)-> data_pr_table
view(data_pr_table)


write.xlsx(list(
  PR_VT_NVT = data_pr,
  PR_VT_NVT_Table = data_pr_table
), "results/07_plots/FigS7_6A PR of VT and NVT_table_20261015.xlsx")

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
    title = "C. Adjusted PRs of PCV10+6A-type",
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
    limits = c(0.1, 1.1),
    breaks = c(0.1, 0.2, 0.5,  1),
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
    title = "D. Adjusted PRs of remaining non-PCV10 serotypes",
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
    limits = c(0.9, 8),
    breaks = c(1, 2, 4, 8),
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

ggsave("results/07_plots/figures/FigS7_6A as vaccine-associated_20261015.tiff", p_all,
       width = 12, height = 7, units = "in", dpi = 600, compression = "lzw"
)

library(eoffice)
topptx(p_all, "results/07_plots/figures/FigS7_6A as vaccine-associated_20261015.pptx",
       width = 12, height = 7)


