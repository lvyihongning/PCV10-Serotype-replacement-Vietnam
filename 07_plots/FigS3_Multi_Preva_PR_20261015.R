rm(list = ls())
library(tidyverse)
library(openxlsx)
library(patchwork)

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
# All four approaches (Method 1-4), unadjusted and adjusted, are stacked in one
# workbook; each sheet already carries `method` and `group` (unadjusted/adjusted).
res_file <- "results/03_multiple_composition/Multi_M1to4_results_re0_20261015.xlsx"

#### 1. Unadjusted Prevalence ####
#### 1.1 Prevalence of VT and NVT ####
##unadjusted and adjusted VT and NVT##
read.xlsx(res_file, sheet = "Share_class") %>%
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
  )-> data_prev

#data_prev$approach <- factor(data_prev$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
data_prev$method <- factor(data_prev$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))
data_prev$year <- factor(data_prev$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

data_prev %>% 
  select(year, class, approach,method,  group, prevalence) %>%
  pivot_wider(names_from = year,
              values_from = prevalence) %>% 
  arrange(group, class, method)-> data_prev_table

write.xlsx(list(
  Prevalence_VT_NVT = data_prev,
  Prevalence_VT_NVT_Table = data_prev_table
), "results/07_plots/FigS3_Multi proportion of VT and NVT_table_20261015.xlsx")

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
    y = "Proportion (%)",
    title = "A. Proportion of PCV10 serotypes"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, 25)) +
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
    title = "B. Proportion of non-PCV10 serotypes"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, 25)) +
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
##unadjusted and adjusted PR of VT and NVT##
read.xlsx(res_file, sheet = "PR_class") %>%
  mutate(
    median_1 = round(median, 2),
    q2.5_1   = round(q2.5, 2),
    q97.5_1  = round(q97.5, 2),
    PR = paste0(median_1, " (", q2.5_1, "-", q97.5_1, ")"),
    approach = case_when(
      method == "Method 4" ~ "Approach 1",
      method == "Method 1" ~ "Approach 2",
      method == "Method 2" ~ "Approach 3",
      method == "Method 3" ~ "Approach 4"
    ) ) %>% 
  mutate(
    q2.5_1 = ifelse(q2.5_1 < 0.01, round(q2.5, 3), q2.5_1),
    year = ifelse(year == "Oct. 2017", "Nov. 2017", year)
  )-> data_pr

#data_pr$approach <- factor(data_pr$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
data_pr$method <- factor(data_pr$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))
data_pr$year <- factor(data_pr$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

data_pr %>% 
  select(year, class, approach, method, group, PR) %>%
  pivot_wider(names_from = year,
              values_from = PR) %>% 
  arrange(group, class, method)-> data_pr_table
view(data_pr_table)

write.xlsx(list(
  PR_VT_NVT       = data_pr,
  PR_VT_NVT_Table = data_pr_table
), "results/07_plots/FigS3_Multi PR of VT and NVT_table_20261015.xlsx")


###Plot of aPR of VT and NVT###
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
    title = "C. Adjusted proportion ratios of PCV10 serotypes",
    x = "Survey time",
    y = "Proportion ratio"
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
    limits = c(0.001, 1.2),
    breaks = c(0.001, 0.1, 0.5,  1),
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
    title = "D. Adjusted proportion ratios of non-PCV10 serotypes",
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
    limits = c(0.8, 2.6),
    breaks = c(0.8, 1, 1.5, 2, 2.5),
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

ggsave("results/07_plots/figures/FigS3_Multiple carriage_VT and NVT proportion and PR_20261015.tiff", p_all,
       width = 12, height = 7, units = "in", dpi = 600, compression = "lzw"
)

library(eoffice)
topptx(p_all, "results/07_plots/figures/FigS3_Multiple carriage_VT and NVT proportion and PR_20261015.pptx",
       width = 12, height = 7)

### 1.2 Serotype prevalence ####
##unadjusted and adjusted serotype prevalence##
read.xlsx(res_file, sheet = "Share_sero") %>%
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
         year = ifelse(year == "Oct. 2017", "Nov. 2017", year)
         )-> data_sero_prev

data_sero_prev$approach <- factor(data_sero_prev$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
#data_sero_prev$method <- factor(data_sero_prev$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))
data_sero_prev$year <- factor(data_sero_prev$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

str(data_sero_prev)

data_sero_prev %>%
  select(year, serotype, sero, approach, group, prevalence) %>%
  pivot_wider(names_from = year,
              values_from = prevalence) %>%
  arrange(group, sero, serotype, approach) -> data_sero_prev_table

write.xlsx(list(
  Serotype_Prevalence       = data_sero_prev,
  Serotype_Prevalence_Table = data_sero_prev_table
), "results/07_plots/NotInManuscript_Multi serotype proportion_table_20261015.xlsx")

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

write.xlsx(data_prev_sero_2016, "results/07_plots/TableS5_Multi proportion in October 2016_20261015.xlsx")


##### Plot unadjusted serotype prevalence ####
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
  # Limits set in coord_cartesian() below, so CIs reaching under 0.1% are
  # clipped at the axis instead of being dropped
  scale_x_log10(
    breaks = c(0.1, 1, 5, 20, 50), # Adjust breaks to log-friendly numbers
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  scale_y_log10(
    breaks = c(0.1, 1, 5,  20, 50), # Adjust breaks to log-friendly numbers
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  coord_cartesian(xlim = c(0.1, 50), ylim = c(0.1, 50)) +
  labs(
    title = "A. Proportion of individual serotypes among co-colonisation",
    x = "Proportion in Oct. 2016 (%)",
    y = "Proportion in Oct. 2019 (%)"
  ) +
  theme_1 +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14),
    legend.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical", 
    axis.text = element_text(color = "black") )

plot1

#### 2. Adjsuted PR ####
### 2.2 Adjusted PR of individual serotypes ####
##unadjusted and adjusted PR of individual serotypes##
read.xlsx(res_file, sheet = "PR_sero") %>%
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
    ),
    year = ifelse(year == "Oct. 2017", "Nov. 2017", year)
  )-> pr_sero

pr_sero$approach <- factor(pr_sero$approach, levels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))
#pr_sero$method <- factor(pr_sero$method, levels = c("Method 1", "Method 2", "Method 3", "Method 4"))
pr_sero$year <- factor(pr_sero$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

pr_sero %>%
  select(year, serotype, sero, approach, group, PR) %>%
  pivot_wider(names_from = year,
              values_from = PR) %>%
  arrange(group, sero, serotype, approach) -> pr_sero_table

write.xlsx(list(
  Serotype_PR       = pr_sero,
  Serotype_PR_Table = pr_sero_table
), "results/07_plots/NotInManuscript_Multi serotype PR_table_20261015.xlsx")

##### Plot adjusted PR of individual serotypes ####
###only Oct. 2019 and include method 1####
pr_sero %>%
  filter(sero != "NO") %>% 
  filter(group == "adjusted") %>% 
  filter(method %in% c("Method 4", "Method 1")) %>% 
  filter(year == "Oct. 2019") %>% 
  select(sero, serotype, year, method, approach, group, PR, median, q2.5, q97.5) -> pr_sero_m4_adj1
#filter(serotype %in% c("19F", "6B", "23F", "15B/C", "6A", "15A", "23A", "NT2"))    
#filter(!(year %in% c("Oct. 2020", "July 2022"))) 

pr_sero_m4_adj1$sero <- factor(pr_sero_m4_adj1$sero, levels =  c("VT", "NVT"))

pr_sero_m4_adj1$serotype <- factor(pr_sero_m4_adj1$serotype,
                                   levels = c(
                                     "19F", "6B", "23F", "14", "9V",  "1",
                                     "15B/C", "6A",  "15A", "19A", "11A", "35B", "23A","NT2",  "34", "13", "3", "6C", "Other NVT"))

lower_cap   <- 0.01
upper_cap   <- 100
dodge_width <- 0.55
n_mod       <- 2

pr_sero_plot_data <- pr_sero_m4_adj1 |>
  mutate(
    serotype = droplevels(serotype),
    method   = factor(method, levels = c("Method 1", "Method 4")),
    # manual dodge: Method 1 on the left, Method 4 on the right
    x_pos    = as.integer(serotype) +
      (as.integer(method) - (n_mod + 1) / 2) * dodge_width / n_mod,
    # unrounded quantiles, so e.g. 0.006 is not rounded up to the 0.01 cap
    below_cap  = !is.na(q2.5)  & q2.5  < lower_cap,
    over_cap   = !is.na(q97.5) & q97.5 > upper_cap,
    lower_plot = pmin(pmax(q2.5,  lower_cap), upper_cap),
    upper_plot = pmin(pmax(q97.5, lower_cap), upper_cap),
    median_plot = ifelse(median < lower_cap | median > upper_cap, NA, median),
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
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
  ) +
  coord_cartesian(ylim = c(lower_cap, upper_cap), clip = "off") +
  labs(
    x = "Serotype",
    y = "Proportion ratio",
    title = "B. Adjusted proportion ratios, Oct. 2019 vs Oct. 2016"
  ) +
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
  theme_1 + 
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 14, vjust = 1),
    legend.title = element_blank() # This will keep hiding the title for both
  )
pr_sero_plot


library(patchwork)
plot2 <- plot1 / pr_sero_plot + plot_layout(guides = "collect")
plot2

ggsave("results/07_plots/figures/NotInManuscript_Multiple carriage_serotype proportion and PR_20261015.tiff", plot2,
       width = 9.5, height = 9, units = "in", dpi = 600, compression = "lzw"
)

library(eoffice)
topptx(plot2, "results/07_plots/figures/NotInManuscript_Multiple carriage_serotype proportion and PR_20261015.pptx",
       width = 9.5, height = 9)




