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

##load data
preva <- read.xlsx("results/01_carriage_multiple/0_co-colonisation_prevalence_PR_unadj_20261015.xlsx" , sheet = "preva_colonistation")

prop_multi <- read.xlsx("results/01_carriage_multiple/0_co-colonisation_proportion_PR_unadj_20261015.xlsx" , sheet = "prev_class")

### 0.1 Overall: prevalence of NO, single, two and three serotypes###
tmp <- preva %>%
  filter(category %in% c( "m1", "m2", "m3")) %>%
  select( year, category, median) %>%
  mutate(median = round(median * 100, 1),
         year = ifelse(year == "Oct. 2017", "Nov. 2017", year)) 
sum(tmp$median)

tmp$category<- factor(tmp$category,
                      levels = c( "m1", "m2", "m3"),
                      labels = c(
                       "Carrying one serotype",
                        "Carrying two serotypes", "Carrying three serotypes"
                      ))
tmp$year <- factor(tmp$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

### stacked bar plot###
# creat a variable to control the order of the stack
tmp$stack_order <- factor(tmp$category, levels = rev(levels(tmp$category)))

tmp <- tmp %>%
  group_by(year) %>%
  arrange(year, desc(as.numeric(stack_order))) %>%
  mutate(label_position = cumsum(median) - 0.5 * median) %>%
  mutate(label_position = ifelse(stack_order == "Carrying three serotypes", label_position +1, label_position)) %>%
  ungroup()
view(tmp)

vline_data <- data.frame(
  xintercept = c(4.7 / 3, 5.4), # Your original xintercept values
  line_name = c("Routine vaccination  begins", "COVID-19 begins"), # Names for your legend
  line_color = c("black", "red") # Explicit colors for each line
)
vline_data$line_name <- factor(vline_data$line_name, levels = c("Routine vaccination  begins", "COVID-19 begins"))

cols_A <- c(
  "One serotype"   = "#A6CEE3",
  "Two serotypes" = "#CAB2D6",
  "Three serotypes" = "#FDBF6F"
)

cols_B <- c(
  "PCV10-only" = "#F4A3A8",
  "Mixed, PCV10-dominant" = "#FDBB84",
  "Mixed, non-PCV10-dominant" = "#B2DFDB",
  "Non-PCV10-only" = "#9ECAE1"
)

t_plot_stack <- ggplot(tmp, aes(x = year, y = median, fill = category)) +
  annotate("rect",
           xmin = 5.4, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3
  ) +
  geom_bar(
    stat = "identity", position = "stack", width = 0.5,
    aes(group = stack_order)
  ) +
  geom_text(
    data = ~ filter(.x, median > 0),   # no labels for 0.0 values
    aes(y = label_position, label = sprintf("%.1f", median)),
    size = 3
  )+
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 300, by = 10),
                     expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "#A6CEE3", "#CAB2D6", "#FDBF6F",
      "#b3cde3", "#decbe4","#fdcdac", "#fbb4ae", 
      "#cbd5e8", "#f4cae4", "#fdcdac", "#f37974",
      "#9CD5E6",  "#D0C5DE","#fdcdac","#f37974", "#fbb4b9" ,
      "#ccebc5", "#fbb4ae", "#b3cde3", "#decbe4"),
    labels = c("Carrying one serotype" = "One serotype",
               "Carrying two serotypes" = "Two serotypes",
               "Carrying three serotypes" = "Three serotypes"),
    name = NULL
  ) +
  labs(
    title = "A. Prevalence",    # 1. Clean title without the letter
    x = "Survey time",
    y = "Prevalence (%)"
  ) +
  theme_1 +
  theme(
    # --- Tag and Title Placement ---
    plot.title.position = "plot",                      # Pushes the title to the left edge of the image
    plot.title = element_text(size = 14, hjust = 0, vjust = 0),
    # --- Legend Sizing ---
    legend.position = "top",
    legend.text = element_text(size = 8.5, face = "bold"),             # 2. Controls the text size in the legend
    legend.key.size = unit(0.3, "cm")                  # 2. Controls the size of the colored boxes
  ) + 
  geom_vline(aes(xintercept = 4.7 / 3),
             linetype = "dashed",
             linewidth = 0.3
  ) 
# +guides(fill = guide_legend(nrow = 1))
t_plot_stack

####proprtion of VT and NVT within co-colonised samples####
prop_multi %>% 
  select(year, category, prevalence, median) %>% 
  mutate(
    prop = round(median *100, 1),
    year = ifelse(year == "Oct. 2017", "Nov. 2017", year)
  ) -> prop_plot

prop_plot$category <- factor(prop_plot$category,
                             levels = c("VT-only", "Mixed, VT-dominant", "Mixed, NVT-dominant", "NVT-only"  ),
                             labels = c("PCV10-only", "Mixed PCV10-dominant", "Mixed non-PCV10-dominant", "Non-PCV10-only" ))
prop_plot$year <- factor(prop_plot$year, levels =  c("Oct. 2016", "June 2017", "Nov. 2017", "Oct. 2018", "Oct. 2019", "Oct. 2020", "July 2022"))

#prop_plot$stack_order <- factor(prop_plot$category, levels = rev(levels(prop_plot$category)))
prop_plot <- prop_plot %>%
  group_by(year) %>%
  arrange(year, desc(as.numeric(category))) %>%
  mutate(label_position = cumsum(prop) - 0.5 * prop) %>%
  mutate(label_position = ifelse(category == "PCV10-only" & year %in% c("Oct. 2018", "Oct. 2019", "Oct. 2020"),
                                 label_position + 2, label_position)) %>%
  mutate(label_position = ifelse(year == "July 2022" & category %in% c("PCV10-only", 
                                                                        "Mixed PCV10-dominant", "Mixed non-PCV10-dominant"), 
                                 NA, label_position)) %>%
  ungroup()


p_co <- ggplot(prop_plot, aes(x = year, y = prop, fill = category)) +
  annotate("rect",
           xmin = 5.4, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3
  ) +
  geom_bar(
    stat = "identity", position = "stack", width = 0.5
  ) +
  geom_text(
    data = ~ filter(.x, prop > 0),     # no labels for 0.0 values
    aes(y = label_position, label = sprintf("%.1f", prop)),
    size = 3
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, 103.5),
                     expand = c(0, 0)) +
  scale_fill_manual(
    values = c( 
      "#F4A3A8", "#FDBB84", "#B2DFDB",   "#9ECAE1",
      "#fbb4ae","#fdcdac","#decbe4", "#b3cde3","#9CD5E6",
               "#fbb4ae","#decbe4", "#b3cde3",  "#ccebc5"),
    name = NULL
  ) +
  labs(
    title = "B. Serotype composition among multiple serotype carriage",    # 1. Clean title without the letter
    x = "Survey time",
    y = "Proportion (%)"
  ) +
  theme_1 +
  theme(
    # --- Tag and Title Placement ---
    plot.title.position = "plot",                      # Pushes the title to the left edge of the image
    plot.title = element_text(size = 14, hjust = 0, vjust = 0),
    
    # --- Legend Sizing ---
    legend.position = "top",
    legend.text = element_text(size = 8, face = "bold"),             # 2. Controls the text size in the legend
    legend.key.size = unit(0.25, "cm")                  # 2. Controls the size of the colored boxes
  ) + 
  geom_vline(aes(xintercept = 4.7 / 3),
             linetype = "dashed",
             linewidth = 0.3
  ) 
#+  guides(fill = guide_legend(nrow = 1))
p_co

library(patchwork)
final_plot <- t_plot_stack + p_co 
final_plot

#
ggsave("results/07_plots/figures/Fig1_Carriage states and types of multiple carriage_20261015.tiff", final_plot,
       width = 11.5, height = 6, units = "in", dpi = 600, compression = "lzw")

library(eoffice)
topptx(final_plot, file = "results/07_plots/figures/Fig1_Carriage states and types of multiple carriage_20261015.pptx",
       width = 11.5, height = 6, units = "in")


