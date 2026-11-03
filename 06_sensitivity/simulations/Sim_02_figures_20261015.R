############################################################
# 06_sensitivity/simulations/Sim_02_figures_20261015.R
#
# High-carriage simulations - step 2: Figures S9 (sim1) and S10 (sim2).
#   A, B  unadjusted prevalence of PCV10 and non-PCV10 serotype carriage,
#         Oct. 2016 and Oct. 2019, by approach
#   C, D  adjusted prevalence ratios (Oct. 2019 vs Oct. 2016) for PCV10
#         and non-PCV10 serotype carriage and each serotype, log10 scale;
#         intervals beyond the axis limits end in an arrow
#
# Input:  results/06_sensitivity/simulations/Sim_results_20261015.rds
# Output: results/06_sensitivity/simulations/figures/Fig_S9_sim1_20261015.{tiff,pptx}
#         results/06_sensitivity/simulations/figures/Fig_S10_sim2_20261015.{tiff,pptx}
############################################################

rm(list = ls())
library(tidyverse)
library(patchwork)
library(eoffice)


############################################################
#### 0) Settings ###########################################
############################################################

date_tag <- "20261015"
out_dir <- "results/06_sensitivity/simulations"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

figures <- c(sim1 = "Fig_S9", sim2 = "Fig_S10")

methods <- c("Dominant serotype", "Equal share", "Relative abundance",
             "Serotype-specific density", "Presence-based")
# Results written before 14 Oct 2026 store approach 2 as "Proportional presence"
method_rename <- c("Proportional presence" = "Equal share")

# One palette for bars and points, so the bar legend also keys panels C and D
method_cols <- set_names(c("#5AB4AC", "#4393C3", "#9970AB", "#D6604D", "#FDBF6F"), methods)

# Display order of serotypes within panels C and D (top to bottom);
# categories absent from a simulation are skipped
sero_order <- list(
  VT = c("PCV10 serotypes", "19F", "6B", "23F", "14", "4", "5", "7F", "18C"),
  NVT = c("Non-PCV10 serotypes", "15B/C", "15A", "23A", "19A", "34", "35B", "11A",
          "NT2", "6A", "Other NVT")
)

# Serotypes not shown in panel C (still included in the PCV10 serotype total)
sero_drop <- c("9V", "1")

# x-axis limits of the PR panels; intervals beyond them end in an arrow
pr_limits <- list(VT = c(0.05, 5), NVT = c(0.4, 300))

theme_1 <- theme(
  legend.position = "right",
  legend.text = element_text(size = 12),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0, size = 14, vjust = 1),
  plot.title.position = "plot",
  axis.line = element_line(colour = "black", linewidth = 0.5),
  axis.text = element_text(colour = "black", size = 10),
  axis.title = element_text(colour = "black", size = 12),
  axis.ticks = element_line(colour = "black", linewidth = 0.5)
)


############################################################
#### 1) Data ###############################################
############################################################

res <- readRDS(file.path(out_dir, paste0("Sim_results_", date_tag, ".rds")))

group_label <- c(VT = "PCV10 serotypes", NVT = "Non-PCV10 serotypes")

prev <- res$prev_class |>
  filter(run_type == "unadj", class %in% c("VT", "NVT")) |>
  mutate(method = factor(recode(method, !!!method_rename), levels = methods),
         year = factor(year, levels = c("Oct. 2016", "Oct. 2019")))
view(prev)

# Adjusted PRs: group rows plus serotype rows
pr <- bind_rows(
  res$pr_class |>
    filter(run_type == "adj", class %in% c("VT", "NVT")) |>
    mutate(group = class, row_label = group_label[class]),
  res$pr_sero |>
    filter(run_type == "adj",
           !serotype %in% sero_drop) |>
    mutate(group = sero_group, row_label = serotype)
) |>
  mutate(method = factor(recode(method, !!!method_rename), levels = methods))
view(pr)

############################################################
#### 2) Panels #############################################
############################################################

plot_prevalence <- function(d, title, y_lab) {
  ggplot(d, aes(x = year, y = median * 100, fill = method)) +
    geom_col(position = position_dodge(width = 0.5), width = 0.5,
             color = "black", linewidth = 0.3, alpha = 0.6) +
    geom_errorbar(aes(ymin = q2.5 * 100, ymax = q97.5 * 100),
                  position = position_dodge(width = 0.5), width = 0.2, linewidth = 0.3) +
    geom_vline(aes(xintercept = 1.5, color = "PCV10 introduction"),
               linetype = "dashed", linewidth = 0.4) +
    scale_fill_manual(values = method_cols, name = NULL) +
    scale_color_manual(values = c("PCV10 introduction" = "black"), name = NULL) +
    scale_y_continuous(limits = c(0, 71), breaks = seq(0, 70, 20), expand = c(0, 0)) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
    labs(title = title, x = NULL, y = y_lab) +
    theme_1
}

plot_pr_forest <- function(d, grp, title, y_lab) {
  lims <- pr_limits[[grp]]
  rows <- intersect(sero_order[[grp]], unique(d$row_label))
  rows <- c(rows, setdiff(unique(d$row_label), rows))   # any unlisted category last
  n_m <- length(methods)
  dodge <- 0.65

  d <- d |>
    filter(group == grp) |>
    mutate(
      row_pos = length(rows) + 1 - match(row_label, rows),
      y = row_pos + (n_m + 1 - 2 * as.integer(method)) / 2 * dodge / n_m,
      lo = pmin(pmax(q2.5, lims[1]), lims[2]),
      hi = pmin(pmax(q97.5, lims[1]), lims[2]),
      mid = if_else(median < lims[1] | median > lims[2], NA_real_, median),
      arrow_end = case_when(q2.5 < lims[1] & q97.5 > lims[2] ~ "both",
                            q2.5 < lims[1] ~ "first",
                            q97.5 > lims[2] ~ "last",
                            TRUE ~ "none")
    )

  ggplot(d, aes(y = y, color = method)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey30") +
    geom_segment(data = filter(d, arrow_end == "none"),
                 aes(x = lo, xend = hi, yend = y),
                 linewidth = 0.45, alpha = 0.75, show.legend = FALSE) +
    map(c("first", "last", "both"), \(e) {
      geom_segment(data = filter(d, arrow_end == e),
                   aes(x = lo, xend = hi, yend = y),
                   linewidth = 0.45, alpha = 0.75, show.legend = FALSE,
                   arrow = arrow(length = unit(1.2, "mm"), type = "closed", ends = e))
    }) +
    geom_point(aes(x = mid), size = 1.5, stroke = 0.65, na.rm = TRUE, show.legend = FALSE) +
    scale_x_log10(labels = \(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
    scale_y_continuous(breaks = seq_along(rows), labels = rev(gsub(" serotypes", "\nserotypes", rows))) +
    coord_cartesian(xlim = lims, clip = "off") +
    scale_color_manual(values = method_cols, guide = "none") +
    labs(title = title, x = "Adjusted prevalence ratio", y = y_lab) +
    theme_1
}

make_figure <- function(sim) {
  p_a <- plot_prevalence(filter(prev, simulation == sim, class == "VT"),
                         "A. Prevalence of PCV10 serotypes", "Prevalence (%)")
  p_b <- plot_prevalence(filter(prev, simulation == sim, class == "NVT"),
                         "B. Prevalence of non-PCV10 serotypes", NULL)
  p_c <- plot_pr_forest(filter(pr, simulation == sim), "VT",
                        "C. Adjusted PRs of PCV10 serotypes", "Serotype")
  p_d <- plot_pr_forest(filter(pr, simulation == sim), "NVT",
                        "D. Adjusted PRs of non-PCV10 serotypes", NULL)
  (p_a + p_b) / (p_c + p_d) + plot_layout(heights = c(0.5, 1.5), guides = "collect")
}


############################################################
#### 3) Save ###############################################
############################################################

iwalk(figures, \(fig, sim) {
  p <- make_figure(sim)
  stem <- file.path(fig_dir, paste0(fig, "_", sim, "_", date_tag))
  ggsave(paste0(stem, ".tiff"), p, width = 10, height = 10, units = "in",
         dpi = 600, compression = "lzw")
  topptx(p, paste0(stem, ".pptx"), width = 10, height = 10)
})

message("Done. Written to ", fig_dir, "/.")
