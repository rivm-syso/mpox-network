# Load required libraries
library(tibble)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(mpoxnetwork)

# source data frames
result_pertub_overtime <- readRDS(file.path(DATA_DIR, "result_perturb_overtime.rds"))

# heatmap over time by outcome and sensitivity measure (elasticity or sensitivity)
age_labels <- c(
  "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50+"
)

plot_inf_elas_overtime <- result_pertub_overtime$inf_elas |>
  ggplot(
    aes(x = time, y = factor(age_group), fill = value)
  ) +
  geom_tile() +
  scale_fill_gradient2(
    # low = "darkblue",
    mid = "white",
    high = "#BB3E03",
    # midpoint = 0,
    name = "Elasticity"
  ) +
  scale_y_discrete(labels = age_labels) +
  labs(
    title = "Elasticity (Infection)",
    x = "Time (days)",
    y = "Age group"
  ) 

plot_inf_sens_overtime <- result_pertub_overtime$inf_sens |>
  ggplot(
    aes(x = time, y = factor(age_group), fill = value)
  ) +
  geom_tile() +
  scale_fill_gradient2(
    # low = "darkblue",
    mid = "white",
    high = "#EE9B00",
    # midpoint = 0,
    name = "Sensitivity"
  ) +
  scale_y_discrete(labels = age_labels) +
  labs(
    title = "Sensitivity (Infection)",
    x = "Time (days)",
    y = "Age group"
  )  

plot_mort_elas_overtime <- result_pertub_overtime$mort_elas |>
  ggplot(
    aes(x = time, y = factor(age_group), fill = value)
  ) +
  geom_tile() +
  scale_fill_gradient2(
    # low = "darkblue",
    mid = "white",
    high = "#005F73",
    # midpoint = 0,
    name = "Elasticity"
  ) +
  scale_y_discrete(labels = age_labels) +
  labs(
    title = "Elasticity (Mortality)",
    x = "Time (days)",
    y = "Age group"
  ) 

plot_mort_sens_overtime <- result_pertub_overtime$mort_sens |>
  ggplot(
    aes(x = time, y = factor(age_group), fill = value)
  ) +
  geom_tile() +
  scale_fill_gradient2(
    # low = "darkblue",
    mid = "white",
    high = "#94D2BD",
    # midpoint = 0,
    name = "Sensitivity"
  ) +
  scale_y_discrete(labels = age_labels) +
  labs(
    title = "Sensitivity (Mortality)",
    x = "Time (days)",
    y = "Age group"
  ) 

figure_S2 <- (plot_inf_elas_overtime / plot_mort_elas_overtime) | (plot_inf_sens_overtime / plot_mort_sens_overtime)
figure_S2 <- figure_S2 + plot_annotation(tag_levels = "A")

ggsave(file.path(DATA_DIR, "figure_S2_elasticity_time.png"),
  plot = figure_S2, width = 17.8, height = 11, units = "cm", dpi = 100, bg = "white"
)
