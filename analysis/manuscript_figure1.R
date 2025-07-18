# Load required libraries
library(tibble)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(reshape2)
library(mpoxnetwork)

# source data frames
result_tibble <- readRDS(file.path(DATA_DIR, "degree_gridsearch_summary.rds"))
source("analysis/initial_ngm_dataframe.R")
source("analysis/LHS_results.R")

# labels for the age groups in plot
age_group_labels <- c(
  "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
  "30-34", "35-39", "40-44", "45-49", "50+"
)

### Figure 1ABC -----
global_NGM_agg_plot <- ggplot(df_global_NGM_agg, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "#AE2012") +
  geom_text(aes(label = round(value, 1)), color = "black", size = 3) + # Add text labels with values rounded to 2 decimal places
  labs(x = "Acquisition route", y = "Onward transmission route", fill = "Infection") +
  coord_fixed() +
  guides(fill = "none")

k_gg_plot <- ggplot(df_k_gg, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "#AE2012") +
  labs(x = "Contactor age", y = "Contactee age", fill = "Infection") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

k_ss_plot <- ggplot(df_k_ss, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "#AE2012") +
  labs(x = "Contactor age", y = "Contactee age", fill = "Infection") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

### Figure 1D -----
# Heatmap
filtered_data <- result_tibble %>%
  filter(N_0_val == N_0_val[which.min(abs(N_0_val - N_0))]) %>% # extract the closest value to the baseline N_0
  filter(k_max_val < 30) %>%
  filter(2 < alpha_val) %>%
  filter(alpha_val < 4)

heatmap_R0_plot <- ggplot(filtered_data, aes(x = alpha_val, y = k_max_val, fill = R_0)) +
  geom_tile() +
  # geom_contour(aes(z = R_0), breaks = 1, color = "black", linewidth = 1) +
  scale_fill_gradient2(low = "#005F73", mid = "white", high = "#AE2012", midpoint = 1) +
  labs(
    x = expression("Pareto exponent " * alpha),
    y = expression("Maximum cutoff " * k[max]),
    fill = expression(R[0])
  ) +
  geom_point(aes(x = alpha, y = k_max), shape = 24, color = "black", fill = "#EE9B00", size = 2.5, stroke = 0.5)

### Figure 1E -----
# LHS plot
LHS_plot <- ggplot() +
  geom_violin(
    data = LHS_results, aes(x = as.factor(age_group), y = value),
    fill = "#94D2BD", alpha = 0.5, trim = FALSE
  ) +
  geom_point(
    data = obs_data_Kami, aes(x = as.factor(age_group), y = value),
    color = "#001219", size = 1
  ) +
  labs(
    # title = "Comparison of model-projected range (blue) and observed data (orange)",
    x = "Age group",
    y = "Proportion"
  ) +
  scale_x_discrete(labels = age_group_labels[1:11]) +
  ylim(c(-0.1, 0.75)) +
  coord_flip()  # Rotate the plot 90 degrees

### Figure 1 (Merge all) ----
figure_1_plot <-(((k_gg_plot | k_ss_plot | global_NGM_agg_plot)) &
                    theme(legend.position = "bottom")) / (heatmap_R0_plot | LHS_plot) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") 

ggsave(fs::path(DATA_DIR, "figure_1_plot.png"),
  plot = figure_1_plot, 
  width = 17.8, height = 12, units = "cm", dpi = 100, bg = "white"
)
