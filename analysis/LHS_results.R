library(lhs)
library(ggplot2)
library(dplyr)
library(mpoxnetwork)

# load parameters to the global environment
source("analysis/parameters.R")

# specify the conditions for LHS sampling
# number of samples
n_samples <- 1000

# parameter range
param_ranges <- list(
  N_0 = c(0.12, 0.2),
  alpha = c(2.9, 3.9),
  k_max = c(10L, 30L), # integer
  gamma = c(1 / 17, 1 / 7),
  SAR_H = c(0, 0.2),
  SAR_S = c(0.3, 0.95)
)

# Generate LHS samples
lhs_samples <- lhs::randomLHS(n_samples, length(param_ranges))

scaled_samples <- sapply(1:length(param_ranges), function(i) {
  low <- param_ranges[[i]][1]
  high <- param_ranges[[i]][2]
  low + (high - low) * lhs_samples[, i]
})

# Compute NGMs given sampled parameter sets
LHS_results <- run_lhs_ngm_simulation(n_samples, scaled_samples, c_ave, age_distribution, contact_matrix, M, N_1, N_2)

# Summarize by sample (repeated values per age_group → collapse to one row per sample)
summary_LHS <- LHS_results |>
  group_by(sample_id) |>
  summarize(
    reproduction_number = first(reproduction_number),
    R_gg = first(R_gg),
    R_gs = first(R_gs),
    R_sg = first(R_sg),
    R_ss = first(R_ss),
    .groups = "drop"
  )

# Compute simulated ranges (e.g. min and max)
range_stats <- summary_LHS |>
  summarize(
    R0_min = min(reproduction_number, na.rm = TRUE),
    R0_max = max(reproduction_number, na.rm = TRUE),
    R_gg_min = min(R_gg, na.rm = TRUE),
    R_gg_max = max(R_gg, na.rm = TRUE),
    R_gs_min = min(R_gs, na.rm = TRUE),
    R_gs_max = max(R_gs, na.rm = TRUE),
    R_sg_min = min(R_sg, na.rm = TRUE),
    R_sg_max = max(R_sg, na.rm = TRUE),
    R_ss_min = min(R_ss, na.rm = TRUE),
    R_ss_max = max(R_ss, na.rm = TRUE)
  )

# Round and print the ranges
round(range_stats, 2)
