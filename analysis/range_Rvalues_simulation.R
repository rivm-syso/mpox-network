library(dplyr)
library(tidyr)
library(purrr)

# load parameters
source("analysis/parameters.R")

# Define ranges for the parameters
N_0_values <- seq(0.1, 0.3, length.out = 100)
alpha_values <- seq(0.5, 4, length.out = 100)
k_max_values <- seq(5, 50, length.out = 100)

# Create a grid of all parameter combinations
parameter_grid <- expand.grid(N_0 = N_0_values, alpha = alpha_values, k_max = k_max_values)

# Compute results using NGM_fn()
results_degree <- parameter_grid %>% # this takes 20-30 min
  rowwise() %>%
  mutate(
    degree_distribution = list(compute_degree_distribution(k_max, N_0, alpha)),
    ngm_result = list(tryCatch(
      NGM_fn(
        beta = beta,
        gamma = gamma,
        p_cont = p_cont,
        degree_distribution = degree_distribution,
        age_distribution = age_distribution,
        contact_matrix = contact_matrix,
        M = M,
        N_1 = N_1,
        N_2 = N_2
      ),
      error = function(e) NULL
    )),
    R_0 = if (!is.null(ngm_result)) ngm_result$reproduction_number else NA_real_,
    R_ss = if (!is.null(ngm_result)) ngm_result$R_ss else NA_real_,
    R_sg = if (!is.null(ngm_result)) ngm_result$R_sg else NA_real_
  ) %>%
  ungroup() %>%
  select(N_0, alpha, k_max, R_0, R_ss, R_sg)

# Convert to a tibble and display
result_tibble <- results_degree |>
  as_tibble() |>
  rename(
    N_0_val = N_0,
    alpha_val = alpha,
    k_max_val = k_max
  )

saveRDS(result_tibble, file = file.path(DATA_DIR, "degree_gridsearch_summary.rds"))
