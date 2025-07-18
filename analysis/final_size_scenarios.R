library(tidyverse)
library(mpoxnetwork)
library(future)
library(furrr)


# load parameters for the model to the global environment
source("analysis/parameters.R")
initial_values <- initialize_infection_agegroup(M, n_active, D_max, seed_infected, degree_zero, initial_age_group, degree_distribution, degree_sizebiased)

########################################################### final size scenarios
reproduction_default <- reproduction_number(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)
r0_default <- reproduction_default$reproduction_number

# keep the reproduction number fixed, decrease beta, match p_cont
df_beta <- purrr::map(c(0.4, 0.5, 0.6, 0.7, 0.8), \(sar_modified) {
  beta_modified <- sar_modified * gamma / (1 - sar_modified)
  result <- uniroot(
    \(x) reproduction_number(beta_modified, gamma, x, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)$reproduction_number - r0_default,
    interval = c(0, 1)
  )
  c_ave <- determine_transmission_parameters(age_distribution = age_distribution, contact_matrix = contact_matrix, gamma)$c_ave
  sar_h_modified <- (c_ave / gamma) * result$root / (1 + (c_ave / gamma) * result$root)
  tibble(
    sar_s = sar_modified,
    sar_h = sar_h_modified,
    beta = beta_modified,
    p_cont = result$root,
    reproduction_number = r0_default,
    name = paste0("beta_", beta_modified, "_pcont_", p_cont, "_sars_", sar_modified, "_sarh_", sar_h_modified)
  )
}) |>
  bind_rows()

# keep the reproduction number fixed, increase p_cont, match beta
df_pcont <- purrr::map(c(0.2, 0.3, 0.4, 0.5), \(sar_modified) {
  c_ave <- determine_transmission_parameters(age_distribution = age_distribution, contact_matrix = contact_matrix, gamma)$c_ave
  p_modified <- sar_modified / (1 - sar_modified) * (gamma / c_ave)
  result <- uniroot(
    \(x) reproduction_number(x, gamma, p_modified, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)$reproduction_number - r0_default,
    interval = c(0, 1)
  )
  sar_s_modified <- result$root / (gamma + result$root)
  tibble(
    sar_s = sar_s_modified,
    sar_h = sar_modified,
    beta = result$root,
    p_cont = p_modified,
    reproduction_number = r0_default,
    name = paste0("beta_", beta, "_pcont_", p_modified, "_sars_", sar_s_modified, "_sarh_", sar_modified)
  )
}) |>
  bind_rows()

df_scenarios <- bind_rows(
  list(
    beta_r0_fixed = df_beta,
    pcont_r0_fixed = df_pcont
  ),
  .id = "id"
)

# parallel map of simulation scenarios
future::plan(multisession, workers = 8)
furrr::future_pmap(
  df_scenarios |> filter(id |> stringr::str_detect("r0_fixed")),
  \(...) {
    source("analysis/parameters.R")
    x <- tibble(...)
    beta_scenario <- x$beta
    p_scenario <- x$p_cont
    name <- x$name
    mpoxnetwork::simulate_model(
      M = M, N_1 = N_1, N_2 = N_2,
      contact_matrix = contact_matrix,
      age_distribution = age_distribution,
      degree_zero = degree_zero,
      degree_distribution = degree_distribution,
      beta = beta_scenario, p_cont = p_scenario, sigma = sigma, gamma = gamma,
      seed_infected = seed_infected,
      time_end = time_end, time_increment = time_increment,
      data_dir = DATA_DIR,
      data_name = name,
      initial = NULL
    )
  }
)

plan(sequential)
