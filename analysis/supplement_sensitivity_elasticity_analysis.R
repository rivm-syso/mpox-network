library(mpoxnetwork)


### load parameters for the model to the global environment
source("analysis/parameters.R")
seed_infected <- 0
initial_values <- initialize_infection_agegroup(M, n_active, D_max, seed_infected, degree_zero, initial_age_group, degree_distribution, degree_sizebiased)
# extract from vector and cast into arrays with correct dimensions
S_initial <- initial_values[grepl("Susceptible", names(initial_values))]
xS_initial <- initial_values[grepl("xS", names(initial_values))]
xS_initial <- array(xS_initial, dim = c(n_active, n_active, D_max, D_max))
xI_initial <- initial_values[grepl("xI", names(initial_values))]
xI_initial <- array(xI_initial, dim = c(n_active, n_active, D_max))
xbar_initial <- initial_values[grepl("xbar", names(initial_values))]
xbar_initial <- array(xbar_initial, dim = c(n_active, n_active, D_max))

seed_infected <- 1 / population_size # make sure to re-define the baseline seed_infected

#### Mortality
# https://worldhealthorg.shinyapps.io/mpx_global/#severity
cfr_by_age_alternative <- c(0.002, rep(0.001, 2), rep(0.002, M - 3)) 
matM_alternative <- diag(rep(cfr_by_age_alternative, 2)) # for 2 modes of transmission (sexual and non-sexual)


### NGM and elasticity at t = 0
ngm_0 <- time_varying_ngm(
  beta, p_cont, gamma, contact_matrix, age_distribution, age_distribution_sa,
  n_active, degree_distribution, degree_sizebiased, D_max, M, N_1, N_2,
  S_initial, xS_initial, xI_initial, xbar_initial
)

S_in_0 <- matrix(1, nrow = n_active, ncol = D_max)
S_in_powered_0 <- compute_S_in_powered(S_in_0)
ngm_vac_0 <- compute_vaccinated_ngm(S_in_powered_0,
  list(S = S_initial, xS = xS_initial, xI = xI_initial, xbar = xbar_initial),
  target = "random", deltaU = 1,
  gamma = gamma, beta = beta
)
ela_0 <- compute_cumulative_elasticity(ngm_0$global_NGM_large, ngm_vac_0$global_NGM_large, M)
sens_0 <- compute_cumulative_sensitivity(ngm_0$global_NGM_large, ngm_vac_0$global_NGM_large, M)
ela_0_matM_alternative <- compute_cumulative_elasticity(matM_alternative %*% ngm_0$global_NGM_large, matM_alternative %*% ngm_vac_0$global_NGM_large, M)
sens_0_matM_alternative <- compute_cumulative_sensitivity(matM_alternative %*% ngm_0$global_NGM_large, matM_alternative %*% ngm_vac_0$global_NGM_large, M)

# At t = 200
sim_output <- readRDS(fs::path(DATA_DIR, "data_elasticity_analysis.rds"))
df_inf_elas_overtime <- compute_time_series_ela_or_sens(sim_output,
  max_t = 200, outcome = "infection", method = "elasticity",
  gamma = gamma, beta = beta, matM = matM_alternative
)
df_inf_sens_overtime <- compute_time_series_ela_or_sens(sim_output,
  max_t = 200, outcome = "infection", method = "sensitivity",
  gamma = gamma, beta = beta, matM = matM_alternative
)
df_mort_elas_overtime <- compute_time_series_ela_or_sens(sim_output,
  max_t = 200, outcome = "death", method = "elasticity",
  gamma = gamma, beta = beta, matM = matM_alternative
)
df_mort_sens_overtime <- compute_time_series_ela_or_sens(sim_output,
  max_t = 200, outcome = "death", method = "sensitivity",
  gamma = gamma, beta = beta, matM = matM_alternative
)

# save results as RDS
result_pertub_overtime_alt <- list(
  inf_elas = df_inf_elas_overtime,
  inf_sens = df_inf_sens_overtime,
  mort_elas = df_mort_elas_overtime,
  mort_sens = df_mort_sens_overtime
)

saveRDS(result_pertub_overtime_alt, fs::path(DATA_DIR, "sensitivity_result_perturb_overtime.rds"))

# save results at t = 0 and 200 as RDS (long data.frame format)
result_pertub_snapshots_long_alt <- format_barplot_data(
  ela_0, sens_0,
  ela_0_matM_alternative, sens_0_matM_alternative,
  df_inf_elas_overtime, df_inf_sens_overtime,
  df_mort_elas_overtime, df_mort_sens_overtime,
  time_obs = 200,
  M,
  age_labels = c(
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
    "30-34", "35-39", "40-44", "45-49", "50+"
  )
)

saveRDS(result_pertub_snapshots_long_alt, fs::path(DATA_DIR, "sensitivity_result_pertub_snapshots_long.rds"))

