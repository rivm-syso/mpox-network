library(tidyverse)
library(mpoxnetwork)


# load parameters for the model to the global environment
source("analysis/parameters.R")


#################################################################
# baseline scenario with initial seed on specific age group
initial_values <- initialize_infection_agegroup(M, n_active, D_max, seed_infected, degree_zero, initial_age_group, degree_distribution, degree_sizebiased)
baseline <- simulate_model(
  M = M, N_1 = N_1, N_2 = N_2,
  contact_matrix = contact_matrix,
  age_distribution = age_distribution,
  degree_zero = degree_zero,
  degree_distribution = degree_distribution,
  beta = beta, p_cont = p_cont, sigma = sigma, gamma = gamma,
  seed_infected = seed_infected,
  time_end = time_end, time_increment = time_increment,
  initial = initial_values, data_dir = DATA_DIR,
  data_name = "baseline"
)


#################################################################
# rerun simulation but now with x-variables
# needed for elasticity analysis up until time 200
initial_values <- initialize_infection_agegroup(M, n_active, D_max, seed_infected, degree_zero, initial_age_group, degree_distribution, degree_sizebiased)
elasticity <- simulate_model(
  M = M, N_1 = N_1, N_2 = N_2,
  contact_matrix = contact_matrix,
  age_distribution = age_distribution,
  degree_zero = degree_zero,
  degree_distribution = degree_distribution,
  beta = beta, p_cont = p_cont, sigma = sigma, gamma = gamma,
  seed_infected = seed_infected,
  time_end = 200, time_increment = time_increment,
  initial = initial_values, data_dir = DATA_DIR, x_output = TRUE,
  data_name = "elasticity_analysis"
)


#################################################################
# contact matrix is the community contact matrix
community_contact_matrix <- prepare_contact_matrix(contact_all, age_dist, N_1, N_2, M)
parms_list <- determine_transmission_parameters(age_distribution, community_contact_matrix, gamma)
beta <- parms_list$beta
p_cont <- parms_list$p_cont
initial_values <- initialize_infection_agegroup(M, n_active, D_max, seed_infected, degree_zero, initial_age_group, degree_distribution, degree_sizebiased)
community <- simulate_model(
  M = M, N_1 = N_1, N_2 = N_2,
  contact_matrix = community_contact_matrix,
  age_distribution = age_distribution,
  degree_zero = degree_zero,
  degree_distribution = degree_distribution,
  beta = beta, p_cont = p_cont, sigma = sigma, gamma = gamma,
  seed_infected = seed_infected,
  time_end = time_end, time_increment = 0.5,
  initial = initial_values, data_dir = DATA_DIR,
  data_name = "community_contact_matrix"
)
