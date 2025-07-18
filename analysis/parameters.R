# baseline_parameters.R
# This script contains all hard-coded default parameter values used in model runs.
# It is sourced by simulation and analysis scripts.

##### 0. Packages and Source Functions ----------------------------------------

library(tidyverse)
library(mpoxnetwork)


# Load and preprocess DRC population age distribution
DRC_pop_age <- DRC_pop2024 %>%
  mutate(total = male + female) %>%
  mutate(prop = total / sum(total))

# Collapse into 11 age groups, with 75+ collapsed
age_dist <- c(DRC_pop_age$prop[1:15], sum(DRC_pop_age$prop[16:21]))
age_distribution <- c(age_dist[1:10], 1 - sum(age_dist[1:10]))

# Sexually active age groups: 15-64 (index 4 to 10)
sexact_age_dist <- c(rep(0, 3), age_dist[4:13] / sum(age_dist[4:13]), rep(0, 3))

##### 1. Model Parameters -----------------------------------------------------

# Age structure
M <- 11 # number of age groups
N_1 <- 4 # index of first sexually active age group
N_2 <- 10 # index of last sexually active age group
n_active <- N_2 - N_1 + 1 # number of sexually active groups

# Degree distribution parameters
k_max <- 25 # maximum degree
N_0 <- 0.14 # proportion with zero partners
alpha <- 2.9 # Pareto exponent

degree_distribution <- compute_degree_distribution(k_max = k_max, N_0 = N_0, alpha = alpha)
degree_zero <- 1 - sum(degree_distribution)
degree_sizebiased <- compute_size_biased_distribution(degree_distribution)

# Contact matrix (age-adjusted)
contact_matrix <- prepare_contact_matrix(contact_home, age_dist, N_1, N_2, M)

# Transmission parameters
sigma <- 1 / 2 # latent period (2 days), Brosius 2023
gamma <- 1 / 8 # infectious period (8 days), Kim 2023
# Mean generation time ~10 days, Miura 2024 J Infect Dis

parms_list <- determine_transmission_parameters(age_distribution, contact_matrix, gamma)
beta <- parms_list$beta
p_cont <- parms_list$p_cont
c_ave <- parms_list$c_ave

# Mortality (case fatality rate by age)
cfr_by_age <- c(0.175, rep(0.05, 2), rep(0.025, M - 3)) # Whittles 2024 medRxiv Fig S1
matM <- diag(rep(cfr_by_age, 2)) # for 2 modes of transmission (sexual and non-sexual)

##### 2. Computational Settings -----------------------------------------------

population_size <- sum(DRC_pop_age$total)
seed_infected <- 1 / population_size # small initial infection
initial_age_group <- 5 # index for age group 20-24
time_end <- 300 # sufficient time for epidemic to end
time_increment <- 0.5 # time step for ODE solver (RK method)


####### helper variables
D_max <- length(degree_distribution)
age_distribution_sa <- age_distribution[N_1:N_2] / sum(age_distribution[N_1:N_2])
