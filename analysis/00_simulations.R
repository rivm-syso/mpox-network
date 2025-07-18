###############################################################################
#
# Main script to reproduce simulation results
# - the below source scripts run simulations, which take some time
# - all results are written to the `results/` directory
#
################################################################################
library(mpoxnetwork)
source("analysis/directory_paths.R")

################################################################################
# Simulated epidemics for the baseline and the community contact matrix scenarios
# - Baseline scenario with initial seed uniform over the population
# - Baseline scenario with initial seed on specific age group
# - Alternative scenario in which contact matrix is the community contact matrix

source("analysis/simulate_scenarios.R")

################################################################################
# Simulated epidemics (to compute final sizes) under different parameter settings

source("analysis/final_size_scenarios.R")

################################################################################
# Projected age-distribution of infection and reproduction numbers
# - Latin hypercube sampling to project the age-distribution of infection
# - Grid search to project the range of reproduction numbers

source("analysis/LHS_results.R")
source("analysis/range_Rvalues_simulation.R")

################################################################################
# Elasticity analysis of next-generation matrices
# - Sensitivity and elasticity of NGMs and mortality-weighted NGMs at t = 0 and 200
# - Calculate sensitivity and elasticitiy under alternative mortality matrix

source("analysis/elasticity_analysis.R")
source("analysis/supplement_sensitivity_elasticity_analysis.R")
