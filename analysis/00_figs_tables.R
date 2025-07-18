################################################################################
#
# Main script to reproduce figures and tables
# - generates all figures and tables used in the manuscript by reading pre-computed simulation outputs from the `results/` directory
# - figures and tables are saved to the results directory
#
################################################################################
library(mpoxnetwork)
library(ggplot2)

source("analysis/directory_paths.R")

# control ggplot font sizes
theme_set(
  theme_minimal() +
  theme(
    strip.text = element_text(size = 7),
    text = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  ))
################################################################################
# Figure 1: projected next-generation matrices and age-distribution of infection
# - Projected NGMs (panel (A)(B)(C))
# - Projected range of R0 (panel (D))
# - Projected age-distribution of infection (panel (E))

source("analysis/manuscript_figure1.R")

################################################################################
# Figure 2: projected epidemic trajectories and elasticity measures
# - Projected epidemics, time-varying R, mean degree (panel (A)(B))
# - Age-specific contribution to infection/death at t = 0 (panel (C)(D))
# - Age-specific contribution to infection/death at t = 200 (panel (E)(F))

source("analysis/manuscript_figure2.R")

################################################################################
# Table 1: final sizes by scenario

source("analysis/manuscript_table1.R")

################################################################################
# Supplementary Figures and Tables
# - Comparison between the baseline and the community contact matrix scenarios
# - Projected next-generation matrix and reproduction numbers
# - Elasticity and sensitivity over time
# - Contact matrices: household and community
# - Incidence by mode of transmission for scenarios of the final size table
# - Elasticity and sensitivity under alternative mortality matrix

source("analysis/supplemental_figure.R")
source("analysis/supp_NGM_figure.R")
source("analysis/elasticity_time_figure.R")
source("analysis/supplement_contactmatrix.R")
source("analysis/supplement_figure_simulations.R")
source("analysis/supplement_figure_sensitivity_elasticity_analysis.R")
