#' Helper function to read RDS file containing simulation data
#'
#' @param scenario_name string - name of the scenario to read in
#' @param data_dir path - path to the directory containing the simulation data
#'
#' @return tibble
#' @export
read_simulation_data <- function(scenario_name, data_dir) {
  data <- readRDS(glue::glue("{data_dir}/data_{scenario_name}.rds"))
  return(data)
}


#' Wrangles data from simulation into incidence data ready for plotting
#'
#' @param population_size int - total population size
#' @param age_distribution vector[dbl] - age distribution in the population, sums to one
#' @param N_1 int - first sexually active age group
#' @param data tibble - simulation data of epidemic outbreak of the mpox network model
#'
#' @return tibble
#' @export
prepare_incidence_data_for_plotting <- function(data, population_size, age_distribution, N_1) {
  # incidence of global transmission
  df_incidence_global <- process_data(
    "incidence_global",
    data,
    population_size,
    age_distribution
  )

  # incidence of sexual transmission
  df_incidence_sexual <- process_data(
    "incidence_sexual",
    data,
    population_size,
    age_distribution,
    sexual_active_only = TRUE,
    N_1 = N_1
  )

  # total incidence over all age groups and modes of transmission
  df_total_incidence <- df_incidence_global |>
    bind_rows(df_incidence_sexual)

  # combine tibbles into one tibble to be able to plot incidence
  df_incidence <- bind_rows(
    list(
      total = df_total_incidence,
      non_sexual = df_incidence_global,
      sexual = df_incidence_sexual
    ),
    .id = "mode"
  )
  return(df_incidence)
}
