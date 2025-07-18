#' Wrapper function to simulate epidemic outbreak for the mpox network model
#'
#' @inheritParams model_mpox_network
#' @param data_name string - user-provided name for the simulation
#' @param data_dir path - user-provided data directory where simulation result is saved
#' @param x_output boolean - whether or not to save x-variables to the results
#'
#' @return tibble
#' @export
simulate_model <- function(M, N_1, N_2,
                           contact_matrix,
                           age_distribution,
                           degree_zero,
                           degree_distribution,
                           beta, p_cont, sigma, gamma,
                           seed_infected,
                           time_end, time_increment, data_dir,
                           data_name, initial, x_output = FALSE) {
  # check that results can be written
  stopifnot(fs::dir_exists(data_dir))

  reproduction <- reproduction_number(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2)
  # check that the reproduction number is larger than one
  r0 <- reproduction$reproduction_number
  message("The reproduction number is ", r0 |> round(2))
  stopifnot(r0 > 1)

  data <- model_mpox_network(
    M = M, N_1 = N_1, N_2 = N_2,
    contact_matrix = contact_matrix,
    age_distribution = age_distribution,
    degree_zero = degree_zero,
    degree_distribution = degree_distribution,
    beta = beta, p_cont = p_cont, sigma = sigma, gamma = gamma,
    seed_infected = seed_infected / population_size,
    time_end = time_end,
    time_increment = time_increment,
    initial = initial,
    x_output = x_output
  )
  data |> saveRDS(glue::glue("{data_dir}/data_{data_name}.rds"))

  return(data)
}
