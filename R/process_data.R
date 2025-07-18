#' Process simulation data to population numbers in long format for further analyses
#'
#' @inheritParams model_mpox_network
#' @param data tibble - epidemic outbreak data, output of `model_mpox_network`
#' @param population_size int - size of the total population
#' @param var_str string - name of variable under consideration, e.g. S, E, I, R, xS, etc.
#' @param sexual_active_only boolean - whether or variables are only applicable for sexually active age groups, default is FALSE
#'
#' @return tibble - subset of data, modified to population numbers, in long format
#' @export
process_data <- function(var_str, data, population_size, age_distribution, sexual_active_only = FALSE, N_1 = NULL) {
  data <- data |>
    select(time, starts_with(var_str)) |>
    tidyr::pivot_longer(cols = !starts_with("t")) |>
    mutate(age_group = stringr::str_extract(name, "(\\d)+$") |> as.numeric()) |>
    mutate(fraction_age_group = value |> as.numeric())

  # if there is only indexing of sexual active age groups, index 1 corresponds to first sexually active age group
  # add first sexually active age group to index to obtain the correct age group
  if (sexual_active_only) {
    data <- data |>
      mutate(age_group = age_group + N_1 - 1)
  }

  # convert fractions to population numbers with the correct indexing
  data <- data |>
    mutate(number = population_size * age_distribution[age_group] * value)

  return(data)
}
