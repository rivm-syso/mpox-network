#' Compute size-biased degree distribution
#'
#' @param degree_distribution A numeric vector representing the degree distribution (excluding zero degree)
#'
#' @return A normalized vector representing the size-biased degree distribution
#' @export
compute_size_biased_distribution <- function(degree_distribution) {
  D_max <- length(degree_distribution)
  degree_sizebiased <- degree_distribution * (1:D_max)
  degree_sizebiased / sum(degree_sizebiased)
}

#' Compute S_in matrix from xbar and age distribution
#'
#' @param xbar A 3D array [i, j, n] representing the susceptible binding sites
#' @param age_distribution_sa A numeric vector for sexually active age distribution
#'
#' @return A matrix [i, n] where each element is the age-weighted sum over j
#' @export
compute_S_in <- function(xbar, age_distribution_sa) {
  apply(xbar, c(1, 3), function(x) sum(x * age_distribution_sa))
}

#' Raise S_in matrix elements to power n (column index)
#'
#' @param S_in A matrix [i, n] representing weighted susceptible binding sites
#'
#' @return A matrix [i, n] with each element raised to the power of its column index
#' @export
compute_S_in_powered <- function(S_in) {
  D_max <- ncol(S_in)
  t(apply(S_in, 1, function(row) row^(1:D_max)))
}

#' Extract reshaped x variables from model output at a specific time point
#'
#' @param model_output A data frame output from model_mpox_network
#' @param time_point A numeric time at which to extract the state
#' @param n_active Number of sexually active age groups
#' @param D_max Maximum degree
#'
#' @return A list containing reshaped x variables (x0, xbar, xE, xI, xS, S)
#' @export
extract_x_arrays_at_time <- function(model_output, time_point, n_active, D_max) {
  time_row <- which.min(abs(model_output$time - time_point))

  extract_array <- function(name, dims) {
    cols <- grep(paste0("^", name), names(model_output), value = TRUE)
    values <- as.numeric(model_output[time_row, cols])
    array(values, dim = dims)
  }

  x0 <- as.numeric(model_output[time_row, grep("^xsingle", names(model_output))])
  xbar <- extract_array("xbar", c(n_active, n_active, D_max))
  xE <- extract_array("xE", c(n_active, n_active, D_max))
  xI <- extract_array("xI", c(n_active, n_active, D_max))
  xS <- extract_array("xS", c(n_active, n_active, D_max, D_max))
  S <- as.numeric(model_output[time_row, grep("^Susceptible", names(model_output))])

  return(list(
    time = model_output$time[time_row],
    S = S,
    x0 = x0,
    xbar = xbar,
    xE = xE,
    xI = xI,
    xS = xS
  ))
}

#' Compute the next-generation matrix (NGM) at a given time
#'
#' @param model_output A data frame output from model_mpox_network
#' @param time_point Time at which to compute the NGM
#' @param age_distribution_sa Age distribution among sexually active groups
#' @param n_active Number of sexually active age groups
#' @param D_max Maximum degree
#' @param ... Additional parameters passed to time_varying_ngm
#'
#' @return A list with `ngm` and the extracted `x_state`
#' @export
compute_ngm_at_time <- function(model_output, time_point, age_distribution_sa,
                                n_active, D_max, ...) {
  x_state <- extract_x_arrays_at_time(model_output, time_point, n_active, D_max)
  ngm <- time_varying_ngm(
    beta, p_cont, gamma,
    contact_matrix,
    age_distribution, age_distribution_sa,
    n_active, degree_distribution, degree_sizebiased, D_max,
    M, N_1, N_2,
    x_state$S, x_state$xS, x_state$xI, x_state$xbar
  )
  list(ngm = ngm, x_state = x_state)
}

#' Compute vaccinated NGM given a S_in matrix and state
#'
#' @param S_in_powered A matrix [i, n] of S_in raised to degree powers
#' @param x_state A list containing S, xS, xI, xbar
#' @param target default is "random"
#' @param deltaU default is 1
#'
#' @return A named list returned by time_varying_ngm
#' @export
compute_vaccinated_ngm <- function(S_in_powered, x_state, target = "random", deltaU = 1,
                                   gamma = gamma, beta = beta) {
  effV <- eff_vaccov_matrix(
    deltaU = deltaU, VE = 0.8, S_in = S_in_powered,
    age_distribution, M, N_1, N_2,
    population_size, degree_distribution,
    target_D = 3, target = target
  )

  susV <- vaccinated_sus_matrix(
    M, N_1, N_2,
    effV$effV_mat, effV$effV_age,
    x_state$S, x_state$xS, x_state$xI, x_state$xbar
  )

  time_varying_ngm(
    beta, p_cont, gamma,
    contact_matrix, age_distribution, age_distribution_sa,
    n_active, degree_distribution, degree_sizebiased, D_max,
    M, N_1, N_2,
    susV$S_vac, susV$xS_vac, susV$xI_vac, susV$xbar_vac
  )
}

#' Compute cumulative elasticity from original and vaccinated NGMs
#'
#' @param ngm matrix[2M,2M] - original NGM
#' @param ngm_vac matrix[2M,2M] - vaccinated NGM
#' @param M int - the number of age groups
#'
#' @return A numeric vector of cumulative elasticities by age group
#' @export
compute_cumulative_elasticity <- function(ngm, ngm_vac, M) {
  ela <- elas_func(
    matA = ngm,
    delta_matA = ngm - ngm_vac
  )

  cum_ela <- colSums(ela[, 1:M] + ela[, (M + 1):(2 * M)])

  return(cum_ela)
}

#' Compute cumulative sensitivity from original and vaccinated NGMs
#'
#' @param ngm matrix[2M,2M] - original NGM
#' @param ngm_vac matrix[2M,2M] - vaccinated NGM
#' @param M int - the number of age groups
#'
#' @return A numeric vector of cumulative elasticities by age group
#' @export
compute_cumulative_sensitivity <- function(ngm, ngm_vac, M) {
  sens <- sens_func(
    matA = ngm,
    delta_matA = ngm - ngm_vac
  )

  sens_ela <- colSums(sens[, 1:M] + sens[, (M + 1):(2 * M)])

  return(sens_ela)
}

#' Compute cumulative elasticity or sensitivity over time
#'
#' @param model_output Output from model_mpox_network with x_output = TRUE
#' @param max_t Integer. Maximum time point (default = 200)
#' @param method Character. Either "elasticity" or "sensitivity"
#' @param outcome Character. Either "infection" or "death"
#' @param gamma Numeric. Recovery rate
#' @param beta Numeric. Transmission rate
#' @param matM Numeric matrix. Mortality matrix
#'
#' @return A data frame with columns: time, age_group, and value
#' @export
compute_time_series_ela_or_sens <- function(model_output, max_t = 200, outcome = "infection", method = "elasticity",
                                            gamma = gamma, beta = beta, matM = matM) {
  pert_data <- lapply(1:max_t, function(t) {
    x_state <- extract_x_arrays_at_time(model_output, t, n_active, D_max)

    ngm <- time_varying_ngm(
      beta, p_cont, gamma,
      contact_matrix,
      age_distribution, age_distribution_sa,
      n_active, degree_distribution, degree_sizebiased, D_max,
      M, N_1, N_2,
      x_state$S, x_state$xS, x_state$xI, x_state$xbar
    )$global_NGM_large

    S_in <- compute_S_in(x_state$xbar, age_distribution_sa)
    S_in_powered <- compute_S_in_powered(S_in)

    ngm_vac <- compute_vaccinated_ngm(S_in_powered, x_state,
      gamma = gamma, beta = beta
    )$global_NGM_large

    # define projection matrix by outcome
    if (outcome == "infection") {
      matP <- ngm
      matP_vac <- ngm_vac
    } else if (outcome == "death") {
      matP <- matM %*% ngm
      matP_vac <- matM %*% ngm_vac
    } else {
      stop("Invalid outcome. Choose 'infection' or 'death'.")
    }

    # compute cumulative elasticity or sensitivity
    if (method == "elasticity") {
      value <- compute_cumulative_elasticity(matP, matP_vac, M)
    } else if (method == "sensitivity") {
      value <- compute_cumulative_sensitivity(matP, matP_vac, M)
    } else {
      stop("Invalid method. Choose 'elasticity' or 'sensitivity'.")
    }

    data.frame(
      time = t,
      age_group = 1:M,
      value = value
    )
  })

  df_pert_all <- do.call(rbind, pert_data)
  return(df_pert_all)
}


#' Format data frames for infection and mortality plots at t = 0 and t = time_obs
#'
#' @param ela_0 Vector of cumulative elasticities at t = 0
#' @param sens_0 Vector of sensitivities at t = 0
#' @param ela_0_matM Vector of elasticities with mortality-weighted matrix at t = 0
#' @param sens_0_matM Vector of sensitivities with mortality-weighted matrix at t = 0
#' @param df_inf_elas_overtime Data frame of infection elasticities over time
#' @param df_inf_sens_overtime Data frame of infection sensitivities over time
#' @param df_mort_elas_overtime Data frame of mortality elasticities over time
#' @param df_mort_sens_overtime Data frame of mortality sensitivities over time
#' @param time_obs Numeric. Time point at which to extract the observed values (default: 200)
#' @param M Integer. Number of age groups
#' @param age_labels Character vector. Age group labels for plotting
#'
#' @return A list of reshaped data frames suitable for ggplot bar charts
#' @export
format_barplot_data <- function(
    ela_0, sens_0,
    ela_0_matM, sens_0_matM,
    df_inf_elas_overtime, df_inf_sens_overtime,
    df_mort_elas_overtime, df_mort_sens_overtime,
    time_obs = 200,
    M,
    age_labels = c(
      "0 to 4", "5 to 9", "10 to 14", "15 to 19", "20 to 24",
      "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50+"
    )) {
  # Validate age labels length
  stopifnot(length(age_labels) == M)

  to_long <- function(vec, label, type_label, time_label) {
    data.frame(
      group = age_labels,
      type = type_label,
      value = vec,
      label = label,
      time = time_label
    )
  }

  # Reshape all relevant vectors/data frames into long format
  df_0_inf <- rbind(
    to_long(df_inf_sens_overtime |> dplyr::filter(time == time_obs) |> dplyr::pull(value), "Infection", "sensitivity", time_obs),
    to_long(df_inf_elas_overtime |> dplyr::filter(time == time_obs) |> dplyr::pull(value), "Infection", "elasticity", time_obs)
  )

  df_0_mort <- rbind(
    to_long(df_mort_sens_overtime |> dplyr::filter(time == time_obs) |> dplyr::pull(value), "Mortality", "sensitivity", time_obs),
    to_long(df_mort_elas_overtime |> dplyr::filter(time == time_obs) |> dplyr::pull(value), "Mortality", "elasticity", time_obs)
  )

  df_t0_inf <- rbind(
    to_long(ela_0, "Infection", "elasticity", 0),
    to_long(sens_0, "Infection", "sensitivity", 0)
  )

  df_t0_mort <- rbind(
    to_long(ela_0_matM, "Mortality", "elasticity", 0),
    to_long(sens_0_matM, "Mortality", "sensitivity", 0)
  )

  # Combine and reshape for plotting
  df_all <- dplyr::bind_rows(df_0_inf, df_0_mort, df_t0_inf, df_t0_mort) |>
    dplyr::mutate(value = ifelse(type == "elasticity", -value, value))

  return(df_all)
}
