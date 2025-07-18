#' Effective vaccination coverage for age group i with degree n
#'
#' @param deltaU double - unit change in the number of vaccines
#' @param VE double - vaccine efficacy
#' @param S_in matrix[n_active, D_max] - proportion susceptible in group (i,n) (the group in age i with degree n)
#' @param age_distribution vector[M] - population age distribution
#' @param M int - total number of age groups
#' @param N_1 int - first sexually active age group
#' @param N_2 int - last sexually active age group, 1 <= N_1 <= N_2 <= M
#' @param population_size double - population size
#' @param degree_distribution vector - defective distribution of degrees 1,..., D_max, excluding degree_zero
#' @param target_D double - integer, minimum degree class to target with vaccination (used in "age-degree" targeting)
#' @param target character - targeting strategy; one of "random", "age", "age-degree"
#'
#' @return list with two elements
#' - effV_mat = matrix[n_active, D_max] - effective vaccine coverage for group (i,n)
#' - effV_age = vector[M] - effective vaccine coverage for age group i (all age groups 1,..,M)
#' @export
eff_vaccov_matrix <- function(deltaU = 1,
                              VE = 0.8,
                              S_in,
                              age_distribution,
                              M,
                              N_1,
                              N_2,
                              population_size,
                              degree_distribution,
                              target_D = NULL,
                              target = c("random", "age", "age-degree")) {
  # sexually active age groups and max degree
  n_active <- N_2 - N_1 + 1
  D_max <- length(S_in[1, ])

  # check if dimension is correct
  stopifnot("Dimension of S_in does not match the number of sexually active age groups" = n_active == length(S_in[, 1]))

  # effective vaccination coverage for group (i,n), i=N_1,N_1+1,...,N_2
  effV_mat <- matrix(0, nrow = n_active, ncol = D_max)

  # effective vaccination coverage for age group i, i=1,2,...,M
  effV_age <- rep(0, M)

  if (target == "random") {
    # sexual transmission
    for (i in 1:n_active) {
      for (n in 1:D_max) {
        effV_mat[i, n] <- VE * deltaU / (S_in[i, n] * population_size)
      }
    }
    # global transmission
    effV_age <- rep(VE * deltaU / population_size, M)
  } else if (target == "age") {
    # sexual transmission)
    for (i in 1:n_active) {
      for (n in 1:D_max) {
        effV_mat[i, n] <- VE * deltaU / (S_in[i, n] * age_distribution[(N_1 - 1 + i)] * population_size)
      }
    }
    # global transmission
    effV_age <- VE * deltaU / (age_distribution * population_size)
  } else if (target == "age-degree") {
    # sexual transmission
    p_T <- sum(degree_distribution[target_D:D_max])
    # for non-target degrees
    for (i in 1:n_active) {
      for (n in 1:(target_D - 1)) {
        effV_mat[i, n] <- 0
      }
    }
    # for non-target degrees
    for (i in 1:n_active) {
      for (n in target_D:D_max) {
        effV_mat[i, n] <- VE * deltaU / (p_T * S_in[i, n] * age_distribution[(N_1 - 1 + i)] * population_size)
      }
    }
    # global transmission
    # effective vaccination coverage for sexually-active group i
    effV_val_i <- c() # to store the result
    for (i in 1:n_active) {
      effV_val_i[i] <- sum(effV_mat[i, ] * S_in[i, ] * degree_distribution) / sum(S_in[i, ] * degree_distribution)
    }
    # for all age groups
    effV_age <- c(
      rep(0, (N_1 - 1)), # sexually inactive groups do not receive vaccines under this strategy
      effV_val_i,
      rep(0, (M - N_2)) # sexually inactive groups
    )
  } else {
    stop("target is not correctly specified")
  }

  # return results
  return(list(
    effV_mat = effV_mat,
    effV_age = effV_age
  ))
}
