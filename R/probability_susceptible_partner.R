#' The probability that a partner of a newly infected individual is susceptible
#'
#' @description refer to manuscript and supplement for explanation and derivation of this quantity
#'
#' @inheritParams time_varying_ngm
#'
#' @return matrix of dimension `n_active`, `D_max` with probabilities
#' @export
probability_susceptible_partner <- function(n_active, D_max, age_distribution_sa, degree_sizebiased, xbar, xS) {
  # helper matrix
  susceptible_rest <- matrix(nrow = n_active, ncol = D_max)
  for (j in 1:n_active) {
    for (k in 1:D_max) {
      susceptible_rest[j, k] <- sum(age_distribution_sa * xbar[j, , k])^(k - 1)
    }
  }

  helper_S <- array(dim = c(n_active, n_active, D_max))
  for (i in 1:n_active) {
    for (j in 1:n_active) {
      for (n in 1:D_max) {
        helper_S[i, j, n] <- sum(xS[i, j, n, ])
      }
    }
  }
  prob_susceptible_partner <- matrix(nrow = n_active, ncol = n_active)
  for (i in 1:n_active) {
    for (j in 1:n_active) {
      prob_susceptible_partner[i, j] <- sum(degree_sizebiased * helper_S[j, i, ] * susceptible_rest[j, ])
    }
  }
  return(prob_susceptible_partner)
}
