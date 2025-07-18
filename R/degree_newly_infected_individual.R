#' Degree distribution of a newly infected individual as a function of time
#'
#' @inheritParams time_varying_ngm
#'
#' @return matrix
#' @export
degree_newly_infected_individual <- function(n_active, D_max, degree_distribution, degree_sizebiased, age_distribution_sa, xI, xbar) {
  qs <- matrix(nrow = n_active, ncol = D_max)
  for (i in 1:n_active) {
    for (n in 1:D_max) {
      qs[i, n] <- n * degree_distribution[n] * sum(age_distribution_sa * xI[i, , n]) * sum(age_distribution_sa * xbar[i, , n])^(n - 1)
    }
  }

  # # normalize degree distribution
  for (i in 1:n_active) {
    if (sum(qs[i, ]) > 0) {
      qs[i, ] <- qs[i, ] / sum(qs[i, ])
    } else {
      qs[i, ] <- degree_sizebiased
    }
  }
  return(qs)
}
