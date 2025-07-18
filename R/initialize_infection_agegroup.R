#' Initialize infection in specified age group
#' 
#' @inheritParams simulate_model
#' @inheritParams time_varying_ngm
#' @param degree_zero double - fraction of individuals with degree zero
#' @param initial_age_group int - index of age group where infection is seeded
#'
#' @returns named vector
#' @export
initialize_infection_agegroup <- function(M, n_active, D_max, seed_infected, degree_zero, initial_age_group, degree_distribution, degree_sizebiased) {
  # setup initial values
  # far past conditions: everyone is susceptible
  S_initial <- rep(1, length = M)
  E_initial <- rep(0, length = M)
  I_initial <- rep(0, length = M)
  x0_initial <- rep(1, length = n_active)
  xbar_initial <- array(1, dim = c(n_active, n_active, D_max))
  xE_initial <- array(0, dim = c(n_active, n_active, D_max))
  xI_initial <- array(0, dim = c(n_active, n_active, D_max))
  xS_initial <- array(1, dim = c(n_active, n_active, D_max, D_max))

  # initial seed infected
  S_initial[initial_age_group] <- 1 - seed_infected
  E_initial[initial_age_group] <- seed_infected

  # perturbation of the disease free steady state:
  if (degree_zero > 0) {
    x0_initial[initial_age_group] <- 1 - seed_infected
  }

  for (i in 1:n_active) {
    for (n in 1:D_max) {
      if (degree_distribution[n] > 0) {
        xE_initial[i, initial_age_group, n] <- seed_infected
        xbar_initial[i, initial_age_group, n] <- 1 - seed_infected
        for (k in 1:D_max) {
          if (degree_distribution[k] > 0) {
            if (i != initial_age_group) {
              xS_initial[i, initial_age_group, n, k] <- degree_sizebiased[k] * (1 - seed_infected)
            } else {
              xS_initial[i, initial_age_group, n, k] <- degree_sizebiased[k] * (1 - 2 * seed_infected)
            }
            for (j in 1:n_active) {
              if (j != initial_age_group) {
                xS_initial[i, j, n, k] <- degree_sizebiased[k]
              }
            }
          }
        }
      }
    }
  }
  # collect initial values into one vector
  initial_values <- c(
    Susceptible = S_initial,
    Exposed = E_initial,
    Infectious = I_initial,
    xsingle = x0_initial,
    xbar = xbar_initial,
    xS = xS_initial,
    xE = xE_initial,
    xI = xI_initial
  )
  return(initial_values)
}
