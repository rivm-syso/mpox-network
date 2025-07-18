#' Mpox network model simulation of epidemic dynamics over time
#'
#' @param M int - total number of age groups
#' @param N_1 int - first sexually active age group
#' @param N_2 int - last sexually active age group, 1 <= N_1 <= N_2 <= M
#' @param contact_matrix matrix[M, M] - global contact rates, contact_matrix[i, j] = contact_matrix[j, i]
#' per-pair global contact rate between age groups i and j
#' @param age_distribution vector[M] - population age distribution
#' @param degree_zero double >=0 - fraction of sexually active population with degree zero
#' @param degree_distribution vector - defective distribution of degrees 1,..., D_max, excluding degree_zero
#' @param beta double >=0 - transmission rate sexual contact per day
#' @param p_cont double between 0 and 1 - transmission probability global contact
#' @param sigma double >= 0 - 1/sigma mean latent period in days
#' @param gamma double >= 0 - recovery rate per day
#' @param seed_infected double between 0 and 1 - initial infected
#' @param time_end int - end time of simulation in days
#' @param initial vector - user provided initial value condition,
#' default is NULL corresponding to uniform distribution of infection over the age groups
#' @param time_increment double - increment time steps in simulation
#'
#' @return tibble
#' @export
#'
#' @examples data <- model_mpox_network(
#'   M = 4, N_1 = 2, N_2 = 3,
#'   contact_matrix = matrix(c(c(1, 0, 0, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 0, 0, 1)), nrow = 4),
#'   age_distribution = c(1 / 8, 3 / 8, 2 / 8, 2 / 8),
#'   degree_zero = 0,
#'   degree_distribution = c(0, 0, 0.3, 0.7),
#'   beta = 0.1, p_cont = 1, sigma = 1 / 5, gamma = 1 / 5,
#'   seed_infected = 1E-3,
#'   time_end = 600,
#'   time_increment = 1
#' )
model_mpox_network <- function(M,
                               N_1,
                               N_2,
                               contact_matrix,
                               age_distribution,
                               degree_zero,
                               degree_distribution,
                               beta,
                               p_cont,
                               sigma,
                               gamma,
                               seed_infected,
                               time_end,
                               time_increment,
                               initial = NULL,
                               x_output = FALSE) {
  #### perform some checks on input
  # check that global transmission probability is between 0 and 1
  stopifnot("Error: global transmission probability should be smaller or equal to 1" = p_cont <= 1)
  # check that degree distribution sums to one
  stopifnot("Error: degree distribution does not sum to one" = degree_zero + sum(degree_distribution) == 1)
  # check that at least one degree > 0
  stopifnot("Error: there should be at least one degree > 0" = sum(degree_distribution) > 0)
  # check that dimensions are correct
  stopifnot("Error: contact matrix dimensions do not match number of age groups" = (nrow(contact_matrix) == M & ncol(contact_matrix) == M))
  stopifnot("Error: age distribution does not match number of age groups" = length(age_distribution) == M)
  # check that age distribution sums to one
  stopifnot("Error: age distribution does not sum to one" = sum(age_distribution) == 1)
  # # age distribution can't have zeroes
  # stopifnot("Error: age distribution cannot have zeroes" = all(age_distribution > 0))
  # check that degree distribution has degree at least degree 2 for giant component
  stopifnot("Error: there should be degree at of at least two" = length(degree_distribution) >= 2)
  # check that the max degree has probability larger than zero
  stopifnot("Error: the maximum degree should have probability larger than zero" = degree_distribution[length(degree_distribution)] > 0)
  # 1 <= N_1 <= N_2 <= M
  stopifnot("Error: check definition of sexually active age groups N_1 and N_2" = (1 <= N_1 & N_1 <= N_2 & N_2 <= M))
  # contact matrix needs to be >0
  stopifnot("Error: all global contacts equal to zero is not allowed" = any(contact_matrix != 0))

  # sexually active age groups
  n_active <- N_2 - N_1 + 1
  # age distribution in sexually active ages
  # only take into account sexually active ages and index from first to last sexually active age group
  age_distribution_sa <- age_distribution[N_1:N_2] / sum(age_distribution[N_1:N_2])
  # check normalization
  stopifnot(
    "Error: normalization of sexually active age distribution is not one" =
      abs(sum(age_distribution_sa) - 1) < 1E-6
  )

  # size biased degree distribution
  D_max <- length(degree_distribution)
  degree_sizebiased <- degree_distribution
  for (k in 1:D_max) {
    degree_sizebiased[k] <- k * degree_distribution[k]
  }
  degree_sizebiased <- degree_sizebiased / sum(degree_sizebiased)

  # named list of all model parameters
  params <- list(
    M = M,
    N_1 = N_1,
    N_2 = N_2,
    n_active = n_active,
    contact_matrix = contact_matrix,
    age_distribution = age_distribution,
    age_distribution_sa = age_distribution_sa,
    degree_zero = degree_zero,
    degree_distribution = degree_distribution,
    degree_sizebiased = degree_sizebiased,
    D_max = D_max,
    beta = beta,
    p_cont = p_cont,
    sigma = sigma,
    gamma = gamma
  )

  if (is.null(initial)) {
    # setup initial values
    # far past conditions: everyone is susceptible
    S_initial <- rep(1 - seed_infected, length = M)
    E_initial <- rep(seed_infected, length = M)
    I_initial <- rep(0, length = M)
    x0_initial <- rep(0, length = n_active)
    xbar_initial <- array(0, dim = c(n_active, n_active, D_max))
    xE_initial <- array(0, dim = c(n_active, n_active, D_max))
    xI_initial <- array(0, dim = c(n_active, n_active, D_max))
    xS_initial <- array(0, dim = c(n_active, n_active, D_max, D_max))

    # perturbation of the disease free steady state:
    for (i in 1:n_active) {
      if (degree_zero > 0) {
        x0_initial[i] <- 1 - seed_infected
      }
    }

    for (i in 1:n_active) {
      for (n in 1:D_max) {
        if (degree_distribution[n] > 0) {
          for (j in 1:n_active) {
            xE_initial[i, j, n] <- seed_infected
            xbar_initial[i, j, n] <- 1 - seed_infected
            for (k in 1:D_max) {
              if (degree_distribution[k] > 0) {
                xS_initial[i, j, n, k] <- degree_sizebiased[k] * (1 - 2 * seed_infected)
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
  }

  # simulate outbreak
  model_output <- deSolve::ode(
    y = initial_values,
    times = seq(0, time_end, by = time_increment),
    func = ode_model_mpox_network,
    parms = params,
    method = "ode23" # "ode23" # "lsoda" #"ode23" # rk, rk4, ode45
  )
  # convert class of output and type of
  model_output <- model_output |>
    as.data.frame() |>
    # change deSolve classes to numeric
    type.convert(as.is = TRUE)
  # rename recovered class and sexual age groups from deSolve output
  model_output <- model_output |>
    dplyr::rename_with(
      .cols = dplyr::starts_with("Recovered"),
      .fn = ~ gsub(".Susceptible", "", .x)
    ) |>
    dplyr::rename_with(
      .cols = dplyr::starts_with("foi_sexual"),
      .fn = ~ gsub(".Susceptible", "", .x)
    )

  if (x_output == FALSE) {
    # remove columns that are not used in subsequent analysis
    model_output <- model_output |>
      select(-dplyr::starts_with("x"))
  }

  return(model_output)
}


#' ODE system for the mpox network model
#'
#' @param times vector of times
#' @param y vector of variables
#' @param parms list of parameters
#'
#' @return list of solution to the ODE system over time
#' @export
ode_model_mpox_network <- function(times, y, parms) {
  with(parms, {
    # get state variables
    S <- y[grepl(x = names(y), pattern = "Susceptible")]
    E <- y[grepl(x = names(y), pattern = "Exposed")]
    I <- y[grepl(x = names(y), pattern = "Infectious")]
    x0 <- y[grepl(x = names(y), pattern = "xsingle")]
    xbar <- y[grepl(x = names(y), pattern = "xbar")]
    xS <- y[grepl(x = names(y), pattern = "xS")]
    xE <- y[grepl(x = names(y), pattern = "xE")]
    xI <- y[grepl(x = names(y), pattern = "xI")]
    # convert back to matrix and 3/4-dim array for the relevant x-variables to
    # access variables through multiple indices for age group and degree
    xbar <- array(xbar, dim = c(n_active, n_active, D_max))
    xS <- array(xS, dim = c(n_active, n_active, D_max, D_max))
    xE <- array(xE, dim = c(n_active, n_active, D_max))
    xI <- array(xI, dim = c(n_active, n_active, D_max))

    # sexually active force of infection
    lambda_s_helper <- matrix(0, nrow = n_active, ncol = D_max)
    for (i in 1:n_active) {
      for (n in 1:D_max) {
        lambda_s_helper[i, n] <- degree_distribution[n] * beta * n * sum(age_distribution_sa * xbar[i, , n])^(n - 1) * sum(age_distribution_sa * xI[i, , n])
      }
    }
    lambda_s <- rep(0, length = n_active)
    for (i in 1:n_active) {
      lambda_s[i] <- sum(lambda_s_helper[i, ])
    }

    # sexually active mean field at distance one
    Lambda_s <- matrix(0, nrow = n_active, ncol = D_max)
    for (j in 1:n_active) {
      for (m in 2:D_max) {
        if (degree_distribution[m] > 0) {
          Lambda_s[j, m] <- beta * (m - 1) * sum(age_distribution_sa * xI[j, , m]) / sum(age_distribution_sa * xbar[j, , m])
        }
      }
    }
    # global force of infection
    lambda_g <- p_cont * age_distribution * contact_matrix %*% I

    # global mean field at distance one
    Lambda_g <- matrix(0, nrow = n_active, ncol = D_max)
    for (j in 1:n_active) {
      for (m in 1:D_max) {
        if (degree_distribution[m] > 0) {
          Lambda_g[j, m] <- sum(age_distribution_sa * xbar[j, , m])^(m - 1) * lambda_g[j + (N_1 - 1)]
        }
      }
    }

    ############################# ODEs ########################################
    # initialize empty arrays with the correct dimensions
    dxS <- array(dim = c(n_active, n_active, D_max, D_max))
    dxE <- array(dim = c(n_active, n_active, D_max))
    dxI <- array(dim = c(n_active, n_active, D_max))
    dxbar <- array(dim = c(n_active, n_active, D_max))
    dx0 <- vector(length = n_active)
    dS <- vector(length = M)
    dE <- vector(length = M)
    dI <- vector(length = M)

    # Note indexing is either from 1:n_active or from 1:M so in some cases the
    # index needs to be moved by N_1-1
    # use i,j for indexing age groups and n, k for indexing degrees
    ########### susceptible binding site dynamics
    for (i in 1:n_active) {
      for (j in 1:n_active) {
        for (n in 1:D_max) {
          if (degree_distribution[n] > 0) {
            dxbar[i, j, n] <- -beta * xI[i, j, n] - 1 / n * xbar[i, j, n] * lambda_g[i + (N_1 - 1)]
            dxE[i, j, n] <- sum(Lambda_s[j, ] * xS[i, j, n, ]) + sum(Lambda_g[j, ]) - sigma * xE[i, j, n] - 1 / n * xE[i, j, n] * lambda_g[i + (N_1 - 1)]
            dxI[i, j, n] <- sigma * xE[i, j, n] - (beta + gamma) * xI[i, j, n] - 1 / n * xI[i, j, n] * lambda_g[i + (N_1 - 1)]
          } else {
            dxbar[i, j, n] <- 0
            dxE[i, j, n] <- 0
            dxI[i, j, n] <- 0
          }
        }
      }
    }

    for (i in 1:n_active) {
      for (j in 1:n_active) {
        for (n in 1:D_max) {
          for (m in 1:D_max) {
            if (degree_distribution[n] > 0 & degree_distribution[m] > 0) {
              dxS[i, j, n, m] <- -Lambda_s[j, m] * xS[i, j, n, m] - 1 / n * xS[i, j, n, m] * lambda_g[i + (N_1 - 1)] - Lambda_g[j, m]
            } else {
              dxS[i, j, n, m] <- 0
            }
          }
        }
      }
    }

    for (i in 1:n_active) {
      if (degree_zero > 0) {
        dx0[i] <- -x0[i] * lambda_g[i + (N_1 - 1)]
      } else {
        dx0[i] <- 0
      }
    }

    ############ individual dynamics
    ## sexually active age groups
    # note lambda_s is indexed from 1:n_active so need to change the index by N_1-1
    for (i in N_1:N_2) {
      dS[i] <- -lambda_s[i - (N_1 - 1)] - lambda_g[i] * S[i]
      dE[i] <- lambda_s[i - (N_1 - 1)] + lambda_g[i] * S[i] - sigma * E[i]
      dI[i] <- sigma * E[i] - gamma * I[i]
    }
    ## check if there are sexually inactive ages and if so calculate epidemic dynamics
    ## sexually inactive: young 1:(N_1-1) and old (N_2+1):M
    # first sexually active age group is older than first age group
    if (N_1 > 1) {
      for (i in 1:(N_1 - 1)) {
        dS[i] <- -lambda_g[i] * S[i]
        dE[i] <- lambda_g[i] * S[i] - sigma * E[i]
        dI[i] <- sigma * E[i] - gamma * I[i]
      }
    }
    # last sexually active age group is younger than last age group
    if (N_2 < M) {
      for (i in (N_2 + 1):M) {
        dS[i] <- -lambda_g[i] * S[i]
        dE[i] <- lambda_g[i] * S[i] - sigma * E[i]
        dI[i] <- sigma * E[i] - gamma * I[i]
      }
    }

    ############################################################################
    # additional non-state variables to track
    # fraction of recovered individuals
    R <- 1 - (S + E + I)

    # incidence of global transmission: index 1 = age group 1
    incidence_global <- S * lambda_g
    # incidence of sexual transmission: note index 1 = age group N_1 etc
    incidence_sexual <- lambda_s

    # force of infection
    foi_global <- incidence_global / S
    foi_sexual <- incidence_sexual / S[N_1:N_2]

    # probability that the degree of a newly infected sexually individual is n
    # first term: infected through sexual transmission route
    # second term: infected through global transmission route
    # summed over all sexually active age groups
    degree_newly_infected_sa <- vector(length = D_max)
    for (n in 1:D_max) {
      degree_newly_infected_sa[n] <- sum(age_distribution_sa * lambda_s_helper[, n]) +
        sum(age_distribution_sa * lambda_g[N_1:N_2]) * degree_distribution[n]
    }
    # normalize into probability distribution
    degree_newly_infected_sa <- degree_newly_infected_sa / sum(degree_newly_infected_sa)
    # calculate mean
    mean_degree_newly_infected <- sum(1:D_max * degree_newly_infected_sa)

    qs <- degree_newly_infected_individual(n_active, D_max, degree_distribution, degree_sizebiased, age_distribution_sa, xI, xbar)
    PS <- probability_susceptible_partner(n_active, D_max, age_distribution_sa, degree_sizebiased, xbar, xS)

    # dominant eigenvalue of the time-varying NGM
    R_t <- time_varying_ngm(
      beta, p_cont, gamma,
      contact_matrix,
      age_distribution, age_distribution_sa, n_active,
      degree_distribution, degree_sizebiased, D_max,
      M, N_1, N_2,
      S, xS, xI, xbar
    )$R_all_t
    # output
    return(list(
      c(dS, dE, dI, dx0, dxbar, dxS, dxE, dxI),
      c(
        Recovered = R,
        incidence_global = incidence_global,
        incidence_sexual = incidence_sexual,
        foi_global = foi_global,
        foi_sexual = foi_sexual,
        mean_degree_newly_infected = mean_degree_newly_infected,
        R_t = R_t,
        qs = qs,
        PS = PS
      )
    ))
  })
}
