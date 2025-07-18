#' Function to compute time varying NGMs
#'
#' @description refer to manuscript and supplement for details on derivation and equations
#'
#'
#' @inheritParams model_mpox_network
#' @param age_distribution_sa vector - age distribution of sexually active age groups, calculated from `age_distribution`
#' @param degree_sizebiased vector - size biased degree distribution, calculated from `degree_distribution`
#' @param S vector of length `M` of susceptible fractions by age group
#' @param xE array of dimension `n_active`, `n_active`, `D_max` of susceptible binding sites with exposed partner
#' @param xI array of dimension `n_active`, `n_active`, `D_max` of susceptible binding sites with infectious partner
#' @param xbar array of dimension `n_active`, `n_active`, `D_max` of susceptible binding sites
#' @param n_active int - number of sexually active age groups
#' @param D_max int - maximum degree in the population
#'
#' @return list with the following elements
#'  - `R_all_t` - time-varying reproduction number
#'  - `R_gg_t` - time-varying secondary cases from global to global transmission
#'  - `R_sg_t` - time-varying secondary cases from global to sexual transmission
#'  - `R_gs_t` - time-varying secondary cases from sexual to global transmission
#'  - `R_ss_t` - time-varying secondary cases from sexual to sexual transmission
#'  - `k_gg_large` - submatrix of the time-varying NGM consisting of secondary cases from global to global transmission
#'  - `k_sg_large` - submatrix of the time-varying NGM consisting of secondary cases from global to sexual transmission
#'  - `k_gs_large` - submatrix of the time-varying NGM consisting of secondary cases from sexual to global transmission
#'  - `k_ss_large` - submatrix of the time-varying NGM consisting of secondary cases from sexual to sexual transmission
#'  - `global_NGM_large` - time-varying NGM combining the four submatrices
#' @export
time_varying_ngm <- function(
    beta, p_cont, gamma,
    contact_matrix,
    age_distribution, age_distribution_sa, n_active,
    degree_distribution, degree_sizebiased, D_max,
    M, N_1, N_2,
    S, xS, xI, xbar) {
  # age-weighted contact matrix
  age_matrix <- matrix(age_distribution, nrow = M, ncol = M, byrow = FALSE)
  age_weighted_contact_matrix <- age_matrix * contact_matrix
  # index of sexual active age groups (0 inactive; 1 active)
  indicator_sa <- matrix(0, nrow = M, ncol = M)
  indicator_sa[, N_1:N_2] <- 1

  #### the probability that a partner of age j of a newly infected individual of age i is susceptible
  prob_susceptible_partner <- probability_susceptible_partner(n_active, D_max, age_distribution_sa, degree_sizebiased, xbar, xS)
  #### the probability that a newly infected individual of age i infected through sexual contact has degree n
  qs <- degree_newly_infected_individual(n_active, D_max, degree_distribution, degree_sizebiased, age_distribution_sa, xI, xbar)

  # rows are zero for non-sexually active age groups
  k_ss <- matrix(0, nrow = n_active, ncol = M)
  for (i in N_1:N_2) {
    for (j in 1:n_active) {
      k_ss[j, i] <- beta / (beta + gamma) * age_distribution_sa[j] * prob_susceptible_partner[i - (N_1 - 1), j] * sum(0:(D_max - 1) * qs[i - (N_1 - 1), ])
    }
  }
  # rows are zero for non-sexually active age groups
  k_sg <- matrix(0, nrow = n_active, ncol = M)
  for (i in N_1:N_2) {
    for (j in 1:n_active) {
      k_sg[j, i] <- beta / (beta + gamma) * age_distribution_sa[j] * prob_susceptible_partner[i - (N_1 - 1), j] * sum(1:D_max * degree_distribution)
    }
  }

  # NGM submatrices
  k_gg_large <- p_cont / gamma * age_weighted_contact_matrix %*% diag(S)
  k_gs_large <- k_gg_large * indicator_sa %*% diag(S)
  k_ss_large <- rbind(
    matrix(0, nrow = N_1 - 1, ncol = M),
    k_ss,
    matrix(0, nrow = M - N_2, ncol = M)
  )
  k_sg_large <- rbind(
    matrix(0, nrow = N_1 - 1, ncol = M),
    k_sg,
    matrix(0, nrow = M - N_2, ncol = M)
  )

  # merge all matrices into one NGM
  global_NGM_large_toprow <- cbind(k_ss_large, k_sg_large) # eq(3)

  global_NGM_large_secondrow <- cbind(k_gs_large, k_gg_large)

  global_NGM_large <- rbind(
    global_NGM_large_toprow,
    global_NGM_large_secondrow
  )

  # calculate the time varying reproduction number as the dominant eigenvalue of the NGM
  R_all_t <- max(Re(eigen(global_NGM_large)$values))

  R_gg_t <- max(Re(eigen(k_gg_large)$values))
  R_ss_t <- max(Re(eigen(k_ss_large)$values))
  R_sg_t <- max(Re(eigen(k_sg_large)$values))
  R_gs_t <- max(Re(eigen(k_gs_large)$values))

  # store all results in a list
  results <- list(
    R_all_t = R_all_t,
    R_gg_t = R_gg_t,
    R_ss_t = R_ss_t,
    R_sg_t = R_sg_t,
    R_gs_t = R_gs_t,
    k_gg_large = k_gg_large,
    k_gs_large = k_gs_large,
    k_sg_large = k_sg_large,
    k_ss_large = k_ss_large,
    global_NGM_large = global_NGM_large
  )

  return(results)
}
