#' Large age-structured next generation matrix - refer to Supplement 1.1.1
#'
#' @inheritParams model_mpox_network
#'
#' @return list
#' @export
NGM_fn <- function(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2) {
  # index of sexual active age groups (0 inactive; 1 active)
  indicator_sa <- matrix(0, nrow = M, ncol = M)
  indicator_sa[, N_1:N_2] <- 1

  # distribution of sexually active age groups
  sexact_age_dist <- rep(0, length = M)
  sexact_age_dist[N_1:N_2] <- age_distribution[N_1:N_2] / sum(age_distribution[N_1:N_2])

  # age-weighted contact matrix
  age_matrix <- matrix(age_distribution, nrow = length(age_distribution), ncol = length(age_distribution), byrow = FALSE)
  sym_contact_home_matrix_ageW <- (age_matrix * contact_matrix)
  # c_ave <- eigen(sym_contact_home_matrix_ageW)$value[1]

  # degree distribution
  first_moment <- degree_distribution_moments(degree_distribution)$first_moment # para[4] #1.9763
  second_moment <- degree_distribution_moments(degree_distribution)$second_moment # para[5] #16.2215

  # next generation matrices (k_gg, k_gs, k_sg, k_ss)
  k_gg_large <- p_cont * (1 / gamma) * sym_contact_home_matrix_ageW

  k_gs_large <- k_gg_large * indicator_sa

  k_sg_large <- rbind(
    matrix(rep(0, M), nrow = N_1 - 1, ncol = M, byrow = TRUE),
    matrix((beta / (beta + gamma)) * first_moment * sexact_age_dist, nrow = N_2 - N_1 + 1, ncol = M, byrow = TRUE),
    matrix(rep(0, M), nrow = M - N_2, ncol = M, byrow = TRUE)
  ) |> t() # this transpose is needed to be in line with K_sg of Appendix

  k_ss_large <- rbind(
    matrix(rep(0, M), nrow = N_1 - 1, ncol = M, byrow = TRUE),
    matrix((beta / (beta + gamma)) * ((second_moment - first_moment) / first_moment) * sexact_age_dist, nrow = N_2 - N_1 + 1, ncol = M, byrow = TRUE),
    matrix(rep(0, M), nrow = M - N_2, ncol = M, byrow = TRUE)
  ) |> t() # this transpose is needed to be in line with K_ss of Appendix

  # merge all matrices into one NGM
  global_NGM_large_toprow <- cbind(k_ss_large, k_sg_large) # eq(3)

  global_NGM_large_secondrow <- cbind(k_gs_large, k_gg_large)

  global_NGM_large <- rbind(global_NGM_large_toprow, global_NGM_large_secondrow)

  t_global_NGM_large <- t(global_NGM_large)

  # Extract the top right eigenvector
  topr_eigenvec_NGM <- abs(Re(eigen((t_global_NGM_large))$vectors[, 1])) # the top right eigenvector of the global NGM

  # Rearrange this vector as the age-distribution of incidence of infection
  topr_eigenvec_age <- (topr_eigenvec_NGM[1:M] + topr_eigenvec_NGM[(M + 1):(2 * M)]) / sum(topr_eigenvec_NGM[1:M] + topr_eigenvec_NGM[(M + 1):(2 * M)]) # incidence of infection by age, normalized
  asym_incdist_age_NGM <- (topr_eigenvec_age * age_distribution) / sum(topr_eigenvec_age * age_distribution)

  reproduction_number <- max(Re(eigen(t_global_NGM_large)$values))
  R_gg <- max(Re(eigen(k_gg_large)$values))
  R_ss <- max(Re(eigen(k_ss_large)$values))
  R_sg <- max(Re(eigen(k_sg_large)$values))
  R_gs <- max(Re(eigen(k_gs_large)$values))


  # store all results in a list
  results <- list(
    reproduction_number = reproduction_number,
    R_gg = R_gg,
    R_ss = R_ss,
    R_sg = R_sg,
    R_gs = R_gs,
    asym_incdist_age_NGM = asym_incdist_age_NGM,
    k_gg_large = k_gg_large,
    k_gs_large = k_gs_large,
    k_sg_large = k_sg_large,
    k_ss_large = k_ss_large,
    t_global_NGM_large = t_global_NGM_large,
    first_moment = first_moment,
    second_moment = second_moment,
    weight = ((second_moment - first_moment) / first_moment)
  )

  return(results)
}
