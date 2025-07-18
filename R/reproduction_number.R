#' The basic reproduction number for the mpox network model
#'
#' @inheritParams model_mpox_network
#'
#' @return list with the following list elements
#' - reproduction_number - the basic reproduction number
#' - R_ss - the expected number of secondary cases from sexual to sexual transmission
#' - R_gg - the expected number of secondary cases from global to global transmission
#'
#' @export
#'
#' @examples reproduction_number(
#'   beta = 0.2, gamma = 1 / 5, p_cont = 0.05,
#'   degree_distribution = c(0.1, 0.3, 0.5, 0.05),
#'   age_distribution = c(0.1, 0.9),
#'   contact_matrix = matrix(1, nrow = 2, ncol = 2),
#'   M = 2, N_1 = 1, N_2 = 1
#' )
reproduction_number <- function(beta, gamma, p_cont, degree_distribution, age_distribution, contact_matrix, M, N_1, N_2) {
  D_max <- length(degree_distribution)
  degree_sizebiased <- 1:D_max * degree_distribution
  degree_sizebiased <- degree_sizebiased / sum(degree_sizebiased)

  age_distribution_sa <- age_distribution[N_1:N_2] / sum(age_distribution[N_1:N_2])

  sexact_age_dist <- rep(0, length = M)
  sexact_age_dist[N_1:N_2] <- age_distribution_sa

  k_ss <- beta / (beta + gamma) * sum(0:(D_max - 1) * degree_sizebiased)

  age_matrix <- matrix(age_distribution, nrow = M, ncol = M, byrow = FALSE)
  k_gg <- p_cont / gamma * age_matrix * contact_matrix

  k_gs <- rep(0, length = M)
  for (i in 1:M) {
    contact_with_i_from_sa <- sum(contact_matrix[i, N_1:N_2] * age_distribution_sa)
    k_gs[i] <- p_cont / gamma * age_distribution[i] * contact_with_i_from_sa
  }

  indicator_sa <- rep(0, length = M)
  indicator_sa[N_1:N_2] <- 1
  k_sg <- beta / (beta + gamma) * sum(1:D_max * degree_distribution) * indicator_sa ### check### #sexact_age_dist?


  k1 <- cbind(k_ss, t(k_sg))
  k2 <- cbind(k_gs, k_gg) 

  ngm_matrix <- rbind(k1, k2) 
  reproduction_number <- max(Re(eigen(ngm_matrix)$values))

  # dominant eigenvalue of k_gg matrix: global contacts - mass action model
  R_gg <- max(Re(eigen(k_gg)$values))

  return(list(
    reproduction_number = reproduction_number,
    R_ss = k_ss,
    R_gg = R_gg
  ))
}
