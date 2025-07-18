#' Grouping of the oldest age groups into one age group in raw contact matrices for the DRC
#'
#' @param raw_contact_matrix matrix - raw contact matrix
#' @param raw_age_dist vector - age distribution corresponding to the age groups in the raw contact matrix
#' @param M int - total number of desired age groups
#' @param N_1 int - first sexually active age group
#' @param N_2 int - last sexually active age group, 1 <= N_1 <= N_2 <= M
#'
#' @return matrix
#' @export
prepare_contact_matrix <- function(raw_contact_matrix, raw_age_dist, N_1, N_2, M) {
  contact_matrix <- (raw_contact_matrix$COG + t(raw_contact_matrix$COG)) / 2 # symmetric contact matrix

  # age distribution from the DRC pop data, adjusted to age groups
  age_distribution <- c(raw_age_dist[1:N_2], 1 - sum(raw_age_dist[1:N_2]))

  # group contact rates for age groups M+
  m <- rep(0, length = M)
  n_age_groups <- length(raw_age_dist) # total number of age groups in raw data
  # total contact rate between individuals in group i and individuals in group j, j = M, ..., n_age_group
  for (i in 1:M) {
    for (j in M:n_age_groups) {
      m[i] <- m[i] + contact_matrix[i, j] * raw_age_dist[i] * raw_age_dist[j]
    }
  }
  # total age probability of groups M:n_age_groups
  age_probability_grouped <- sum(raw_age_dist[M:n_age_groups])
  # symmetric contact matrix with rates c[i,M], c[M,i] representing the grouped contact rates
  contact_matrix <- contact_matrix[1:M, 1:M]
  for (i in 1:M) {
    contact_matrix[i, M] <- m[i] / (raw_age_dist[i] * age_probability_grouped)
    contact_matrix[M, i] <- contact_matrix[i, M]
  }
  # check that contact matrix is symmetric
  stopifnot("contact matrix is not symmetric" = all(contact_matrix == t(contact_matrix)))
  return(contact_matrix)
}
