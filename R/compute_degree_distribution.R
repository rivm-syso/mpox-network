#' Degree distribution excluding zero-degree individuals
#'
#' Computes the normalized degree distribution for individuals with one or more sexual partners,
#' based on a power-law (Pareto) form, with a proportion N_0 of individuals having zero partners.
#'
#' @param k_max Integer. Maximum number of sexual partners (degree).
#' @param N_0 Double in [0, 1). Proportion of individuals with zero partners.
#' @param alpha Double > 0. Pareto exponent for power-law degree distribution.
#'
#' @return A numeric vector of length `k_max` with probabilities for degrees 1 to `k_max`.
#' @export
#'
#' @examples
#' compute_degree_distribution(k_max = 25, N_0 = 0.2, alpha = 3.0)
compute_degree_distribution <- function(k_max, N_0, alpha) {
  # Zero-degree probability
  p_0 <- N_0

  # Power-law distribution for k >= 1
  p_k <- (1:k_max)^(-alpha)
  p_k <- p_k / sum(p_k) # Normalize

  # Scale non-zero part by (1 - N_0)
  scaled_p_k <- (1 - p_0) * p_k

  return(scaled_p_k) # degrees 1 to k_max
}
