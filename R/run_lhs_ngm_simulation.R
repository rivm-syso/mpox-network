#' Run LHS sampling and compute NGM outcomes
#'
#' @param n_samples Integer. Number of samples to simulate.
#' @param scaled_samples A numeric matrix of LHS samples (scaled to param ranges).
#' @param c_ave Double. Average household contact rate.
#' @param age_distribution Numeric vector. Population age distribution.
#' @param contact_matrix Numeric matrix. Age-specific contact matrix.
#' @param M Integer. Number of age groups.
#' @param N_1 Integer. Index of first sexually active age group.
#' @param N_2 Integer. Index of last sexually active age group.
#'
#' @return A data frame of simulation outcomes by sample and age group.
#' @export
run_lhs_ngm_simulation <- function(n_samples, scaled_samples, c_ave, age_distribution, contact_matrix, M, N_1, N_2) {
  compute_p_cont <- function(SAR_h, gamma, c_ave) {
    SAR_h / (1 - SAR_h) * (gamma / c_ave)
  }

  compute_beta <- function(SAR_s, gamma) {
    SAR_s * gamma / (1 - SAR_s)
  }

  results <- data.frame()

  for (i in 1:n_samples) {
    N_0_sample <- scaled_samples[i, 1]
    alpha_sample <- scaled_samples[i, 2]
    k_max_sample <- round(scaled_samples[i, 3])
    gamma_sample <- scaled_samples[i, 4]
    SAR_h_sample <- scaled_samples[i, 5]
    SAR_s_sample <- scaled_samples[i, 6]

    p_cont_sample <- compute_p_cont(SAR_h_sample, gamma_sample, c_ave)
    beta_sample <- compute_beta(SAR_s_sample, gamma_sample)

    degree_distribution_sample <- compute_degree_distribution(k_max_sample, N_0_sample, alpha_sample)

    ngm_result <- tryCatch(
      {
        NGM_fn(
          beta_sample, gamma_sample, p_cont_sample, degree_distribution_sample,
          age_distribution, contact_matrix, M, N_1, N_2
        )
      },
      error = function(e) NULL
    )

    if (!is.null(ngm_result)) {
      inc_age_dist <- ngm_result$asym_incdist_age_NGM
      results <- rbind(results, data.frame(
        sample_id = i,
        age_group = seq_along(inc_age_dist),
        value = inc_age_dist,
        reproduction_number = ngm_result$reproduction_number,
        R_gg = ngm_result$R_gg,
        R_gs = ngm_result$R_gs,
        R_sg = ngm_result$R_sg,
        R_ss = ngm_result$R_ss
      ))
    } else {
      results <- rbind(results, data.frame(
        sample_id = i,
        age_group = 1:M,
        value = NA,
        reproduction_number = NA,
        R_gg = NA,
        R_gs = NA,
        R_sg = NA,
        R_ss = NA
      ))
    }
  }

  LHS_results <- na.omit(results)
  return(LHS_results)
}
