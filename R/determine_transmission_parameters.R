#' Determine transmission parameters from SAR, age distribution, contact matrix and recovery rate
#'
#' @inheritParams model_mpox_network
#'
#' @return list of transmission parameters
#'  - `p_cont` - transmission probability in global transmission
#'  - `c_ave` - average contact rate
#'  - `beta` - transmission rate sexual contacts
#' @export
determine_transmission_parameters <- function(age_distribution, contact_matrix, gamma) {
  # per-contact probability of infection: community contacts
  ### Age-weighted household contact matrix, and its eigenvalue
  age_matrix <- matrix(age_distribution, nrow = length(age_distribution), ncol = length(age_distribution), byrow = FALSE)
  # elementwise multiplication
  sym_contact_home_matrix_ageWeighted <- (age_matrix * contact_matrix)
  # average household contact rate
  c_ave <- eigen(sym_contact_home_matrix_ageWeighted)$value[1]

  ### SAR to per-contact probability of infection
  SAR_h <- 0.1 # household SAR seq(0,0.2,by=0.01)
  p_cont <- SAR_h / (1 - SAR_h) * (gamma / c_ave)

  # transmission rate sexual contacts: beta
  SAR_s <- 0.8
  beta <- SAR_s * gamma / (1 - SAR_s)
  return(list(
    p_cont = p_cont,
    beta = beta,
    c_ave = c_ave
  ))
}
