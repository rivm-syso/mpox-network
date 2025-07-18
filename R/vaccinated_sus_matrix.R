#' Function to compute effective vaccination coverage for age group i with degree n
#'
#' @param M int - total number of age groups
#' @param N_1 int - first sexually active age group
#' @param N_2 int - last sexually active age group, 1 <= N_1 <= N_2 <= M
#' @param effV_mat matrix[n_active, D_max] - effective vaccine coverage for group (i,n)
#' @param effV_age vector[M] - effective vaccine coverage for all age groups
#' @param S vector[M] - proportion susceptible for all age groups
#' @param xS array[n_active, n_active, D_max, D_max] - probability that a binding site is with group (i,n) and having a susceptible partner of group (j,k). Array index = [i, j, n, k]
#' @param xI array[n_active, n_active, D_max] - probability that a binding site is with group (i,n) and having an infectious partner of group j. Array index = [i, j, n]
#' @param xbar array[n_active, n_active, D_max] - probability that a binding site is with group (i,n) and having a partner of group j. Array index = [i, j, n]
#'
#' @return list with four elements
#' - xS_vac array[n_active, n_active, D_max, D_max] - xS after vaccination
#' - xI_vac array[n_active, n_active, D_max] - xI after vaccination
#' - xbar_vac array[n_active, n_active, D_max] - xbar after vaccination
#' - S_vac vector[M] - susceptible proportion after vaccination for all age groups
#' @export
vaccinated_sus_matrix <- function(
    M,
    N_1,
    N_2,
    effV_mat,
    effV_age,
    S,
    xS,
    xI,
    xbar) {
  # sexually active age groups and max degree
  n_active <- N_2 - N_1 + 1
  D_max <- length(effV_mat[1, ])

  # check if dimension is correct
  stopifnot("Dimension of effV_mat does not match the number of sexually active age groups" = n_active == length(effV_mat[, 1]))

  ### sexual transmission: susceptible proportion after vaccination
  # Prepare arrays to store the results
  xS_vac <- array(0, dim = dim(xS))
  xI_vac <- array(0, dim = dim(xI))
  xbar_vac <- array(0, dim = dim(xbar))

  # multiply effV_val[i,n] for all j and m
  for (i in 1:n_active) {
    for (n in 1:D_max) {
      # effective vaccination coverage for group (i,n)
      effV_val <- effV_mat[i, n]

      # change in susceptible proportions for group (i,n)
      xS_vac[i, , n, ] <- xS[i, , n, ] * (1 - effV_val)
      xI_vac[i, , n] <- xI[i, , n] * (1 - effV_val)
      xbar_vac[i, , n] <- xbar[i, , n] * (1 - effV_val)
    }
  }

  ### global transmission: susceptible proportion after vaccination
  S_vac <- S * (1 - effV_age)

  # return results as a list
  return(
    list(
      xS_vac = xS_vac,
      xI_vac = xI_vac,
      xbar_vac = xbar_vac,
      S_vac = S_vac
    )
  )
}
