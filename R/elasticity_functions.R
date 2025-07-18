#' top eigenvalue of the matrix
#'
#' @param x matrix
#'
#' @return largest real eigenvalue
#' @export
lambda <- function(x) {
  return(Re(RSpectra::eigs(x, k = 1, which = "LM")$values)) # Extracts the largest real eigenvalue
}


#' Sensitivity
#'
#' @param matA original matrix A
#' @param delta_matA matrix with perturbed elements
#'
#' @return matrix
#' @export
#'
#' @examples matA <- matrix(1:16, nrow = 4, ncol = 4)
#' delta_matA <- matrix(1e-06, nrow = 4, ncol = 4)
#' sens_func(matA, delta_matA)
sens_func <- function(matA,
                      delta_matA) {
  # matrix dimension
  m <- nrow(matA)

  # eigenvalue
  stat <- lambda(matA)

  # initialize sensitivity matrix
  sensA <- matrix(NA, ncol = m, nrow = m)

  # matrix perturbation
  for (i in 1:m) {
    for (j in 1:m) {
      pertA <- matA
      pertA[i, j] <- pertA[i, j] + delta_matA[i, j]
      statPert <- lambda(pertA)
      if ((matA[i, j] - pertA[i, j]) == 0) {
        sensA[i, j] <- 0
      } else {
        sensA[i, j] <- (stat - statPert) / (matA[i, j] - pertA[i, j])
      }
    }
  }
  return(sensA)
}


#' Elasticity
#'
#' @param matA original matrix A
#' @param delta_matA matrix with perturbed elements
#'
#' @return matrix
#' @export
#'
#' @examples matA <- matrix(1:16, nrow = 4, ncol = 4)
#' delta_matA <- matrix(1e-06, nrow = 4, ncol = 4)
#' elas_func(matA, delta_matA) # elasticity
elas_func <- function(matA,
                      delta_matA) {
  # eigenvalue
  stat <- lambda(matA)

  # sensitivity matrix
  sensA <- sens_func(matA, delta_matA)

  # return
  return(sensA * matA / stat)
}
