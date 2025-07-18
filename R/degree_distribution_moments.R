#' First and second moment of a degree distribution
#'
#' @param degree_distribution vector of probabilities, where the index represents
#' the degree
#'
#' @return list of `first_moment` and `second_moment` of the degree distribution
#' @export
degree_distribution_moments <- function(degree_distribution) {
  D_max <- length(degree_distribution)
  first_moment <- sum(1:D_max * degree_distribution)
  second_moment <- sum((1:D_max)^2 * degree_distribution)

  return(list(
    first_moment = first_moment,
    second_moment = second_moment
  ))
}
