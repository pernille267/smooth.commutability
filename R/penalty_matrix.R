#' O'Sullivan penalized splines penalty matrix
#'
#' @param x_min A \code{double} indicating the lower boundary knot. The knot sequence for the B-spline basis will be augmented using this value.
#' @param x_max A \code{double} indicating the upper boundary knot. The knot sequence for the B-spline basis will be augmented using this value.
#' @param interior_knots A \code{numeric} vector containing the interior knots for the B-spline basis. Cannot include \code{x_min} or \code{x_max}.
#'
#' @return The O'Sullivan B-spline penalty matrix
#' @export
#'
#' @examples print(1)

penalty_matrix <- function(x_min, x_max, interior_knots){
  all_knots <- c(rep(x_min, 4L), interior_knots, rep(x_max, 4L))
  n_interior_knots <- length(interior_knots)
  n_rows <- 3 * (n_interior_knots + 8)
  xtilde <- (rep(all_knots, each = 3)[-c(1, (n_rows-1), n_rows)] + rep(all_knots, each = 3)[-c(1, 2, n_rows)]) / 2
  wts <- rep(diff(all_knots), each = 3) * rep(c(1, 4, 1) / 6, n_interior_knots + 7)
  Bdd <- spline.des(knots = all_knots, x = xtilde, derivs = rep(2, length(xtilde)), outer.ok = TRUE)$design
  Omega <- t(Bdd) %*% diag(wts) %*% Bdd
  return(Omega)
}