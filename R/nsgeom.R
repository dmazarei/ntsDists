#' Neutrosophic Geometric Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Geometric distribution with
#' parameter \code{prob} = \eqn{p_N}.
#'
#' The neutrosophic Geometric distribution with parameter \eqn{p_N}
#' has the density
#' \deqn{f_X(x)=p_N\left(1-p_N\right)^x}
#' for \eqn{p_N \in (p_L, p_U)} which must be \eqn{0<p_N<1}
#' and \eqn{x \in \{0, 1, 2, \ldots\}}.
#'
#' @name Neutrosophic Geometric
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param prob probability of success on each trial, \eqn{0 < prob <= 1}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsgeom} gives the distribution function,
#'  \code{dnsgeom} gives the density,
#'  \code{qnsgeom} gives the quantile function and
#'  \code{rnsgeom} generates random variables from the Geometric Distribution.
#' @references
#'        Granados, C. (2022).
#'        Some discrete neutrosophic distributions with neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5), 1442-1457.
#' @importFrom stats runif dgeom pgeom qgeom
#' @examples
#' # One person participates each week with a ticket in a lottery game, where
#' # the probability of winning the first prize is (10^(-8), 10^(-6)).
#' # Probability of one persons wins at the fifth year?
#' dnsgeom(x = 5, prob = c(1e-8,1e-6))
#' @export
dnsgeom <- function(x, prob) {
  if (any(prob <= 0) || any(prob > 1) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
  }

  if (any(x < 0) && any(x - floor(x) == 0)) {
    stop(message = "Warning: x should be a positive integer.")
  }

  prob <- rep(prob, length.out = 2)

  if (is.vector(x)) {
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dgeom(x[, i], prob = prob[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Geometric
#' @examples
#' # Probability of one persons wins after 10 years?
#' pnsgeom(x = 10, prob = c(1e-8,1e-6), lower.tail = FALSE)
#' @export

pnsgeom <- function(q, prob, lower.tail = TRUE) {
  if (any(prob <= 0) || any(prob > 1) || any(q < 0)) {
    stop(message = "Arguments are incompatible.")
  }
  if (any(q < 0) && any(q - floor(q) == 0)) {
    stop(message = "Warning: q should be a  positive integer.")
  }

  prob <- rep(prob, length.out = 2)
  if (is.vector(q)) {
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pgeom(q, prob = prob)

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name Neutrosophic Geometric
#' @examples
#' # Calculate the quantiles
#' qnsgeom(q = c(0.25,0.5,0.75), prob = c(1e-8,1e-6))
#' @export
qnsgeom <- function(p, prob) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  prob <- rep(prob, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qgeom(p[, i], prob = prob[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Geometric
#' @examples
#' # Simulate 10 numbers
#' rnsgeom(n = 10, prob = 0.5)
#' @export
rnsgeom <- function(n, prob) {
  if (any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  prob <- rep(prob, length.out = 2)
  u <- matrix(runif(n), ncol = 2)
  X <- qnsgeom(u, prob)

  return(X)
}
