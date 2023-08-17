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
#' @name NGEOMD
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param prob probability of success on each trial, \eqn{0 < prob <= 1}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pngeom} gives the distribution function,
#'  \code{dngeom} gives the density,
#'  \code{qngeom} gives the quantile function and
#'  \code{rngeom} generates random variables from the Geometric Distribution.
#' @references
#'        Granados, C. (2022).
#'        Some discrete neutrosophic distributions with neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5), 1442-1457.
#' @importFrom stats runif dgeom pgeom qgeom
#' @examples
#' dngeom(x, prob = 0.5)
#' dngeom(x2, lambda = c(2, 2))
#' @export
dngeom <- function(x, prob, log = FALSE) {
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
    pdf[, i] <- stats::dgeon(x[, i], prob = prob[i], log = log)
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name NGEOMD
#' @examples
#' x <- 1:10
#' x2 <- matrix(1:20, ncol = 2)
#' pngeom(x, prob = 0.5)
#' pngeom(x2, prob = c(.3, .6))
#' @export

pngeom <- function(q, prob, lower.tail = TRUE, log.p = FALSE) {
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

  cdf <- stats::pgeom(q, prob = prob, log.p = log.p)

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name NGEOMD
#' @examples
#' q1 <- seq(0.1, 1, length.out = 40)
#' qngeom(q1, prob = 0.5)
#' q2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qngeom(q2, lambda = c(2, 2))
#' @export
qngeom <- function(p, prob, log.p = FALSE) {
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
    quantiles[, i] <- stats::qgeom(p[, i], prob = prob[i], log.p = log.p)
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name NGEOMD
#' @examples
#' n <- 10
#' rngeom(n, prob = 0.5)
#' rngeom(n, lambda = c(1, 2))
#' @export
rngeom <- function(n, prob) {
  if (any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  prob <- rep(prob, length.out = 2)
  u <- matrix(runif(n), ncol = 2)
  X <- qngeom(u, prob)

  return(X)
}
