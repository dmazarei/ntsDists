#' Neutrosophic Uniform Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Uniform distribution of a continuous
#' variable \eqn{X} with parameters \eqn{a_N} and  \eqn{b_N}.
#'
#' The neutrosophic Uniform distribution with parameters
#' \eqn{\max_N} and \eqn{\min_N} has the density
#' \deqn{f_N(x)=\frac{1}{b_N-a_N}}
#' for \eqn{a_N \in (a_L, a_U)}  lower parameter interval, \eqn{b_N \in (b_L,b_U)},
#'  upper parameter interval.
#'
#' @name Neutrosophic Uniform
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param min lower limits of the distribution. Must be finite.
#' @param max upper limits of the distribution. Must be finite.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsunif} gives the distribution function,
#'  \code{dnsunif} gives the density,
#'  \code{qnsunif} gives the quantile function and
#'  \code{rnsunif} generates random variables from the neutrosophic Uniform Distribution.
#' @references
#'    Alhabib, R., Ranna, M. M., Farah, H., & Salama, A. A. (2018).
#'     Some neutrosophic probability distributions.
#'     \emph{Neutrosophic Sets and Systems}, 22, 30-38.
#'
#' @importFrom stats runif dunif punif qunif
#' @examples
#'
#' dnsunif(x = 1, min = c(0,5), max = c(15,20))
#' dnsunif(x = c(6,10), min = c(0,5), max = c(15,20))
#'
#' punif(x = 1, min = c(0,5), max = c(15,20))
#' punif(x = c(6,10), min = c(0,5), max = c(15,20))
#'
#' qnsunif(x, min = c(0,5), max = c(15,20))
#'
#' n <- 10
#' rnsunif(n, min = c(0,5), max = c(15,20))
#'
#' @export
dnsunif <- function(x, min, max) {
  if (any(max <= min)) {
    stop(message = "Arguments are incompatible.")
  }
  max <- rep(max, length.out = 2)
  min <- rep(min, length.out = 2)

  if (is.vector(x)) {
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dunif(x[, i], min = min[i], max = max[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Uniform
#' @export
pnsunif <- function(q, min, max, lower.tail = TRUE) {
  if (any(max <= min)) {
    stop(message = "Arguments are incompatible.")
  }

  max <- rep(max, length.out = 2)
  min <- rep(min, length.out = 2)

  if (is.vector(q)) {
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::punif(q, min = min, max = max)

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name Neutrosophic Uniform
#' @export
qnsunif <- function(p, min, max) {
  if (any(max <= min)) {
    stop(message = "Arguments are incompatible.")
  }
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  max <- rep(max, length.out = 2)
  min <- rep(min, length.out = 2)

  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qunif(p[, i], min = min[i], max = max[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}
#' @name Neutrosophic Uniform
#' @export
rnsunif <- function(n, min, max) {
  if (any(max <= min)) {
    stop(message = "Arguments are incompatible.")
  }
  max <- rep(max, length.out = 2)
  min <- rep(min, length.out = 2)

  u <- matrix(runif(n * 2), nrow = n, ncol = 2)
  X <- qnsunif(u, max, min)

  return(X)
}
