#' Neutrosophic Uniform Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic Uniform distribution of a continuous
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
#' dnsunif(x = 1, min = c(0, 5), max = c(15, 20))
#' dnsunif(x = c(6, 10), min = c(0, 5), max = c(15, 20))
#'
#' punif(q = 1, min = c(0, 5), max = c(15, 20))
#' punif(q = c(6, 10), min = c(0, 5), max = c(15, 20))
#'
#' qnsunif(p = c(0.25, 0.5, 0.75), min = c(0, 5), max = c(15, 20))
#'
#' rnsunif(n = 10, min = c(0, 5), max = c(15, 20))
#'
#' @export
dnsunif <- function(x, min, max) {
  if (any(max <= min)) {
    stop(message = "Arguments are incompatible.")
  }
  max <- rep(max, length.out = 2)
  min <- rep(min, length.out = 2)

  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- stats::dunif(x[, i], min = min[i], max = max[i])
  }



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

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- stats::punif(q[, i], min = min[i], max = max[i])
  }



  if (!lower.tail) {
    cdf <- 1 - cdf
  }


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

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- stats::qunif(p[, i], min = min[i], max = max[i])
  }
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

  X <- qnsunif(runif(n), min, max)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
