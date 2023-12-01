#' Neutrosophic Discrete Uniform Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic discrete uniform distribution with
#' parameter \eqn{k_N}.
#'
#' Let \eqn{X_N} be a neutrosophic random variable and denote
#' \eqn{X_N \sim \mathcal{U}(1,2,\ldots,k_N)} as neutrosophic discrete
#' uniform distribution with parameter \eqn{k_N} has the density
#' \deqn{f_N(x)=\frac{1}{k_N}}
#' for \eqn{k_N \in (k_L, k_U)}.
#'
#' @name Neutrosophic Discrete Uniform
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param k parameter of the distribution that must be a positive finite interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsdunif} gives the distribution function,
#'  \code{dnsdunif} gives the density,
#'  \code{qnsdunif} gives the quantile function and
#'  \code{rnsdunif} generates random variables from the neutrosophic Discrete Uniform Distribution.
#' @references
#'        Granados, C. (2022). Some discrete neutrosophic distributions with
#'        neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5),
#'          1442-1457.
#'
#' @examples
#' dnsdunif(x = 8, k = c(10, 11))
#' dnsdunif(x = c(8, 9), k = c(10, 11))
#' @export
dnsdunif <- function(x, k) {
  if (any(k <= 0)) {
    stop("Arguments are incompatible.")
  }


  if (any(round(x) != x)) {
    warning(paste("non-integer"))
  }

  k <- rep(k, length.out = 2)

  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- ifelse(x[, i] >= 1 & x[, i] <= k[i] & round(x[, i]) == x[, i], 1 / k[i], 0)
  }


  return(pdf)
}
#' @name Neutrosophic Discrete Uniform
#' @examples
#'
#' pnsdunif(q = 2, k = c(10, 11))
#' @export
pnsdunif <- function(q, k, lower.tail = TRUE) {
  if (any(k <= 0)) {
    stop("Arguments are incompatible.")
  }


  q <- floor(q)
  k <- rep(k, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }


  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- ifelse(q[, i] < 1, 0, ifelse(q[, i] <= k[i], floor(q[, i]) / k[i], 1))
  }

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  return(cdf)
}
#' @name Neutrosophic Discrete Uniform
#' @examples
#'
#' qnsdunif(p = 0.2, k = c(10, 11))
#'
#' @export
qnsdunif <- function(p, k) {
  if (any(k <= 0)) {
    stop("Arguments are incompatible.")
  }


  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  k <- rep(k, length.out = 2)

  if (is.vector(p)) {
    p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- ceiling(k[i] * p[, i])
  }


  return(quantiles)
}
#' @name Neutrosophic Discrete Uniform
#' @examples
#'
#' # Simulate 10 numbers
#' rnsdunif(n = 10, k = c(10, 11))
#'
#' @export
rnsdunif <- function(n, k) {
  if (any(k < 0) || any(k == 0)) {
    stop("Arguments are incompatible.")
  }

  k <- rep(k, length.out = 2)

  X <- qnsdunif(runif(n), k)

  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
