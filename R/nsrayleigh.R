#' Neutrosophic Rayleigh Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic Rayleigh distribution with
#' parameter \eqn{\theta_N}.
#'
#' The neutrosophic Rayleigh distribution with parameter \eqn{\theta}
#' has the density
#' \deqn{f_N(x)=\frac{x}{\theta_N^2} e^{-\frac{1}{2}\left(\frac{x}{\theta_N}\right)^2}}
#' for  \eqn{\theta_N \in (\theta_L, \theta_U)}, which must be a positive
#' interval and \eqn{x \ge 0}.
#'
#' @name Neutrosophic Rayleigh
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param theta the shape parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{dnsrayleigh} gives the density,
#'  \code{pnsrayleigh} gives the distribution function,
#'  \code{qnsrayleigh} gives the quantile function and
#'  \code{rnsrayleigh} generates random variables from the Neutrosophic Rayleigh Distribution.
#' @references
#' Khan, Z., Gulistan, M., Kausar, N. and Park, C. (2021).
#' Neutrosophic Rayleigh Model With Some Basic Characteristics and
#' Engineering Applications, in \emph{IEEE Access}, 9, 71277-71283.
#'
#' @importFrom stats runif
#' @examples
#' data(remission)
#' dnsrayleigh(x = remission, theta = c(9.6432, 9.8702))
#'
#' pnsrayleigh(q = 20, theta = c(9.6432, 9.8702))
#'
#' # Calculate quantiles
#' qnsrayleigh(p = c(0.25, 0.5, 0.75), theta = c(9.6432, 9.8702))
#'
#' # Simulate 10 values
#' rnsrayleigh(n = 10, theta = c(9.6432, 9.8702))
#'
#' @export
dnsrayleigh <- function(x, theta) {
  if (any(theta <= 0)) {
    stop("Arguments are incompatible.")
  }

  theta <- rep(theta, length.out = 2)
  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }


  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- (x[, i] / theta[i]^2) * exp((-1 / 2) * (x[, i] / theta[i])^2)
  }



  return(pdf)
}
#' @name Neutrosophic Rayleigh
#' @export
pnsrayleigh <- function(q, theta, lower.tail = TRUE) {
  if (any(theta <= 0) || any(q < 0)) {
    stop("Arguments are incompatible.")
  }

  theta <- rep(theta, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }



  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- 1 - exp((-1 / 2) * (q[, i] / theta[i])^2)
  }



  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  return(cdf)
}
#' @name Neutrosophic Rayleigh
#' @export
qnsrayleigh <- function(p, theta) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(theta <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  theta <- rep(theta, length.out = 2)

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- theta[i] * sqrt(-2 * log(1 - p[, i]))
  }


  return(quantiles)
}

#' @name Neutrosophic Rayleigh
#' @export
rnsrayleigh <- function(n, theta) {
  if (any(theta <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  theta <- rep(theta, length.out = 2)

  X <- qnsrayleigh(runif(n), theta)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
