#' Neutrosophic Rayleigh Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Rayleigh distribution with
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
#' dnsrayleigh(x = remission, theta = c(9.6432,9.8702))
#'
#' pnsrayleigh(q = 20, theta = c(9.6432,9.8702))
#'
#' # Calculate quantiles
#' qnsrayleigh(p = c(0.25,0.5,0.75), theta = c(9.6432,9.8702))
#'
#' # Simulate 10 values
#' rnsrayleigh(n, theta = c(9.6432,9.8702))
#'
#' @export
dnsrayleigh <- function(x, theta) {
  if (any(theta <= 0))
    stop("Arguments are incompatible.")

  theta <- rep(theta, length.out = 2)
  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- (x[, i] / theta[i]^2) * exp((-1 / 2) * (x[, i] / theta[i])^2)
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Rayleigh
#' @export
pnsrayleigh <- function(q, theta, lower.tail = TRUE) {
  if (any(theta <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  theta <- rep(theta, length.out = 2)
  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- 1 - exp((-1 / 2) * (q / theta)^2)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

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
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- theta[i] * sqrt(-2 * log(1 - p[, i]))
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Rayleigh
#' @export
rnsrayleigh <- function(n, theta) {
  if (any(theta <= 0))
    stop(message = "Arguments are incompatible.")

  theta <- rep(theta, length.out = 2)
  u <- matrix(runif(n * length(theta)), nrow = n, ncol = length(theta))
  X <- qnsrayleigh(u, theta)
  return(X)
}