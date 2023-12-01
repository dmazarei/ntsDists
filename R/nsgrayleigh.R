#' Neutrosophic Generalized Rayleigh Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic generalized Rayleigh distribution with
#' parameters \code{shape} = \eqn{\nu_N} and \code{scale} = \eqn{\sigma_N}.
#'
#' The neutrosophic generalized Rayleigh distribution with parameters \eqn{\nu_N} and
#' \eqn{\sigma_N} has the density
#' \deqn{f_N(x)=\frac{2\nu_N}{\sigma_N^2}x \exp\{-\left(\frac{x}{\sigma_N} \right)^2\}\left[1-\exp\{-\left(\frac{x}{\sigma_N} \right)^2\}\right]^{\nu_N-1}}
#' for \eqn{x > 0}, \eqn{\nu_N \in (\nu_L, \nu_U)}, the shape
#' parameter which must be a positive interval and
#' \eqn{\sigma_n \in (\sigma_L, \sigma_U)}, the scale parameter which
#' must be a positive interval.
#'
#' @name Neutrosophic Generalized Rayleigh
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param shape the shape parameter, which must be a positive interval.
#' @param scale the scale parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \leq x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{dnsgrayleigh} gives the density,
#'  \code{pnsgrayleigh} gives the distribution function,
#'  \code{qnsgrayleigh} gives the quantile function and
#'  \code{rnsgrayleigh} generates random variables from the Neutrosophic Generalized Rayleigh Distribution.
#' @references
#' Norouzirad, M., Rao, G. S., & Mazarei, D. (2023).
#' Neutrosophic Generalized Rayleigh Distribution with Application.
#' \emph{Neutrosophic Sets and Systems}, 58(1), 250-262.
#'
#' @importFrom stats runif
#' @examples
#' data(remission)
#' dnsgrayleigh(x = remission,shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' pnsgrayleigh(q = 20, shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' # Calculate quantiles
#' qnsgrayleigh(p = c(0.25, 0.5, 0.75), shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' # Simulate 10 values
#' rnsgrayleigh(n = 10, shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' @export
dnsgrayleigh <- function(x, shape, scale) {
  if (any(shape <= 0) || any(scale <= 0) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }


  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- (2 * shape[i]) / (scale[i]^2) * x[,i] * exp(-(x[,i] / scale[i])^2) * (1 - exp(-(x[,i] / scale[i])^2))^(shape[i] - 1)
  }



  return(pdf)
}
#' @name Neutrosophic Generalized Rayleigh
#' @export
pnsgrayleigh <- function(q, shape, scale, lower.tail = TRUE) {
  if (any(shape <= 0) || any(scale <= 0) || any(q < 0)) {
    stop(message = "incompatible arguments.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }



  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- (1 - exp(-(q[,i] / scale[i])^2))^(shape[i])
  }

  if (!lower.tail) {
    cdf <- 1 - cdf
  }


  return(cdf)
}
#' @name Neutrosophic Generalized Rayleigh
#' @export
qnsgrayleigh <- function(p, shape, scale) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(shape <= 0) || any(scale <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- sqrt(-scale[i]^2 * log(1 - p[,i]^(1 / shape[i])))
  }

  return(quantiles)
}

#' @name Neutrosophic Generalized Rayleigh
#' @export
rnsgrayleigh <- function(n, shape, scale) {
  if (any(shape <= 0) || any(scale <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  X <- qnsgpd(runif(n), shape, scale)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
