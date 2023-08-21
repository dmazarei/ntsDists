#' Neutrosophic Weibull Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Weibull distribution with \code{scale}
#' parameter \eqn{\alpha_N} and \code{shape} parameter \eqn{\beta_N}.
#'
#' The neutrosophic Rayleigh distribution with parameters
#' \eqn{\alpha_N} and \eqn{\beta_N} has the density
#' \deqn{f_N(x)=\frac{\beta_N}{\alpha_N^{\beta_N}} x^{\beta_N-1}
#'     \exp\{-\left(x / \alpha_N\right)^{\beta_N}\}}
#' for \eqn{\beta_N \in (\beta_L, \beta_U)} the shape parameter must
#' be a positive interval, \eqn{\alpha_N \in (\alpha_L,\alpha_U)},
#' the scale parameter which be a positive interval, and \eqn{x > 0}.
#'
#' @name Neutrosophic Weibull
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param shape shape parameter, which must be a positive interval.
#' @param scale scale parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsweibull} gives the distribution function,
#'  \code{dnsweibull} gives the density,
#'  \code{qnsweibull} gives the quantile function and
#'  \code{rnsweibull} generates random variables from the neutrosophic Weibull dDistribution.
#' @references
#'    Alhasan, K. F. H. and Smarandache, F. (2019). Neutrosophic Weibull
#'    distribution and Neutrosophic Family Weibull Distribution,
#'    \emph{Neutrosophic Sets and Systems}, 28, 191-199.
#'
#' @importFrom stats runif dweibull pweibull qweibull
#' @examples
#' data(remission)
#' dnsweibull(x = remission, shape = c(1.0519,1.0553), scale = c(9.3370,9.4544))
#'
#' pnsweibull(q = 20, shape = c(1.0519,1.0553), scale = c(9.3370,9.4544))
#'
#' # Calculate quantiles
#' qnsweibull(p = c(0.25,0.5,0.75), shape = c(1.0519,1.0553), scale = c(9.3370,9.4544))
#'
#' # Simulate 10 numbers
#' rnsweibull(n = 10, shape = c(1.0519,1.0553), scale = c(9.3370,9.4544))
#'
#' @export
dnsweibull <- function(x, shape, scale) {
  if (any(scale <= 0) || any(shape <= 0) || any(x < 0))
    stop("Arguments are incompatible.")

  scale <- rep(scale, length.out = 2)
  shape  <- rep(shape, length.out = 2)

  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dweibull(x[, i], shape = shape[i], scale = scale[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Weibull
#' @export
pnsweibull <- function(q, shape, scale, lower.tail = TRUE) {
  if (any(scale <= 0) || any(shape <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  scale <- rep(scale, length.out = 2)
  shape  <- rep(shape, length.out = 2)

  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pweibull(q, shape = shape, scale = scale)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
  }
#' @name Neutrosophic Weibull
#' @export
qnsweibull <- function(p, shape, scale) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(scale <= 0) || any(shape <= 0)){
    stop(message = "Arguments are incompatible.")
  }

  scale <- rep(scale, length.out = 2)
  shape  <- rep(shape, length.out = 2)

  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qweibull(p[, i], shape = shape[i], scale = scale[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}
#' @name Neutrosophic Weibull
#' @export
rnsweibull <- function(n, shape, scale) {
  if (any(scale <= 0) || any(shape <= 0))
    stop(message = "Arguments are incompatible.")

  scale <- rep(scale, length.out = 2)
  shape  <- rep(shape, length.out = 2)

  u <- matrix(runif(n), ncol = 2)
  X <- qnsweibull(u, scale, shape)

  return(X)
}
