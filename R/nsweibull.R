#' Neutrosophic Weibull Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic Weibull distribution with \code{scale}
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
#' dnsweibull(x = remission, shape = c(1.0519, 1.0553), scale = c(9.3370, 9.4544))
#'
#' pnsweibull(q = 20, shape = c(1.0519, 1.0553), scale = c(9.3370, 9.4544))
#'
#' # Calculate quantiles
#' qnsweibull(p = c(0.25, 0.5, 0.75), shape = c(1.0519, 1.0553), scale = c(9.3370, 9.4544))
#'
#' # Simulate 10 numbers
#' rnsweibull(n = 10, shape = c(1.0519, 1.0553), scale = c(9.3370, 9.4544))
#'
#' @export
dnsweibull <- function(x, shape, scale) {
  if (any(scale <= 0) || any(shape <= 0) || any(x < 0)) {
    stop("Arguments are incompatible.")
  }

  scale <- rep(scale, length.out = 2)
  shape <- rep(shape, length.out = 2)

  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- stats::dweibull(x[, i], shape = shape[i], scale = scale[i])
  }

  return(pdf)
}
#' @name Neutrosophic Weibull
#' @export
pnsweibull <- function(q, shape, scale, lower.tail = TRUE) {
  if (any(scale <= 0) || any(shape <= 0) || any(q < 0)) {
    stop("Arguments are incompatible.")
  }

  scale <- rep(scale, length.out = 2)
  shape <- rep(shape, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- stats::pweibull(q[, i], shape = shape[i], scale = scale[i])
  }



  if (!lower.tail) {
    cdf <- 1 - cdf
  }


  return(cdf)
}
#' @name Neutrosophic Weibull
#' @export
qnsweibull <- function(p, shape, scale) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(scale <= 0) || any(shape <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  scale <- rep(scale, length.out = 2)
  shape <- rep(shape, length.out = 2)

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- stats::qweibull(p[, i], shape = shape[i], scale = scale[i])
  }

  return(quantiles)
}
#' @name Neutrosophic Weibull
#' @export
rnsweibull <- function(n, shape, scale) {
  if (any(scale <= 0) || any(shape <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  scale <- rep(scale, length.out = 2)
  shape <- rep(shape, length.out = 2)

  X <- qnsweibull(runif(n), scale, shape)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
