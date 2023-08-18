#' Neutrosophic Weibull Distribution (NWD)
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
#' @name NWD
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
#'  \code{pnweibull} gives the distribution function,
#'  \code{dnweibull} gives the density,
#'  \code{qnweibull} gives the quantile function and
#'  \code{rnweibull} generates random variables from the neutrosophic Weibull dDistribution.
#' @references
#'    Alhasan, K. F. H. and Smarandache, F. (2019). Neutrosophic Weibull
#'    distribution and Neutrosophic Family Weibull Distribution,
#'    \emph{Neutrosophic Sets and Systems}, 28, 191-199.
#'
#' @importFrom stats runif dweibull pweibull qweibull
#' @examples
#' dnweibull(x, shape = 1, scale = 2)
#'
#' dnweibull(x2, shape = c(1, 2), scale = c(2, 2))
#' @export
dnweibull <- function(x, shape, scale, log = FALSE) {
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
    pdf[, i] <- stats::dweibull(x[, i], shape = shape[i], scale = scale[i], log = log)
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnweibull(x, shape = 1, scale = 2)
#' pnweibull(x2, shape = c(1, 2), scale = c(2, 2))
#' @export
pnweibull <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(scale <= 0) || any(shape <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  scale <- rep(scale, length.out = 2)
  shape  <- rep(shape, length.out = 2)

  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pweibull(q, shape = shape, scale = scale, lower.tail = lower.tail, log.p =log.p)


  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
  }
#' @name NWD
#' @name NWD
#' @examples
#' qnweibull(x, shape = 1, scale = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnweibull(x2, shape = c(1, 2), scale = c(2, 2))
#' @export
qnweibull <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE) {
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
    quantiles[, i] <- stats::qweibull(p[, i], shape = shape[i], scale = scale[i], lower.tail = lower.tail, log.p =log.p)
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}
#' @name NWD
#' @examples
#' n <- 10
#' rnweibull(n, shape = 1, scale = 2)
#' rnweibull(n, shape = c(1, 2), scale = c(1, 1))
#' @export
rnweibull <- function(n, shape, scale) {
  if (any(scale <= 0) || any(shape <= 0))
    stop(message = "Arguments are incompatible.")

  scale <- rep(scale, length.out = 2)
  shape  <- rep(shape, length.out = 2)

  u <- matrix(runif(n), ncol = 2)
  X <- qnweibull(u, scale, shape)

  return(X)
}
