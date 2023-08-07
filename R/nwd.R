#' Neutrosophic Weibull Distribution (NWD)
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Weibull distribution with scale
#' parameter \eqn{\alpha_N} and shape parameter \eqn{\beta_N}.
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
#' @param beta shape parameter, which must be a positive interval.
#' @param alpha scale parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnwd} gives the distribution function,
#'  \code{dnwd} gives the density,
#'  \code{qnwd} gives the quantile function and
#'  \code{rnwd} generates random variables from the neutrosophic Weibull dDistribution.
#' @references
#'    Alhasan, K. F. H. and Smarandache, F. (2019). Neutrosophic Weibull
#'    distribution and Neutrosophic Family Weibull Distribution,
#'    \emph{Neutrosophic Sets and Systems}, 28, 191-199.
#'
#' @importFrom stats runif dweibull pweibull qweibull
#' @examples
#' dnwd(x, beta = 1, alpha = 2)
#'
#' dnwd(x2, beta = c(1, 2), alpha = c(2, 2))
#' @export
dnwd <- function(x, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0) || any(x < 0))
    stop("Arguments are incompatible.")

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)

  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dweibull(x[, i], shape = beta[i], scale = alpha[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnwd(x, beta = 1, alpha = 2)
#' pnwd(x2, beta = c(1, 2), alpha = c(2, 2))
#' @export
pnwd <- function(q, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)

  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pweibull(q, shape = beta, scale = alpha)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
  }
#' @name NWD
#' @name NWD
#' @examples
#' qnwd(x, beta = 1, alpha = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnwd(x2, beta = c(1, 2), alpha = c(2, 2))
#' @export
qnwd <- function(p, beta, alpha) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(alpha <= 0) || any(beta <= 0)){
    stop(message = "Arguments are incompatible.")
  }

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)

  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qweibull(p[, i], shape = beta[i], scale = alpha[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}
#' @name NWD
#' @examples
#' n <- 10
#' rnwd(n, beta = 1, alpha = 2)
#' rnwd(n, beta = c(1, 2), alpha = c(1, 1))
#' @export
rnwd <- function(n, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0))
    stop(message = "Arguments are incompatible.")

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)

  u <- matrix(runif(n * length(theta)), nrow = n, ncol = length(theta))
  X <- qnwd(u, alpha, beta)

  return(X)
}
