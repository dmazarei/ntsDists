#' Neutrosophic Generalized Exponential Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic generalized exponential
#' distribution with shape parameter \eqn{\delta_N} and scale parameter
#' \eqn{\nu_N} or equ.
#'
#' The neutrosophic generalized exponential distribution with parameters
#' \eqn{\delta} and \eqn{\nu} has density
#' \deqn{f_n(x)=\frac{\delta_N}{\nu_N}\left(1-\exp \left\{-\frac{x_N}{\nu_N}\right\}\right)^{\delta_N-1} e^{\left\{-\frac{x_N}{\nu_N}\right\}}}
#' for \eqn{\delta_N \in (\delta_L, \delta_U)}, the shape parameter
#' which must be a positive interval, and \eqn{\nu_N \in (\nu_L, \nu_U)}, the
#' scale parameter which must also be a positive interval, and \eqn{x \ge 0}.
#'
#' @name Neutrosophic Generalized Exponential
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param nu the scale parameter, which must be a positive interval.
#' @param delta  the shape parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsgexp} gives the distribution function,
#'  \code{dnsgexp} gives the density,
#'  \code{qnsgexp} gives the quantile function and
#'  \code{rnsgexp} generates random variables
#'  from the neutrosophic generalized exponential distribution.
#'
#' @references
#'    Rao, G. S., Norouzirad, M., and Mazarei . D. (2023). Neutrosophic
#'    Generalized Exponential Distribution with Application.
#'    \emph{Neutrosophic Sets and Systems}, 55, 471-485.
#'
#' @importFrom stats runif
#' @examples
#'
#' data(remission)
#' dnsgexp(x = remission, nu = c(7.9506,8.0568), delta = c(1.2390,1.2397))
#'
#' @export
dnsgexp <- function(x, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0) || any(x < 0))
    stop(message = "Arguments are incompatible.")

  nu     <- rep(nu, length.out = 2)
  delta  <- rep(delta, length.out = 2)


  if (is.vector(x) && length(x) == 1) {
    x <- matrix(rep(x, each = 2), ncol = 2, byrow = TRUE)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- (delta[i] / nu[i]) * (1 - exp(-x[, i] / nu[i]))^(delta[i] - 1) * exp(-x[, i] / nu[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}

#' @name Neutrosophic Generalized Exponential
#' @examples
#' pnsgexp(q = 20, nu = c(7.9506,8.0568), delta = c(1.2390,1.2397))
#'
#' @export
pnsgexp <- function(q, nu, delta, lower.tail = TRUE) {
  if (any(nu <= 0) || any(delta <= 0) || any(q < 0))
    stop(message = "incompatible arguments.")

  nu     <- rep(nu, length.out = 2)
  delta  <- rep(delta, length.out = 2)

  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- (1 - exp(-q / nu))^delta

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)
  # Identify rows where col1 > col2
  swap_rows <- cdf[, 1] > cdf[, 2]
  # Swap values using logical indexing
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}

#' @name Neutrosophic Generalized Exponential
#' @examples
#'
#' # Calcluate quantiles
#' qnsgexp(c(0.25,0.5,0.75), nu = c(7.9506,8.0568), delta = c(1.2390,1.2397))
#'
#' @export
qnsgexp <- function(p, nu, delta) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(nu <= 0) || any(delta <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  nu     <- rep(nu, length.out = 2)
  delta  <- rep(delta, length.out = 2)

  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- log(-p[, i]^(1 / delta[i]) + 1) * (-nu[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Generalized Exponential
#' @examples
#' # Simulate 10 values
#' rnsgexp(n = 10, nu = c(7.9506,8.0568), delta = c(1.2390,1.2397))
#'
#' @export
rnsgexp <- function(n, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0))
    stop(message = "Arguments are incompatible.")

  nu     <- rep(nu, length.out = 2)
  delta  <- rep(delta, length.out = 2)

  u <- matrix(runif(n), ncol = 2)
  X <- qnsgexp(u, nu, delta)


  return(X)
}
