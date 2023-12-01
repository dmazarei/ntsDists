#' Neutrosophic Generalized Exponential Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic generalized exponential
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
#' dnsgexp(x = remission, nu = c(7.9506, 8.0568), delta = c(1.2390, 1.2397))
#'
#' @export
dnsgexp <- function(x, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
  }

  nu <- rep(nu, length.out = 2)
  delta <- rep(delta, length.out = 2)


  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }


  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- (delta[i] / nu[i]) * (1 - exp(-x[, i] / nu[i]))^(delta[i] - 1) * exp(-x[, i] / nu[i])
  }


  return(pdf)
}

#' @name Neutrosophic Generalized Exponential
#' @examples
#' pnsgexp(q = 20, nu = c(7.9506, 8.0568), delta = c(1.2390, 1.2397))
#'
#' @export
pnsgexp <- function(q, nu, delta, lower.tail = TRUE) {
  if (any(nu <= 0) || any(delta <= 0) || any(q < 0)) {
    stop(message = "incompatible arguments.")
  }

  nu <- rep(nu, length.out = 2)
  delta <- rep(delta, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }



  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- (1 - exp(-q[, i] / nu[i]))^delta[i]
  }


  if (!lower.tail) {
    cdf <- 1 - cdf
  }


  return(cdf)
}

#' @name Neutrosophic Generalized Exponential
#' @examples
#'
#' # Calcluate quantiles
#' qnsgexp(c(0.25, 0.5, 0.75), nu = c(7.9506, 8.0568), delta = c(1.2390, 1.2397))
#'
#' @export
qnsgexp <- function(p, nu, delta) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(nu <= 0) || any(delta <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  nu <- rep(nu, length.out = 2)
  delta <- rep(delta, length.out = 2)

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- log(-p[, i]^(1 / delta[i]) + 1) * (-nu[i])
  }

  return(quantiles)
}

#' @name Neutrosophic Generalized Exponential
#' @examples
#' # Simulate 10 values
#' rnsgexp(n = 10, nu = c(7.9506, 8.0568), delta = c(1.2390, 1.2397))
#'
#' @export
rnsgexp <- function(n, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  nu <- rep(nu, length.out = 2)
  delta <- rep(delta, length.out = 2)

  X <- qnsgexp(runif(n), nu, delta)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]


  return(X)
}
