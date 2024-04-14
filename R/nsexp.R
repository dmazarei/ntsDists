#' Neutrosophic Exponential Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic exponential distribution with the
#' parameter \code{rate} =  \eqn{\theta_N}.
#'
#' The neutrosophic exponential distribution with parameter \eqn{\theta_N}
#' has density
#' \deqn{f_N(x)=\theta_N \exp \left(-x \theta_N\right)}
#' for \eqn{x \ge 0} and \eqn{\theta_N \in (\theta_L, \theta_U)},
#' the rate parameter must be a positive interval and \eqn{x \ge 0}.
#' @name Neutrosophic Exponential
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param rate the shape parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \leq x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{dnsExp} gives the density function
#'
#'  \code{pnsExp} gives the distribution function
#'
#'  \code{qnsExp} gives the quantile function
#'
#'  \code{rnsExp} generates random values from the neutrosophic exponential distribution.
#'
#' @references
#' Duan, W., Q., Khan, Z., Gulistan, M., Khurshid, A. (2021). Neutrosophic
#' Exponential Distribution: Modeling and Applications for Complex Data Analysis,
#' \emph{Complexity}, 2021, 1-8.
#'
#' @importFrom stats runif
#'
#' @examples
#' # Example 4 of Duan et al. (2021)
#'
#' data <- matrix(c(4, 4, 3.5, 3.5, 3.9, 4.1, 4.2, 4.2, 4.3, 4.6, 4.7, 4.7),
#'  nrow = 6, ncol = 2, byrow = TRUE)
#'
#' dnsExp(data, rate = c(0.23, 0.24))
#' dnsExp(x = c(4, 4.1), rate = c(0.23, 0.24))
#'
#' dnsExp(4, rate = c(0.23, 0.23))
#' @export
dnsExp <- function(x, rate) {
  if (any(rate <= 0) || any(x < 0)) {
    stop("Arguments are incompatible.")
  }

  rate <- rep(rate, length.out = 2)
  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- rate[i] * exp(-x[, i] * rate[i])
  }

  return(pdf)
}
#' @name Neutrosophic Exponential
#' @examples
#'
#' # The cumulative distribution function for the nuetrosophic observation (4,4.1)
#' pnsExp(c(4, 4.1), rate = c(0.23, 0.24), lower.tail = TRUE)
#'
#' pnsExp(4, rate = c(0.23, 0.24))
#' @export
pnsExp <- function(q, rate, lower.tail = TRUE) {
  if (any(rate <= 0) || any(q < 0)) {
    stop("Arguments are incompatible.")
  }

  rate <- rep(rate, length.out = 2)
  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- 1 - exp(-q[, i] * rate[i])
  }


  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  return(cdf)
}
#' @name Neutrosophic Exponential
#' @examples
#' # The first percentile
#' qnsExp(p = 0.1, rate = 0.25)
#'
#' # The quantiles
#' qnsExp(p = c(0.25, 0.5, 0.75), rate = c(0.24, 0.25))
#'
#' @export
qnsExp <- function(p, rate) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(rate <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  rate <- rep(rate, length.out = 2)
  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- stats::qexp(p[, i], rate = rate[i])
  }

  return(quantiles)
}

#' @name Neutrosophic Exponential
#' @examples
#' # Simulate 10 numbers
#' rnsExp(n = 10, rate = c(0.23, 0.24))
#' @export
#'
rnsExp <- function(n, rate) {
  if (any(rate <= 0)) {
    stop(message = "Arguments are incompatible.")
  }
  rate <- rep(rate, length.out = 2)

  X <- qnsExp(runif(n), rate)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
