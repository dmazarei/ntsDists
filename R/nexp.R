#' Neutrosophic Exponential Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic exponential distribution with the
#' parameter \code{rate} =  \eqn{\theta_N}.
#'
#' The neutrosophic exponential distribution with parameter \eqn{\theta_N}
#' has density
#' \deqn{f_N(x)=\theta_N \exp \left(-x \theta_N\right)}
#' for \eqn{x \ge 0} and \eqn{\theta_N \in (\theta_L, \theta_U)},
#' the rate parameter must be a positive interval and \eqn{x \ge 0}.
#' @name NED
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param rate the shape parameter, which must be a positive interval.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnexp} gives the distribution function,
#'  \code{dnexp} gives the density,
#'  \code{qnexp} gives the quantile function and
#'  \code{rnexp} generates random values from the neutrosophic exponential distribution.
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
#' data <- matrix(c(4,4,3.5,3.5,3.9,4.1,4.2,4.2,4.3,4.6, 4.7, 4.7), nrow = 6, ncol = 2, byrow = TRUE)
#' # Note that we have a matrix of nuetrosophic observations
#'
#' # The density function of data with the estimated value for parameters based on Duan et al. (2021)
#' dnexp(data, rate = c(0.23, 0.24))
#' # The density function for the nuetrosophic observation (4,4.1)
#' dnexp(x = c(4,4.1), rate = c(0.23, 0.24))
#'
#' # The density function for the nuetrosophic observation 4
#' # Here, 4 is equivalent to c(4,4).
#' dnexp(4, rate = c(0.23, 0.23))
#' @export
dnexp <- function(x, rate, log = FALSE) {
  if (any(rate <= 0) || any(x < 0))
    stop("Arguments are incompatible.")

  rate <- rep(rate, length.out = 2)
  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- rate[i] * exp(-x[, i] * rate[i])
  }
  if(log == TRUE){
    pdf <- log(pdf)
  }
  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name NED
#' @examples
#'
#' # The cumulative distribution function for the nuetrosophic observation (4,4.1)
#' pnexp(c(4,4.1), rate = c(0.23, 0.24), lower.tail = TRUE)
#'
#' pnexp(4, rate = c(0.23, 0.24))
#' @export
pnexp <- function(q, rate, lower.tail = TRUE, log.p = FALSE) {
  if (any(rate <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  rate <- rep(rate, length.out = 2)
  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- 1 - exp(-q * rate)

  if (!lower.tail)
    cdf <- 1 - cdf
  if(log.p == TRUE){
    cdf <- log(cdf)
  }
  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name NED
#' @examples
#' # The first percentile
#' qnexp(p = 0.1, rate = 0.25)
#'
#' # The quantiles
#' qnexp(p = c(0.25, 0.5, 0.75), rate = c(0.24, 0.25))
#'
#' @export
qnexp <- function(p, rate, log.p = FALSE) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(rate <= 0)){
    stop(message = "Arguments are incompatible.")
  }

  rate <- rep(rate, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- -log(1 - p[, i]) / rate[i]
  }
  if(log.p == TRUE){
    quantiles <- log(quantiles)
  }
  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name NED
#' @examples
#' Simulate 10 numbers
#' rnexp(n = 10, rate = c(0.23, 0.24))
#' @export
#'
rnexp <- function(n, rate) {
  if (any(rate <= 0))
    stop(message = "Arguments are incompatible.")
  rate <- rep(rate, length.out = 2)
  u <- matrix(runif(n * length(rate)), nrow = n, ncol = length(rate))
  X <- qnexp(u, rate)

  return(X)
}
