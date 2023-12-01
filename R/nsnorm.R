#' Neutrosophic Normal Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic generalized exponential
#' distribution with parameters \code{mean} =  \eqn{\mu_N} and standard deviation
#' \code{sd} = \eqn{\sigma_N}.
#'
#' The neutrosophic normal distribution with parameters mean
#' \eqn{\mu_N} and standard deviation \eqn{\sigma_N} has density function
#' \deqn{f_N(x) = \frac{1}{\sigma_N \sqrt{2 \pi}} \exp\{\left(\frac{\left(X-\mu_N\right)^2}{2 \sigma_N^2}\right)}\}
#' for \eqn{\mu_N \in (\mu_L, \mu_U)}, the mean which must be an interval, and
#' \eqn{\sigma_N \in (\sigma_L, \sigma_U)}, the standard deviation which must
#' also be a positive interval, and \eqn{-\infty < x < \infty}.
#'
#'
#' @name Neutrosophic Normal
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param mean the mean, which must be an interval.
#' @param sd the standard deviations that must be positive.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#' \code{pnsnorm} gives the distribution function,
#' \code{dnsnorm} gives the density,
#' \code{qnsnorm} gives the quantile function and
#' \code{rnsnorm} generates random variables from the neutrosophic normal distribution.
#'
#' @references
#'    Patro, S. and Smarandache, F. (2016). The Neutrosophic Statistical Distribution, More Problems, More Solutions. Infinite Study.
#'
#' @importFrom stats runif dnorm pnorm qnorm
#' @examples
#' data(balls)
#' dnsnorm(x = balls, mean = c(72.14087, 72.94087), sd = c(37.44544, 37.29067))
#'
#' pnsnorm(q = 5, mean = c(72.14087, 72.94087), sd = c(37.44544, 37.29067))
#'
#' # Calculate quantiles
#' qnsnorm(p = c(0.25, 0.5, 0.75), mean = c(9.1196, 9.2453), sd = c(10.1397, 10.4577))
#'
#' # Simulate 10 values
#' rnsnorm(n = 10, mean = c(4.141, 4.180), sd = c(0.513, 0.521))
#'
#' @export
dnsnorm <- function(x, mean, sd) {
  if (any(sd <= 0)) {
    stop("Arguments are incompatible.")
  }

  mean <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }


  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- stats::dnorm(x[, i], mean = mean[i], sd = sd[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}

#' @name Neutrosophic Normal
#' @export

pnsnorm <- function(q, mean, sd, lower.tail = TRUE) {
  if (any(sd <= 0)) {
    stop("Arguments are incompatible.")
  }

  mean <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }



  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- stats::pnorm(q[, i], mean = mean[i], sd = sd[i])
  }


  if (!lower.tail) {
    cdf <- 1 - cdf
  }


  return(cdf)
}

#' @name Neutrosophic Normal
#' @export
qnsnorm <- function(p, mean, sd) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(sd <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  mean <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- stats::qnorm(p[, i], mean = mean[i], sd = sd[i])
  }

  return(quantiles)
}

#' @name Neutrosophic Normal
#' @export
rnsnorm <- function(n, mean, sd) {
  if (any(sd <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  mean <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  X <- qnsnorm(runif(n), mean, sd)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
