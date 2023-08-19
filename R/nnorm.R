#' Neutrosophic Normal Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic generalized exponential
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
#' @name NND
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param mean the mean, which must be an interval.
#' @param sd the standard deviations that must be positive.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#' \code{pnnorm} gives the distribution function,
#' \code{dnnorm} gives the density,
#' \code{qnnorm} gives the quantile function and
#' \code{rnnorm} generates random variables from the neutrosophic normal distribution.
#'
#' @references
#'    Patro, S. and Smarandache, F. (2016). The Neutrosophic Statistical Distribution, More Problems, More Solutions. Infinite Study.
#'
#' @importFrom stats runif dnorm pnorm qnorm
#' @examples
#' dnnorm(x = 0.5, mean = 1, sd = 2)
#'
#' x1 <- c(-0.8, 0.2, 1.6, 3.9)
#' dnnorm(x = 0.5, mean = c(1, 1), sd = c(2, 2))
#'
#' x2 <- matrix(seq(-3, 3, length.out = 10), nrow = 5, ncol = 2)
#' dnnorm(x2, mean = c(1, 2), sd = c(2, 2))
#'
#' @export
dnnorm <- function(x, mean, sd, log = FALSE) {
  if (any(sd <= 0))
    stop("Arguments are incompatible.")

  mean    <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  if (is.vector(x) && length(x) == 1) {
    x <- matrix(rep(x, each = 2), ncol = 2, byrow = TRUE)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dnorm(x[, i], mean = mean[i], sd = sd[i], log = log)
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)}

#' @name NND
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' pnnorm(x, mean = c(2, 2), sd = c(1, 1))
#'
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnnorm(x2, mean = c(1, 2), sd = c(2, 2))
#' @export

pnnorm <- function(q, mean, sd, lower.tail = TRUE, log.p = FALSE) {
  if (any(sd <= 0))
    stop("Arguments are incompatible.")

  mean    <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pnorm(q, mean = mean, sd = sd, lower.tail = lower.tail, log.p =log.p)

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)
  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}

#' @name NND
#' @examples
#' q1 <- c(0.01)
#' qnnorm(q1, mean = 1, sd = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnnorm(x2, mean = c(1, 2), sd = c(2, 2))
#' @export
qnnorm <- function(p, mean, sd, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(sd <= 0))
    stop(message = "Arguments are incompatible.")

  mean    <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qnorm(p[, i], mean = mean[i], sd = sd[i], lower.tail = lower.tail, log.p =log.p)
  }
  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name NND
#' @examples
#' n <- 100
#' rnnorm(n, mean = c(1, 1), sd = c(2, 2))
#' rnnorm(n, mean = c(1, 2), sd = c(1, 1))
#' @export
rnnorm <- function(n, mean, sd) {
  if (any(sd <= 0))
    stop(message = "Arguments are incompatible.")

  mean    <- rep(mean, length.out = 2)
  sd <- rep(sd, length.out = 2)

  u <- matrix(runif(n), ncol = 2)
  X <- qnnorm(u, mean, sd)

  return(X)
}
