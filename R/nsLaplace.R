#' Neutrosophic Laplace (Double Exponential) Distribution
#'
#' Density, distribution function, quantile function, and random
#' generation for the neutrosophic Laplace (Double Exponential)
#' distribution with parameters \code{location} =  \eqn{\theta_N} and
#' \code{scale} = \eqn{\beta_N}.
#'
#' The neutrosophic Laplace distribution with parameters \eqn{\theta_N}
#' and \eqn{\beta_N} has density
#' \deqn{f_N(x) = \frac{1}{2\beta_N} \exp\left\{-\frac{|x-\theta_N|}{\beta_N}\right\}}
#' for \eqn{-\infty < x < \infty}, \eqn{\theta_N \in (\theta_L, \theta_U)}, the location parameter,
#' \eqn{\beta_N \in (\beta_L, \beta_U)}, the scale parameter which be a positive interval.
#'
#'
#' @name Neutrosophic Laplace
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param location the location parameter, which is the mean.
#' @param scale the scale parameter, Must consist of positive values.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#' \code{dnsLaplace} gives the density function
#'
#' \code{pnsLaplace} gives the distribution function
#'
#' \code{qnsLaplace} gives the quantile function
#'
#' \code{rnsLaplace} generates random values from the neutrosophic Laplace distribution.
#'
#' @references
#' Rahul, T., Malik, S. C., Raj, M. (2023). Neutrosophic Laplace Distribution
#' with Application in Financial Data Analysis, \emph{Neutrosophic Sets and Systems},
#' 57(1), 224-233.
#'
#' @importFrom stats runif
#'
#' @examples
#' dnsLaplace(x = c(4, 4.1), location = c(0.23, 0.24), scale = c(1, 2))
#' dnsLaplace(4, location = c(0.23, 0.24), scale = c(1, 2))
#'
#' @export
dnsLaplace <- function(x, location, scale) {
  if (any(scale < 0) || any(x < 0)) {
    stop("Arguments are incompatible.")
  }

  location <- rep(location, length.out = 2)
  scale <- rep(scale, length.out = 2)
  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- ((1 / (2 * scale[i]) * exp(-abs(x[, i] - location[i]) / scale[i])))
  }
  return(pdf)
}

#' @name Neutrosophic Laplace
#' @examples
#'
#' # The cumulative distribution function for the neutrosophic observation (4,4.1)
#' pnsLaplace(q = c(4, 4.1), location = c(0.23, 0.24), scale = c(1, 2))
#' @export
pnsLaplace <- function(q, location, scale, lower.tail = TRUE) {
  if (any(scale < 0) || any(q < 0)) {
    stop("Arguments are incompatible.")
  }

  location <- rep(location, length.out = 2)
  scale <- rep(scale, length.out = 2)
  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- ifelse(q[, i] < location[i], 0.5 * exp((q[, i] - location[i]) / scale[i]), 1 - 0.5 * exp((q[, i] - location[i])))
  }

  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  return(cdf)
}

#' @name Neutrosophic Laplace
#' @examples
#' # The first percentile
#' qnsLaplace(p = 0.1, location = 0.24, scale = 2)
#'
#' # The quantiles
#' qnsLaplace(p = c(0.25, 0.5, 0.75), location = c(0.23, 0.24), scale = c(1, 2))
#'
#' @export
qnsLaplace <- function(p, location, scale) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }
  if (any(scale < 0)) {
    stop(message = "Arguments are incompatible.")
  }
  location <- rep(location, length.out = 2)
  scale <- rep(scale, length.out = 2)
  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- ifelse(p[, i] <= 0.5, location[i] - scale[i] * log(2 * p[, i]), location[i] + scale[i] * log(2 * p[, i]))
  }

  return(quantiles)
}

#' @name Neutrosophic Laplace
#' @examples
#' # Simulate 10 numbers
#' rnsLaplace(n = 10, location = c(0.23, 0.24), scale = c(1, 2))
#' @export
#'
rnsLaplace <- function(n, location, scale) {
  if (any(scale < 0)) {
    stop(message = "Arguments are incompatible.")
  }
  location <- rep(location, length.out = 2)
  scale <- rep(scale, length.out = 2)

  X <- qnsLaplace(runif(n), location, scale)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
