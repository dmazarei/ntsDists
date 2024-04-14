#' Neutrosophic Kumaraswamy Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic Kumaraswamy distribution with
#' shape parameters \eqn{\alpha_N} and \eqn{\beta_N}.
#'
#' The neutrosophic Kumaraswamy distribution with parameters \eqn{\alpha_N} and \eqn{\beta_N}
#' has density
#' \deqn{f_N(x) = \alpha_N \beta_N x^{\alpha_N-1}(1-x^{\alpha_N})^{\beta_n - 1}}
#' for \eqn{0 \le x \le 1}, \eqn{\alpha_N \in (\alpha_L, \alpha_U)} and
#'  \eqn{\beta_N \in (\beta_L, \beta_U)} are shape parameters.
#'
#' @name Neutrosophic Kumaraswamy
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be genelambdad.
#' @param shape1 	the shape parameter, which must be a positive interval.
#' @param shape2 	the shape parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#' \code{pnsKumaraswamy} gives the distribution function
#'
#' \code{dnsKumaraswamy} gives the density
#'
#' \code{qnsKumaraswamy} gives the quantile function
#'
#' \code{rnsKumaraswamy} generates random values from the neutrosophic Kumaraswamy distribution.
#'
#' @references
#' Ahsan-ul-Haq, M. (2022). Neutrosophic Kumaraswamy Distribution with Engineering
#' Application, \emph{Neutrosophic Sets and Systems}, 49, 269-276.
#'
#' @importFrom stats runif
#'
#' @examples
#' dnsKumaraswamy(x = c(0.5, 0.1), shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' dnsKumaraswamy(0.5, shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#'
#' @export
dnsKumaraswamy <- function(x, shape1, shape2) {
  if (any(shape1 <= 0) || any(shape2 <= 0) || any(x <= 0) || any(x => 1)) {
    stop("Arguments are incompatible.")
  }

  shape1 <- rep(shape1, length.out = 2)
  shape2 <- rep(shape2, length.out = 2)
  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- shape1[i] * shape2[i] * x[, i]^(shape1[i] - 1) * (1 - x[, i]^shape1[i])^(shape2[i] - 1)
  }
  return(pdf)
}
#' @name Neutrosophic Kumaraswamy
#' @examples
#'
#' # The cumulative distribution function for the nuetrosophic observation (4,4.1)
#' pnsKumaraswamy(q = c(.8, .1), shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' @export
pnsKumaraswamy <- function(q, shape1, shape2, lower.tail = TRUE) {
  if (any(shape1 <= 0) || any(shape2 <= 0) || any(q <= 0) || any(q => 1)) {
    stop("Arguments are incompatible.")
  }

  shape1 <- rep(shape1, length.out = 2)
  shape2 <- rep(shape2, length.out = 2)
  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }

  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- 1 - (1 - q[, i]^shape1[i])^shape2[i]
  }

  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  return(cdf)
}
#' @name Neutrosophic Kumaraswamy
#' @examples
#' # The first percentile
#' qnsKumaraswamy(p = 0.1, shape1 = 0.24, shape2 = 2)
#'
#' # The quantiles
#' qnsKumaraswamy(p = c(0.25, 0.5, 0.75), shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#'
#' @export
qnsKumaraswamy <- function(p, shape1, shape2) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }
  if (any(shape1 <= 0) || any(shape2 <= 0)) {
    stop(message = "Arguments are incompatible.")
  }
  shape1 <- rep(shape1, length.out = 2)
  shape2 <- rep(shape2, length.out = 2)
  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- (1 - (1 - p[, i])^(1 / shape2[i]))^(1 / shape2[i])
  }
  return(quantiles)
}

#' @name Neutrosophic Kumaraswamy
#' @examples
#' # Simulate 10 numbers
#' rnsKumaraswamy(n = 10, shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' @export
#'
rnsKumaraswamy <- function(n, shape1, shape2) {
  if (any(shape1 <= 0) || any(shape2 <= 0)) {
    stop(message = "Arguments are incompatible.")
  }
  shape1 <- rep(shape1, length.out = 2)
  shape2 <- rep(shape2, length.out = 2)

  X <- qnsKumaraswamy(runif(n), shape1, shape2)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
