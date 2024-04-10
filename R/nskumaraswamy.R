#' Neutrosophic Kumaraswamy Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic Kumaraswamy distribution with the
#' parameter \code{lambda} =  \eqn{\shape1_N}.
#'
#' The neutrosophic Kumaraswamy distribution with parameter \eqn{\shape1_N}
#' has density
#' \deqn{f_N(x)=\shape1_N \exp \left(-x \shape1_N\right)}
#' for \eqn{x \ge 0} and \eqn{\shape1_N \in (\shape1_L, \shape1_U)},
#' the lambda parameter must be a positive interval and \eqn{x \ge 0}.
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
#'  \code{pnskumaraswamy} gives the distribution function,
#'  \code{dnskumaraswamy} gives the density,
#'  \code{qnskumaraswamy} gives the quantile function and
#'  \code{rnskumaraswamy} genelambdas random values from the neutrosophic Kumaraswamy distribution.
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
#' data <- matrix(c(4, 4, 3.5, 3.5, 3.9, 4.1, 4.2, 4.2, 4.3, 4.6, 4.7, 4.7),
#'   nrow = 6, ncol = 2, byrow = TRUE
#' )
#'
#' dnskumaraswamy(x = c(4, 4.1), shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' dnskumaraswamy(4, shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' @export
dnskumaraswamy <- function(x, shape1, shape2) {
  if (any(shape1 <= 0) || any(shape2 <= 0) || any(x < 0)) {
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
#' pnskumaraswamy(q = c(4, 4.1), shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' @export
pnskumaraswamy <- function(q, shape1, shape2, lower.tail = TRUE) {
  if (any(shape1 <= 0) || any(shape2 <= 0) || any(q < 0)) {
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
#' qnskumaraswamy(p = 0.1, shape1 = 0.24, shape2 = 2)
#'
#' # The quantiles
#' qnskumaraswamy(p = c(0.25, 0.5, 0.75), shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#'
#' @export
qnskumaraswamy <- function(p, shape1, shape2) {
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
#' rnskumaraswamy(n = 10, shape1 = c(0.23, 0.24), shape2 = c(1, 2))
#' @export
#'
rnskumaraswamy <- function(n, shape1, shape2) {
  if (any(shape1 <= 0) || any(shape2 <= 0)) {
    stop(message = "Arguments are incompatible.")
  }
  shape1 <- rep(shape1, length.out = 2)
  shape2 <- rep(shape2, length.out = 2)

  X <- qnskumaraswamy(runif(n), shape1, shape2)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
