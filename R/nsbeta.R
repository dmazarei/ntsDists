#' Neutrosophic Beta Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the neutrosophic Beta distribution with shape parameters
#' \code{shape1} = \eqn{\alpha_N} and \code{shape2} = \eqn{\beta_N}.
#'
#' The neutrosophic beta distribution with parameters \eqn{\alpha_N} and
#' \eqn{\beta_N} has the probability density function
#' \deqn{f_N(X) = \frac{1}{B(\alpha_N, \beta_N)} X^{\alpha_N - 1} (1 - X)^{\beta_N - 1}}
#' for \eqn{\alpha_N \in (\alpha_L, \alpha_U)}, the first shape parameter which
#' must be a positive interval, and \eqn{\beta_N \in (\beta_L, \beta_U)},
#' the second shape parameter which must also be a positive interval, and
#' \eqn{0 \le x \le 1}. The function \eqn{B(a, b)}
#' returns the beta function and can be calculated using \code{\link{beta}}.
#'
#' @name Neutrosophic Beta
#'
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param shape1 the first shape parameter, which must be a positive interval.
#' @param shape2 the second shape parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#' \code{pnsbeta} gives the distribution function,
#' \code{dnsbeta} gives the density,
#' \code{qnsbeta} gives the quantile function and
#' \code{rnsbeta} generates random values from the neutrosophic Beta distribution.
#'
#' @references
#'  Sherwani, R. Ah. K., Naeem, M., Aslam, M., Reza, M. A., Abid, M., Abbas, S. (2021).
#'     Neutrosophic beta distribution with properties and applications.
#'     \emph{Neutrosophic Sets and Systems}, 41, 209-214.
#' @importFrom stats runif dbeta pbeta qbeta
#' @examples
#'
#' dnsbeta(x = c(0.1, 0.2), shape1 = c(1, 1), shape2 = c(2, 2))
#' dnsbeta(x = 0.1, shape1 = c(1, 1), shape2 = c(2, 2))
#'
#' x <- matrix(c(0.1, 0.1, 0.2, 0.3, 0.5, 0.5), ncol = 2, byrow = TRUE)
#' dnsbeta(x, shape1 = c(1, 2), shape2 = c(2, 3))
#'
#' @export
dnsbeta <- function(x, shape1, shape2) {
  if (any(shape1 <= 0) || any(shape2 <= 0) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
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
    pdf[, i] <- stats::dbeta(x[, i], shape1 = shape1[i], shape2 = shape2[i])
  }

  return(pdf)
}
#' @name Neutrosophic Beta
#' @examples
#'
#' pnsbeta(q = c(0.1, 0.1), shape1 = c(3, 1), shape2 = c(1, 3), lower.tail = FALSE)
#' pnsbeta(x, shape1 = c(1, 2), shape2 = c(2, 2))
#'
#' @export
pnsbeta <- function(q, shape1, shape2, lower.tail = TRUE) {
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
    cdf[, i] <- stats::pbeta(q[, i], shape1 = shape1[i], shape2 = shape2[i])
  }

  if (!lower.tail) {
    cdf <- 1 - cdf
  }


  return(cdf)
}
#' @name Neutrosophic Beta
#' @examples
#'
#' qnsbeta(p = 0.1, shape1 = c(1, 1), shape2 = c(2, 2))
#' qnsbeta(p = c(0.25, 0.5, 0.75), shape1 = c(1, 2), shape2 = c(2, 2))
#'
#' @export
qnsbeta <- function(p, shape1, shape2) {
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
    quantiles[, i] <- stats::qbeta(p[, i], shape1 = shape1[i], shape2 = shape2[i])
  }

  return(quantiles)
}
#' @name Neutrosophic Beta
#' @examples
#' # Simulate 10 numbers
#' rnsbeta(n = 10, shape1 = c(1, 2), shape2 = c(1, 1))
#' @export
rnsbeta <- function(n, shape1, shape2) {
  if (any(shape1 <= 0) || any(shape2 <= 0)) {
    stop(message = "Arguments are incompatible.")
  }
  shape1 <- rep(shape1, length.out = 2)
  shape2 <- rep(shape2, length.out = 2)

  X <- qnsbeta(runif(n), shape1, shape2)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
