#' Neutrosophic Beta Distribution (NBD)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic Beta distribution with parameters \code{shape1}=\eqn{\alpha_N}
#' and \code{shape2}=\eqn{\beta_N}.
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
#' @name NBD
#'
#' @param x a vector or matrix of quantiles for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param alpha the first shape parameter, which must be a positive interval.
#' @param beta the second shape parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#' \code{pnbd} gives the distribution function,
#' \code{dnbd} gives the density,
#' \code{qnbd} gives the quantile function and
#' \code{rnbd} generates random values from the neutrosophic Beta distribution.
#'
#' @references
#'  Sherwani, R. Ah. K., Naeem, M., Aslam, M., Reza, M. A., Abid, M., Abbas, S. (2021).
#'     Neutrosophic beta distribution with properties and applications.
#'     \emph{Neutrosophic Sets and Systems}, 41, 209-214.
#' @importFrom stats runif dbeta pbeta qbeta
#' @examples
#' dnbd(x = c(0.1, 0.2), alpha = c(1,1), beta = c(2, 2))
#'
#' dnbd(x = c(0.1, 0.1), alpha = c(0.5,0.7), beta = c(0.2, 2))
#'
#' x <- matrix(c(0.1, 0.1, 0.2, 0.3, 0.5, 0.5), ncol = 2, byrow = TRUE)
#' dnbd(x, alpha = c(1, 2), beta = c(2, 3))
#'
#' @export
dnbd <- function(x, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0) || any(x < 0))
    stop(message = "Arguments are incompatible.")

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)


  if (is.vector(x) && length(x) == 1) {
    x <- matrix(rep(x, each = 2), ncol = 2, byrow = TRUE)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dbeta(x[, i], shape1 = alpha[i], shape2 = beta[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name NBD
#' @examples
#' pnbd(q = c(0.1, 0.1), alpha = c(3, 1), beta = c(1,3), lower.tail = FALSE)
#'
#' pnbd(x, alpha = c(1, 2), beta = c(2, 2))
#'
#' @export
pnbd <- function(q, alpha, beta, lower.tail = TRUE) {
  if (any(theta <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)
  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pbeta(q, shape1 = alpha, shape2 = beta)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)
  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name NBD
#' @examples
#'
#' qnbd(p = 0.1, alpha = c(1,1), beta = c(2, 2))
#'
#' qnbd(p = c(0.25, 0.5, 0.75), alpha = c(1, 2), beta = c(2, 2))
#'
#' @export
qnbd <- function(p, alpha, beta) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(alpha <= 0) || any(beta <= 0))
    stop(message = "Arguments are incompatible.")

  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)


  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qbeta(p[, i], shape1 = alpha[i], shape2 = beta[i])
  }
  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}
#' @name NBD
#' @examples
#' Simulate 10 numbers
#' rnbd(n = 10, alpha = c(1, 2), beta = c(1, 1))
#' @export
rnbd <- function(n, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0))
    stop(message = "Arguments are incompatible.")
  alpha <- rep(alpha, length.out = 2)
  beta  <- rep(beta, length.out = 2)

  u <- matrix(runif(n), ncol = 2)
  X <- qnbd(u, alpha, beta)
  return(X)
}
