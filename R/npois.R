#' Neutrosophic Poisson Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Poisson distribution with
#' parameter \eqn{\lambda_N}.
#'
#' The neutrosophic Poisson distribution with parameter \eqn{\lambda_N}
#' has the density
#' \deqn{f_N(x)= \exp\{-\lambda_N\}  \frac{\left(\lambda_N\right)^x}{x !}}
#' for \eqn{\lambda_N \in (\lambda_L, \lambda_U)} which must be a positive
#' interval and \eqn{x \in \{0, 1, 2, \ldots\}}.
#'
#' @name NPD
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param lambda the mean, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnpois} gives the distribution function,
#'  \code{dnpois} gives the density,
#'  \code{qnpois} gives the quantile function and
#'  \code{rnpois} generates random variables from the neutrosophic Poisson Distribution.
#' @references
#'        Alhabib, R., Ranna, M. M., Farah, H., Salama, A. A. (2018).
#'        Some neutrosophic probability distributions.
#'        \emph{Neutrosophic Sets and Systems},  22, 30-38.
#' @importFrom stats runif dpois ppois qpois
#' @examples
#' dnpois(x, lambda = 2)
#' dnpois(x2, lambda = c(2, 2))
#' @export
dnpois <- function(x, lambda) {
  if (any(lambda < 0))
    stop(message = "Arguments are incompatible.")

  if (any(x < 0) && any(x - floor(x) == 0))
    stop(message = "Warning: x should be a positive integer.")

  lambda <- rep(lambda, length.out = 2)
  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dpois(x[, i], lambda = lambda[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name NPD
#' @examples
#' x <- 1:10
#' x2 <- matrix(1:20, ncol = 2)
#' pnpois(x, lambda = 1)
#' pnpois(x2, lambda = c(2, 3))
#' @export
pnpois <- function(q, lambda, lower.tail = TRUE) {
  if (any(lambda < 0))
    stop("Arguments are incompatible.")

  if (any(q < 0) && any(q - floor(q) == 0))
    stop(message = "Warning: q should be a  positive integer.")

  lambda <- rep(lambda, length.out = 2)
  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::ppois(q, lambda)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name NPD
#' @examples
#' q1 <- seq(0.1, 1, length.out = 40)
#' qnpois(q1, lambda = 2)
#' q2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnpois(q2, lambda = c(2, 2))
#' @export
qnpois <- function(p, lambda) {
  if (any(p < 0) || any(p > 1))
    stop(message = "Warning: p should be in the interval [0,1].")

  if (any(lambda < 0))
    stop(message = "Arguments are incompatible.")

  lambda <- rep(lambda, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qpois(p[, i], lambda = lambda[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name NPD
#' @examples
#' n <- 10
#' rnpois(n, lambda = 1)
#' rnpois(n, lambda = c(1, 2))
#' @export
rnpois <- function(n, lambda) {
  if (any(lambda < 0))
    stop(message = "Arguments are incompatible.")

  lambda <- rep(lambda, length.out = 2)
  u <- matrix(runif(n), ncol = 2)
  X <- qnpois(u, lambda)

  return(X)
}
