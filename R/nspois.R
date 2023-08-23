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
#' @name Neutrosophic Poisson
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param lambda the mean, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnspois} gives the distribution function,
#'  \code{dnspois} gives the density,
#'  \code{qnspois} gives the quantile function and
#'  \code{rnspois} generates random variables from the neutrosophic Poisson Distribution.
#' @references
#'        Alhabib, R., Ranna, M. M., Farah, H., Salama, A. A. (2018).
#'        Some neutrosophic probability distributions.
#'        \emph{Neutrosophic Sets and Systems},  22, 30-38.
#' @importFrom stats runif dpois ppois qpois
#' @examples
#' # In a company, Phone employee receives phone calls, the calls arrive with
#' # rate of [1 , 3] calls per minute, we will calculate
#' # the probability that the employee will not receive any call within a minute
#' dnspois(x = 0, lambda = c(1,3))
#'
#' # the probability that employee would not receive any call within 5 minutes
#' dnspois(x = 0, lambda = c(5,15))
#' @export
dnspois <- function(x, lambda) {
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
#' @name Neutrosophic Poisson
#' @examples
#' # the probability that the employee will receive at least one call within a minute
#' pnspois(q = 1, lambda = c(1,3), lower.tail = FALSE)
#' # the probability that the employee will receive at most three calls within 5 minutes
#' pnspois(q = 3, lambda = c(5,15), lower.tail = TRUE)
#' @export
pnspois <- function(q, lambda, lower.tail = TRUE) {
  if (any(lambda < 0))
    stop("Arguments are incompatible.")

  if (any(q < 0) && any(q - floor(q) == 0))
    stop(message = "Warning: q should be a  positive integer.")

  lambda <- rep(lambda, length.out = 2)

  if (is.vector(q)){
    q <- matrix(rep(q, each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q)>2){
    stop(message = "Arguments are incompatible.")
  }

  cdf <- stats::ppois(q, lambda)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name Neutrosophic Poisson
#' @examples
#' # Calcaute the quantiles
#' qnspois(p = c(0.25,0.5,0.75), lambda = c(1,3))
#' @export
qnspois <- function(p, lambda) {
  if (any(p < 0) || any(p > 1))
    stop(message = "Warning: p should be in the interval [0,1].")

  if (any(lambda < 0))
    stop(message = "Arguments are incompatible.")

  lambda <- rep(lambda, length.out = 2)

  if (is.vector(p)){
    p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p)>2){
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- stats::qpois(p, lambda = lambda)
  quantiles <- matrix(quantiles, ncol = 2, byrow = TRUE)

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Poisson
#' @examples
#' # Simulate 10 values
#' rnspois(n = 10, lambda = 1)
#' @export
rnspois <- function(n, lambda) {
  if (any(lambda < 0))
    stop(message = "Arguments are incompatible.")

  lambda <- rep(lambda, length.out = 2)

  X <- qnspois(runif(n), lambda)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
