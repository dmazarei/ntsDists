#' Neutrosophic Hypergeometric Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic hypergeometric distribution with
#' parameters \eqn{m_N}, \eqn{n_N}, and \eqn{k_N}.
#'
#' The neutrosophic hypergeometric distribution with parameters \eqn{k_N}
#' and \eqn{n_N} and \eqn{m_N}  has the density
#' \deqn{\frac{\left(\begin{array}{c}
#' m_N \\
#' x
#' \end{array}\right)\left(\begin{array}{c}
#'                        n_N \\
#'                        k_N-x
#'                        \end{array}\right)}{\left(\begin{array}{c}
#'                                                  m_N + n_N \\
#'                                                  k_N
#'                                                  \end{array}\right)}}
#' for \eqn{n_N \in (n_L, n_U)} which must be a positive interval and
#'  \eqn{k_N \in (k_L, k_U)} which must be a positive interval and
#'  \eqn{k \in \{0, 1, 2, \ldots\}} and
#' and \eqn{x \in \{0, 1, 2, \ldots, k\}}.
#'
#' @name Neutrosophic Hypergeometric
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param nn number of random values to be generated.
#' @param N number of population size, which must be a positive interval.
#' @param K number of success states in the population, which must be a positive interval.
#' @param k number of draws (i.e. quantity drawn in each trial), , which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnshyper} gives the distribution function,
#'  \code{dnshyper} gives the density,
#'  \code{qnshyper} gives the quantile function and
#'  \code{rnshyper} generates random variables from the HyperGeometric Distribution.
#' @references
#'        Granados, C. (2022).
#'        Some discrete neutrosophic distributions with neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5), 1442-1457.
#' @importFrom stats runif dhyper phyper qhyper
#' @examples
#' dnshyper(x, N, K, k)
#' dnshyper(x2, N, K, k)
#' @export
dnshyper <- function(x, N, K, k) {
  if (any(K <= 0) || any(N <= 0) || any(k <= 0) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
  }

  if (any(x < 0) && any(x - floor(x) == 0)) {
    stop(message = "Warning: x should be a positive integer.")
  }

  K <- rep(K, length.out = 2)
  N <- rep(N, length.out = 2)
  k <- rep(k, length.out = 2)

  if (is.vector(x)) {
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dhyper(x[, i], m = K[i], n = N[i] - K[i], k = k[i])
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Hypergeometric
#' @examples
#' x <- 1:10
#' x2 <- matrix(1:20, ncol = 2)
#' pnshyper(x, N, K, k)
#' pnshyper(x2, N, K, k)
#' @export

pnshyper <- function(q, N, K, k, lower.tail = TRUE) {
  if (any(K <= 0) || any(N <= 0) || any(k <= 0) || any(q < 0)) {
    stop(message = "Arguments are incompatible.")
  }
  if (any(q < 0) && any(q - floor(q) == 0)) {
    stop(message = "Warning: q should be a  positive integer.")
  }

  K <- rep(K, length.out = 2)
  N <- rep(N, length.out = 2)
  k <- rep(k, length.out = 2)

  if (is.vector(q)) {
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::phyper(q, m = K, n = N - K, k = k)

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name Neutrosophic Hypergeometric
#' @examples
#' q1 <- seq(0.1, 1, length.out = 40)
#' qnshyper(q1, N, K, k)
#' q2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnshyper(q2, N, K, k)
#' @export
qnshyper <- function(p, N, K, k) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(K <= 0) || any(N <= 0) || any(k <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  K <- rep(K, length.out = 2)
  N <- rep(N, length.out = 2)
  k <- rep(k, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qhyper(p[, i], m = K[i], n = N[i] - K[i], k = k[i])
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Hypergeometric
#' @examples
#' n <- 10
#' rnshyper(n, N, K, k)
#' rnshyper(n, N, K, k)
#' @export
rnshyper <- function(nn, m, n, k) {
  if (any(k <= 0) || any(N <= 0) || any(k <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  K <- rep(K, length.out = 2)
  N <- rep(N, length.out = 2)
  k <- rep(k, length.out = 2)
  u <- matrix(runif(n), ncol = 2)
  X <- qnshyper(u, N, K, k)

  return(X)
}
