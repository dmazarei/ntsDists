#' Neutrosophic Binomial Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic binomial distribution with
#' parameters \code{size} = \eqn{n} and \code{prob} = \eqn{p_N}.
#'
#' The neutrosophic binomial distribution with parameters \eqn{n} and \eqn{p_N}
#' has the density
#' \deqn{f_X(x)=\bigg(\begin{array}{c}n \\ x\end{array}\bigg) p_N^{x}\left(1-p_N\right)^{n-x}}
#' for \eqn{n \in \{1, 2, \ldots\}} and \eqn{p_N \in (p_L, p_U)} which must be \eqn{0<p_N<1}
#' and \eqn{x \in \{0, 1, 2, \ldots\}}.
#'
#' @name Neutrosophic Binomial
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param size number of trials (zero or more), which must be a positive interval.
#' @param prob probability of success on each trial, \eqn{0 \le prob \le 1}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsbinom} gives the distribution function,
#'  \code{dnsbinom} gives the density,
#'  \code{qnsbinom} gives the quantile function and
#'  \code{rnsbinom} generates random variables from the Binomial Distribution.
#' @references
#'        Granados, C. (2022). Some discrete neutrosophic distributions with
#'         neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5),
#'         1442-1457.
#' @importFrom stats runif dbinom pbinom qbinom
#' @examples
#' # Probability of X = 17 when X follows bin(n = 20, p = [0.8,0.9])
#' dnsbinom(x = 17, size = 20, prob = c(0.8,0.9))
#'
#' @export
dnsbinom <- function(x, size, prob) {
  if (any(size < 1) || any(prob <= 0) || any(prob > 1) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
  }

  if (any(x < 0) && any(x - floor(x) == 0)) {
    stop(message = "Warning: x should be a positive integer.")
  }

  size <- rep(size, length.out = 2)
  prob <- rep(prob, length.out = 2)

  if (is.vector(x)) {
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- stats::dbinom(x, size = size, prob = prob)
  pdf <- matrix(pdf, ncol = 2, byrow = TRUE)

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Binomial
#' @examples
#'
#' pnsbinom(q = 17, size = 20, prob = c(0.8,0.9))
#'
#' @export

pnsbinom <- function(q, size, prob, lower.tail = TRUE) {
  if (any(size < 1) || any(prob <= 0) || any(prob > 1) || any(q < 0)) {
    stop(message = "Arguments are incompatible.")
  }
  if (any(q < 0) && any(q - floor(q) == 0)) {
    stop(message = "Warning: q should be a  positive integer.")
  }

  size <- rep(size, length.out = 2)
  prob <- rep(prob, length.out = 2)
  if (is.vector(q)) {
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- stats::pbinom(q, size = size, prob = prob)

  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name Neutrosophic Binomial
#' @examples
#'
#' qnsbinom(p = 0.5, size = 20, prob = c(0.8,0.9))
#'
#' @export
qnsbinom <- function(p, size, prob) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(size < 1) || any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  size <- rep(size, length.out = 2)
  prob <- rep(prob, length.out = 2)

  if (is.vector(p)){
    p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(p)>2){
    stop(message = "Arguments are incompatible.")
  }

  quantiles <- stats::qbinom(p, size = size, prob = prob)
  quantiles <- matrix(quantiles, ncol = 2, byrow = TRUE)

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Binomial
#' @examples
#'
#' # Simulate 10 numbers
#' rnsbinom(n = 10, size = 20, prob = c(0.8,0.9))
#'
#' @export
rnsbinom <- function(n, size, prob) {
  if (any(size < 1) || any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  size <- rep(size, length.out = 2)
  prob <- rep(prob, length.out = 2)

  u <- matrix(runif(n), ncol = 2)
  X <- qnsbinom(u, size, prob)

  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
