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
#' @name NBIND
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param size number of trials (zero or more), which must be a positive interval.
#' @param prob probability of success on each trial, \eqn{0 < prob <= 1}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnbind} gives the distribution function,
#'  \code{dnbind} gives the density,
#'  \code{qnbind} gives the quantile function and
#'  \code{rnbind} generates random variables from the Binomial Distribution.
#' @references
#'        Granados, C. (2022). Some discrete neutrosophic distributions with
#'         neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5),
#'         1442-1457.
#' @importFrom stats runif dbinom pbinom qbinom
#' @examples
#' dnbind(x, size = 2, prob = 0.5)
#' @export
dnbind <- function(x, size, prob, log = FALSE) {
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

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dbinom(x[, i], size = size[i], prob = prob[i], log = log)
  }

  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name NBIND
#' @examples
#' x <- 1:10
#' x2 <- matrix(1:20, ncol = 2)
#' pnbind(x, size = 2, prob = 0.5)
#' pnbind(x2, size = c(2, 2), prob = c(.3, .6))
#' @export

pnbind <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE) {
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

  cdf <- stats::pbinom(q, size = size, prob = prob, lower.tail = lower.tail, log.p =log.p)


  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name NBIND
#' @examples
#' q1 <- seq(0.1, 1, length.out = 40)
#' qnbind(q1, size = 2, prob = 0.5)
#' q2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnbind(q2, lambda = c(2, 2))
#' @export
qnbind <- function(p, size, prob, lower.tail = TRUE, log.p = FALSE) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(size < 1) || any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  size <- rep(size, length.out = 2)
  prob <- rep(prob, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- stats::qbinom(p[, i], size = size[i], prob = prob[i], lower.tail = lower.tail, log.p =log.p)
  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name NBIND
#' @examples
#' n <- 10
#' rnbind(n, size = 2, prob = 0.5)
#' rnbind(n, lambda = c(1, 2))
#' @export
rnbind <- function(n, size, prob) {
  if (any(size < 1) || any(prob <= 0) || any(prob > 1)) {
    stop(message = "Arguments are incompatible.")
  }

  size <- rep(size, length.out = 2)
  prob <- rep(prob, length.out = 2)
  u <- matrix(runif(n), ncol = 2)
  X <- qnbind(u, size, prob)

  return(X)
}
