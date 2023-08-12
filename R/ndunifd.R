#' Neutrosophic Discrete Uniform Distribution (NDUNIFD)
#'
#' Density, distribution function, quantile function and random
#' generation for the nuetrosophic Discrete Uniform distribution with
#' parameters \eqn{a_N} and  \eqn{b_N}.
#'
#' The neutrosophic Discrete Uniform distribution with parameters
#' \eqn{\max_N} and \eqn{\min_N} has the density
#' \deqn{f_N(x)=\frac{1}{b_N-a_N+1}}
#' for \eqn{a_N \in (a_L, a_U)}  lower parameter interval, \eqn{b_N \in (b_L,b_U)},
#'  upper parameter interval.
#'
#' @name NDUNIFD
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param min lower limits of the distribution. Must be finite.
#' @param max upper limits of the distribution. Must be finite.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pndunifd} gives the distribution function,
#'  \code{dndunifd} gives the density,
#'  \code{qndunifd} gives the quantile function and
#'  \code{rndunifd} generates random variables from the neutrosophic Discrete Uniform Distribution.
#' @references
#'        Granados, C. (2022).
#'        Some discrete neutrosophic distributions with neutrosophic parameters based on neutrosophic random variables.
#'         \emph{Hacettepe Journal of Mathematics and Statistics}, 51(5), 1442-1457.
#'
#' @examples
#' dndunifd(x, min = 1, max = 2)
#'
#' dndunifd(x2, min = c(1, 2), max = c(2, 2))
#' @export
dndunifd <- function(x, min, max) {
  if (any(max <= min))
    stop(message = "Arguments are incompatible.")

  if (any(round(x) == x))
    warning(paste("non-integer"))

  max <- rep(max, length.out = 2)
  min  <- rep(min, length.out = 2)

  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- ifelse(x[,i]>=min[i] & x[,i]<=max[i] & round(x[,i])==x[,i], 1/(max[i]-min[i]+1), 0)
  }


  swap_rows <- pdf[, 1] > pdf[, 2]
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pndunifd(x, min = 1, max = 2)
#' pndunifd(x2, min = c(1, 2), max = c(2, 2))
#' @export
pndunifd <- function(q, min, max) {
  if (any(max <= min))
    stop(message = "Arguments are incompatible.")

  q <- floor(q)
  max <- rep(max, length.out = 2)
  min  <- rep(min, length.out = 2)

  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  q <- matrix(q, ncol = 2)

  cdf <- ifelse(q<min,0,ifelse(q>max),1,(q-min+1)/(max-min+1))

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)

  swap_rows <- cdf[, 1] > cdf[, 2]
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}
#' @name NDUNIFD
#' @name NDUNIFD
#' @examples
#' qndunifd(x, min = 1, max = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qndunifd(x2, min = c(1, 2), max = c(2, 2))
#' @export
qndunifd <- function(p, min, max) {
  if (any(max <= min))
    stop(message = "Arguments are incompatible.")
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  max <- rep(max, length.out = 2)
  min  <- rep(min, length.out = 2)

  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- ceiling((max[i]-min[i]+1)*p[,i]+min[i]-1)

  }

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}
#' @name NDUNIFD
#' @examples
#' n <- 10
#' rndunifd(n, min = 1, max = 2)
#' rndunifd(n, min = c(1, 2), max = c(1, 1))
#' @export
rndunifd <- function(n, min, max) {
  if (any(max <= min))
    stop(message = "Arguments are incompatible.")
  max <- rep(max, length.out = 2)
  min  <- rep(min, length.out = 2)

  X <- matrix(NA, nrow = n, ncol = 2)
  for(i in 1:nrow(X))
  X[,i] <- sample(min[i]:max[i],size=n,replace = TRUE)

  return(X)
}
