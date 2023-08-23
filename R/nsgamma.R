#' Neutrosophic Gamma Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic gamma distribution with parameter \code{shape} = \eqn{\alpha_N}
#' and \code{scale}=\eqn{\lambda_N}.
#'
#' The neutrosophic gamma distribution with parameters \eqn{\alpha_N} and
#' \eqn{\lambda_N} has density
#' \deqn{f_n(x)=\frac{1}{\Gamma(\alpha_n) \lambda_n^{\alpha_n}} x^{\alpha_n-1} \exp\{-\left(x / \lambda_n\right)\}}
#' for \eqn{x \ge 0}, \eqn{\alpha_N \in (\alpha_L, \alpha_U)}, the shape
#' parameter which must be a positive interval and
#' \eqn{\lambda_n \in (\lambda_L, \lambda_U)}, the scale parameter which
#' must be a positive interval. Here, \eqn{\Gamma(\cdot)} is gamma
#' function implemented by \code{\link{gamma}}.
#'
#' @name Neutrosophic Gamma
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param shape the shape parameter, which must be a positive interval.
#' @param scale the scale parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsgamma} gives the distribution function,
#'  \code{dnsgamma} gives the density,
#'  \code{qnsgamma} gives the quantile function and
#'  \code{rnsgamma} generates random variables from the neutrosophic gamma distribution.
#' @references
#'    Khan, Z., Al-Bossly, A., Almazah, M. M. A., and Alduais, F. S. (2021).
#'    On statistical development of neutrosophic gamma distribution with
#'    applications to complex data analysis, \emph{Complexity}, 2021, Article ID 3701236.
#' @importFrom stats runif dgamma pgamma qgamma
#' @examples
#' data(remission)
#' dnsgamma(x = remission, shape = c(1.1884,1.1896), scale = c(7.6658,7.7796))
#'
#' pnsgamma(q = 20, shape = c(1.1884,1.1896), scale = c(7.6658,7.7796))
#'
#' # Calculate quantiles
#' qnsgamma(p = c(0.25,0.5,0.75), shape = c(1.1884,1.1896), scale = c(7.6658,7.7796))
#'
#' # Simulate 10 numbers
#' rnsgamma(n=10, shape = c(1.1884,1.1896), scale = c(7.6658,7.7796))
#' @export
dnsgamma <- function(x, shape, scale) {
  if (any(shape <= 0) || any(scale <= 0) || any(x < 0))
    stop(message = "Arguments are incompatible.")

  shape <- rep(shape, length.out = 2)
  scale  <- rep(scale, length.out = 2)

  if (is.vector(x) && length(x) == 1) {
    x <- matrix(rep(x, each = 2), ncol = 2, byrow = TRUE)
  }

  x <- matrix(x, ncol = 2)

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- stats::dgamma(x[, i], shape = shape[i], scale = scale[i])
  }

  # Identify rows where col1 > col2
  swap_rows <- pdf[, 1] > pdf[, 2]
  # Swap values using logical indexing
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name Neutrosophic Gamma
#' @export
pnsgamma <- function(q, shape, scale, lower.tail = TRUE) {
  if (any(shape <= 0) || any(scale <= 0) || any(q < 0))
    stop(message = "incompatible arguments.")

  shape <- rep(shape, length.out = 2)
  scale  <- rep(scale, length.out = 2)

  if (is.vector(q)){
    q <- matrix(rep(q, each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q)>2){
    stop(message = "Arguments are incompatible.")
  }

  cdf <- stats::pgamma(q, shape = shape, scale = scale)

  if (!lower.tail)
    cdf <- 1 - cdf

  cdf <- matrix(cdf, ncol = 2, byrow = TRUE)
  # Identify rows where col1 > col2
  swap_rows <- cdf[, 1] > cdf[, 2]
  # Swap values using logical indexing
  cdf[swap_rows, c(1, 2)] <- cdf[swap_rows, c(2, 1)]

  return(cdf)
}

#' @name Neutrosophic Gamma
#' @export
qnsgamma <- function(p, shape, scale) {
  if (any(p < 0) || any(p > 1))
    stop(message = "Warning: p should be in the interval [0,1].")

  if (any(shape <= 0) || any(scale <= 0))
    stop(message = "Arguments are incompatible.")

  shape <- rep(shape, length.out = 2)
  scale  <- rep(scale, length.out = 2)

  if (is.vector(p)){
    p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(p)>2){
    stop(message = "Arguments are incompatible.")
  }

  quantiles <- stats::qgamma(p, shape = shape, scale = scale)
  quantiles <- matrix(quantiles, ncol = 2, byrow = TRUE)

  swap_rows <- quantiles[, 1] > quantiles[, 2]
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

  return(quantiles)
}

#' @name Neutrosophic Gamma
#' @export
rnsgamma <- function(n, shape, scale) {
  if (any(shape <= 0) || any(scale <= 0))
    stop(message = "Arguments are incompatible.")

  shape <- rep(shape, length.out = 2)
  scale  <- rep(scale, length.out = 2)

  X <- qnsgamma(runif(n), shape, scale)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
