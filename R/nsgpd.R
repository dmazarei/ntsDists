#' Neutrosophic Generalized Pareto Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the neutrosophic generalized pareto distribution with parameters \code{shape} = \eqn{\alpha_N}
#' and \code{scale}=\eqn{\beta_N}.
#'
#' The neutrosophic generalized pareto distribution with parameters \eqn{\alpha_N} and
#' \eqn{\beta_N} has density
#' \deqn{f(x; \alpha_N, \beta_N)=\frac{1}{\beta_N}\left(1+\frac{\alpha_Nx}{\beta_N} \right)^{-\frac{1}{\alpha_N}-1},}
#' for \eqn{x \ge 0}, \eqn{\alpha_N \in (\alpha_L, \alpha_U)}, the shape
#' parameter which must be a positive interval and
#' \eqn{\beta_N \in (\beta_L, \beta_U)}, the scale parameter which
#' must be a positive interval.
#'
#' @name Neutrosophic Generalized Pareto
#' @param x a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param shape the shape parameter, which must be a positive interval.
#' @param scale the scale parameter, which must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \leq x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pnsgpd} gives the distribution function,
#'  \code{dnsgpd} gives the density,
#'  \code{qnsgpd} gives the quantile function and
#'  \code{rnsgpd} generates random variables from the neutrosophic generalized pareto distribution.
#' @references
#'    Eassa, N. I., Zaher, H. M., & El-Magd, N. A. A. (2023).
#'    Neutrosophic Generalized Pareto Distribution, \emph{Mathematics and Statistics},
#'    11(5), 827--833.
#' @importFrom stats runif dgamma pgamma qgamma
#' @examples
#' data(remission)
#' dnsgpd(x = remission, shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' pnsgpd(q = 20, shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' # Calculate quantiles
#' qnsgpd(p = c(0.25, 0.5, 0.75), shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#'
#' # Simulate 10 numbers
#' rnsgpd(n = 10, shape = c(1.1884, 1.1896), scale = c(7.6658, 7.7796))
#' @export
dnsgpd <- function(x, shape, scale) {
  if (any(shape <= 0) || any(scale <= 0) || any(x < 0)) {
    stop(message = "Arguments are incompatible.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  if (is.vector(x) || ncol(x) == 1) {
    x <- matrix(rep(as.numeric(x), each = 2), ncol = 2, byrow = TRUE)
  }

  if (ncol(x) > 2) {
    stop(message = "Arguments are incompatible.")
  }


  pdf <- matrix(NA, nrow = nrow(x), ncol = 2)
  for (i in 1:2) {
    pdf[, i] <- (1 / scale[i]) * (1 + (shape[i] * x[,i]) / scale[i])^(-1 / shape[i] - 1)
  }



  return(pdf)
}
#' @name Neutrosophic Generalized Pareto
#' @export
pnsgpd <- function(q, shape, scale, lower.tail = TRUE) {
  if (any(shape <= 0) || any(scale <= 0) || any(q < 0)) {
    stop(message = "incompatible arguments.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  if (is.vector(q) || ncol(q) == 1) {
    q <- matrix(rep(as.numeric(q), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(q) > 2) {
    stop(message = "Arguments are incompatible.")
  }



  cdf <- matrix(NA, nrow = nrow(q), ncol = 2)
  for (i in 1:2) {
    cdf[, i] <- 1 - (1 + (shape[i] * q[,i]) / scale[i])^(-1 / shape[i])
  }

  if (!lower.tail) {
    cdf <- 1 - cdf
  }




  return(cdf)
}

#' @name Neutrosophic Generalized Pareto
#' @export
qnsgpd <- function(p, shape, scale) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(shape <= 0) || any(scale <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  if (is.vector(p) || ncol(p) == 1) {
    p <- matrix(rep(as.numeric(p), each = 2), ncol = 2, byrow = TRUE)
  }
  if (ncol(p) > 2) {
    stop(message = "Arguments are incompatible.")
  }
  quantiles <- matrix(NA, nrow = nrow(p), ncol = 2)
  for (i in 1:2) {
    quantiles[, i] <- scale[i] / shape[i] * ((1 - p[, i])^(-shape[i]) - 1)
  }

  return(quantiles)
}

#' @name Neutrosophic Generalized Pareto
#' @export
rnsgpd <- function(n, shape, scale) {
  if (any(shape <= 0) || any(scale <= 0)) {
    stop(message = "Arguments are incompatible.")
  }

  shape <- rep(shape, length.out = 2)
  scale <- rep(scale, length.out = 2)

  X <- qnsgpd(runif(n), shape, scale)
  condition <- X[, 1] > X[, 2]
  X[condition, 1:2] <- X[condition, 2:1]

  return(X)
}
