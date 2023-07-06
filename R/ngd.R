#' Neutrosophic Gamma Distribution (NGD)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the nuetrosophic generalized exponential distribution.
#'\deqn{f\xi_n(z)=\frac{1}{\Gamma p_n \lambda_n^{p_n}} z^{p_n-1} e^{-\left(z / \lambda_n\right)}}
#' for   \eqn{p_n > 0}, the first shape parameter, and \eqn{\lambda_n > 0}, the second scale parameter.
#' @name NGD
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param p the value or vector lower and upper of the first shape parameter, must be positive.
#' @param lambda the value or vector lower and upper of the second scale parameter, must be positive.
#'
#' @return  \code{pngd} gives the distribution function,
#'  \code{dngd} gives the density,
#'  \code{qngd} gives the quantile function and
#'  \code{rngd} generates random variables from the Neutrosophic Gamma Distribution (NGD).
#' @references
#'    Khan, Zahid, et al. "On statistical development of neutrosophic gamma distribution with applications to complex data analysis." Complexity 2021 (2021): 1-8.
#' @importFrom stats runif dgamma pgamma qgamma
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' x2 <- matrix(seq(0, 2, length.out = 40), ncol = 2)
#' pngd(x)
#' pngd(x, p = 2, lambda = 1)
#' pngd(x2, p = c(1, 2), lambda = c(2, 2))
#' @export

pngd <- function(x, p = 1, lambda = 2) {
  if (any(p <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (any(x <= 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    F0 <- stats::pgamma(x, shape = p[1], scale = lambda[1])
  } else {
    if (length(p) < 2 || length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- stats::pgamma(x[, 1], shape = p[1], scale = lambda[1])
      F0[, 2] <- stats::pgamma(x[, 2], shape = p[2], scale = lambda[2])
    }
  }
  return(F0)
}

#' @name NGD
#' @examples
#' dngd(x, p = 1, lambda = 2)
#' curve(dngd, .1, 3)
#' dngd(x2, p = c(1, 2), lambda = c(2, 2))
#' @export
dngd <- function(x, p = 1, lambda = 2) {
  if (any(p <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (any(x <= 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- stats::dgamma(x, shape = p[1], scale = lambda[1])
  } else {
    if (length(p) < 2 || length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- stats::dgamma(x[, 1], shape = p[1], scale = lambda[1])
      df[, 2] <- stats::dgamma(x[, 2], shape = p[2], scale = lambda[2])
    }
  }
  return(df)
}

#' @name NGD
#' @examples
#' qngd(x, p = 1, lambda = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qngd(x2, p = c(1, 2), lambda = c(2, 2))
#' @export
qngd <- function(q, p = 1, lambda = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(p <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(p) < 2 || length(lambda) < 2) {
    qf <- stats::qgamma(q, shape = p[1], scale = lambda[1])
  } else {
    if (length(p) < 2 || length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- stats::qgamma(q[, 1], shape = p[1], scale = lambda[1])
      qf[, 2] <- stats::qgamma(q[, 2], shape = p[2], scale = lambda[2])
    }
  }
  return(qf)
}

#' @name NGD
#' @examples
#' n <- 10
#' rngd(n, p = 2, lambda = 1)
#' rngd(n, p = c(1, 2), lambda = c(1, 1))
#' @export
rngd <- function(n, p = 1, lambda = 2) {
  if (any(p <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (length(p) < 2 || length(lambda) < 2) {
    u <- runif(n)
    X <- qngd(u, p, lambda)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qngd(u, p, lambda)
  }
  return(X)
}
