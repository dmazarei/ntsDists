#' Neutrosophic Gamma Distribution (NGD)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic gamma distribution with parameter \code{shape}=\eqn{\alpha_N} and
#' \code{scale}=\eqn{\lambda_N}.
#'
#' The neutrosophic gamma distribution with parameters \code{shape}=\eqn{\alpha_N} and
#' \code{scale}=\eqn{\lambda_N} has density
#' \deqn{f_n(x)=\frac{1}{\Gamma(\alpha_n) \lambda_n^{\alpha_n}} x^{\alpha_n-1} e^{-\left(x / \lambda_n\right)}}
#' for \eqn{x \ge 0}, \eqn{\alpha_n > 0}, the shape parameter, and
#' \eqn{\lambda_n > 0}, the scale parameter. Here, \eqn{\Gamma(\cdot)} is gamma
#' function implemented by \code{\link{gamma}}.
#'
#' @name NGD
#' @param x scaler or vector or matrix lower and upper of values at which the
#' pdf or cdf needs to be computed.
#' @param p scaler or vector of probabilities at which the quantile needs to be computed.
#' @param q scaler or vector of quantiles.
#' @param n number of random numbers to be generated.
#' @param alpha the positive value or vector lower and upper of the first shape parameter.
#' @param lambda the positive value or vector lower and upper of the second scale parameter.
#'
#' @return
#'  \code{pngd} gives the distribution function,
#'  \code{dngd} gives the density,
#'  \code{qngd} gives the quantile function and
#'  \code{rngd} generates random variables from the neutrosophic gamma distribution.
#' @references
#'    Khan, Z., Al-Bossly, A., Almazah, M. M. A., and Alduais, F. S. (2021).
#'    On statistical development of neutrosophic gamma distribution with
#'    applications to complex data analysis, \emph{Complexity}, 2021, Article ID 3701236.
#' @importFrom stats runif dgamma pgamma qgamma
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' pngd(x, alpha = 2, lambda = 1)
#'
#' x2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pngd(x2, alpha = c(1, 2), lambda = c(2, 2))
#' @export

pngd <- function(p, alpha = 1, lambda = 2) {
  if (any(alpha <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    F0 <- stats::pgamma(p, shape = alpha[1], scale = lambda[1])
  } else {
    if (length(alpha) < 2 || length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- stats::pgamma(x[, 1], shape = alpha[1], scale = lambda[1])
      F0[, 2] <- stats::pgamma(x[, 2], shape = alpha[2], scale = lambda[2])
    }
  }
  return(F0)
}

#' @name NGD
#' @examples
#' dngd(x, alpha = 1, lambda = 2)
#'
#' dngd(x2, alpha = c(1, 2), lambda = c(2, 2))
#' @export
dngd <- function(x, alpha = 1, lambda = 2) {
  if (any(alpha <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- stats::dgamma(x, shape = alpha[1], scale = lambda[1])
  } else {
    if (length(alpha) < 2 || length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- stats::dgamma(x[, 1], shape = alpha[1], scale = lambda[1])
      df[, 2] <- stats::dgamma(x[, 2], shape = alpha[2], scale = lambda[2])
    }
  }
  return(df)
}

#' @name NGD
#' @examples
#' qngd(x, alpha = 1, lambda = 2)
#'
#' qngd(x2, alpha = c(1, 2), lambda = c(2, 2))
#' @export
qngd <- function(q, alpha = 1, lambda = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(alpha <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(alpha) < 2 || length(lambda) < 2) {
    qf <- stats::qgamma(q, shape = alpha[1], scale = lambda[1])
  } else {
    if (length(alpha) < 2 || length(lambda) < 2) {
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
#' rngd(n, alpha = 2, lambda = 1)
#' rngd(n, alpha = c(1, 2), lambda = c(1, 1))
#' @export
rngd <- function(n, alpha = 1, lambda = 2) {
  if (any(alpha <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (length(alpha) < 2 || length(lambda) < 2) {
    u <- runif(n)
    X <- qngd(u, alpha, lambda)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qngd(u, alpha, lambda)
  }
  return(X)
}
