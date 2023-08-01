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
#' @param x,q a vector or matrix of quantiles at which the pdf or cdf needs to be computed.
#' @param p a vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha shape parameter, must be positive.
#' @param lambda scale parameter, must be positive.
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
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pngd(x2, alpha = c(1, 2), lambda = c(2, 2))
#' @export

pngd <- function(q, alpha = 1, lambda = 2) {
  if (any(alpha <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (any(q < 0)) stop(message = "[Warning] 0 < q ")
  if (is.vector(q)) {
    F0 <- stats::pgamma(q, shape = alpha[1], scale = lambda[1])
  } else {
    if (length(alpha) < ncol(q) || length(lambda) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- stats::pgamma(q[, i], shape = alpha[i], scale = lambda[i])
      }
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
    if (length(alpha) < ncol(x) || length(lambda) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- stats::dgamma(x[, i], shape = alpha[i], scale = lambda[i])
      }
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
qngd <- function(p, alpha = 1, lambda = 2) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(alpha <= 0) || any(lambda <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(alpha) < 2 || length(lambda) < 2) {
    qf <- stats::qgamma(p, shape = alpha[1], scale = lambda[1])
  } else {
    if (length(alpha) < ncol(p) || length(lambda) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == ncol(p)) {
        p <- matrix(p, nrow = 1, ncol = ncol(p))
      }
      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- stats::qgamma(p[, i], shape = alpha[i], scale = lambda[i])
      }
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
    a <- min(length(alpha), length(lambda))
    u <- matrix(runif(n * a), nrow = n, ncol = a)
    X <- qngd(u, alpha[1:a], lambda[1:a])
  }
  return(X)
}
