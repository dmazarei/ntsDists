#' Neutrosophic Weibull Distribution (NWD)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic Weibull distribution with scale parameter \eqn{\alpha_N} and
#' shape parameter \eqn{\beta_N}.
#'
#' The neutrosophic Rayleigh distribution with parameters \eqn{\alpha_N} and
#' \eqn{\beta_N} has the density
#' \deqn{f_N(x)=\frac{\beta_N}{\alpha_N^{\beta_N}} x^{\beta_N-1}
#'     e^{-\left(x / \alpha_N\right)^{\beta_N}}}
#' for \eqn{x > 0}, \eqn{\beta_N > 0}, shape parameter, and \eqn{\alpha_N > 0},
#' scale parameter.
#'
#' @name NWD
#' @param x,q vector or matrix of quantiles at which the pdf or cdf needs to be computed.
#' @param p vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param beta the value or vector of the first shape parameter. must be positive.
#' @param alpha the value or vector of the second scale parameter, must be positive.
#'
#' @return
#'  \code{pnwd} gives the distribution function,
#'  \code{dnwd} gives the density,
#'  \code{qnwd} gives the quantile function and
#'  \code{rnwd} generates random variables from the Neutrosophic Weibull Distribution (NWD).
#' @references
#'    Alhasan, K. F. H. and Smarandache, F. (2019). Neutrosophic Weibull
#'    distribution and Neutrosophic Family Weibull Distribution,
#'    \emph{Neutrosophic Sets and Systems}, 28, 191-199.
#'
#' @importFrom stats runif dweibull pweibull qweibull
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnwd(x, beta = 1, alpha = 2)
#' pnwd(x2, beta = c(1, 2), alpha = c(2, 2))
#' @export

pnwd <- function(q, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(q < 0)) stop(message = "[Warning] 0 < q ")
  if (is.vector(q)) {
    F0 <- stats::pweibull(q, shape = beta[1], scale = alpha[1])
  } else {
    if (length(alpha) < ncol(q) || length(beta) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- stats::pweibull(q[, i], shape = beta[i], scale = alpha[i])
      }
    }
  }
  return(F0)
}

#' @name NWD
#' @examples
#' dnwd(x, beta = 1, alpha = 2)
#'
#' dnwd(x2, beta = c(1, 2), alpha = c(2, 2))
#' @export
dnwd <- function(x, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- stats::dweibull(x, shape = beta[1], scale = alpha[1])
  } else {
    if (length(alpha) < ncol(x) || length(beta) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- stats::dweibull(x[, i], shape = beta[i], scale = alpha[i])
      }
    }
  }
  return(df)
}

#' @name NWD
#' @examples
#' qnwd(x, beta = 1, alpha = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnwd(x2, beta = c(1, 2), alpha = c(2, 2))
#' @export
qnwd <- function(p, beta, alpha) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(alpha) < 2 || length(beta) < 2) {
    qf <- stats::qweibull(p, shape = beta[1], scale = alpha[1])
  } else {
    if (length(alpha) < ncol(p) || length(beta) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == ncol(p)) {
        p <- matrix(p, nrow = 1, ncol = ncol(p))
      }
      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- stats::qweibull(p[, i], shape = beta[i], scale = alpha[i])
      }
    }
  }
  return(qf)
}

#' @name NWD
#' @examples
#' n <- 10
#' rnwd(n, beta = 1, alpha = 2)
#' rnwd(n, beta = c(1, 2), alpha = c(1, 1))
#' @export
rnwd <- function(n, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (length(alpha) < 2 || length(beta) < 2) {
    u <- runif(n)
    X <- qnwd(u, alpha, beta)
  } else {
    a <- min(length(alpha), length(beta))
    u <- matrix(runif(n * a), nrow = n, ncol = a)
    X <- qnwd(u, alpha[1:a], beta[1:a])
  }
  return(X)
}
