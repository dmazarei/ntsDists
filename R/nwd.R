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
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param beta the value or vector lower and upper of the first shape parameter. must be positive.
#' @param alpha the value or vector lower and upper of the second scale parameter, must be positive.
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
#' p <- seq(0.1, 1, length.out = 21)
#' p2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pnwd(p, beta = 1, alpha = 2)
#' pnwd(p2, beta = c(1, 2), alpha = c(2, 2))
#' @export

pnwd <- function(p, beta, alpha) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(p < 0)) stop(message = "[Warning] 0 < p ")
  if (is.vector(p)) {
    F0 <- stats::pweibull(p, shape = beta[1], scale = alpha[1])
  } else {
    if (length(alpha) < 2 || length(beta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(p), ncol = 2)
      F0[, 1] <- stats::pweibull(p[, 1], shape = beta[1], scale = alpha[1])
      F0[, 2] <- stats::pweibull(p[, 2], shape = beta[2], scale = alpha[2])
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
    if (length(alpha) < 2 || length(beta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- stats::dweibull(x[, 1], shape = beta[1], scale = alpha[1])
      df[, 2] <- stats::dweibull(x[, 2], shape = beta[2], scale = alpha[2])
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
qnwd <- function(q, beta, alpha) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(alpha) < 2 || length(beta) < 2) {
    qf <- stats::qweibull(q, shape = beta[1], scale = alpha[1])
  } else {
    if (length(alpha) < 2 || length(beta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- stats::qweibull(q[, 1], shape = beta[1], scale = alpha[1])
      qf[, 2] <- stats::qweibull(q[, 2], shape = beta[2], scale = alpha[2])
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
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnwd(u, alpha, beta)
  }
  return(X)
}
