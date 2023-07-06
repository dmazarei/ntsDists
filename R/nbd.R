#' Neutrosophic Beta Distribution (NBD)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the nuetrosophic generalized exponential distribution.
#'\deqn{f_N(X)=\frac{X^{\alpha_N-1}(1-X)^{\beta_N-1}}{\beta\left(\alpha_N, \beta_N\right)}}
#' for   \eqn{\alpha_N > 0}, the first shape parameter, and \eqn{\beta_N > 0}, the second shape parameter.
#' @name NBD
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value or vector lower and upper of the first shape parameter. must be positive.
#' @param beta the value or vector lower and upper of the second shape parameter, must be positive.
#'
#' @return  \code{pnbd} gives the distribution function,
#'  \code{dnbd} gives the density,
#'  \code{qnbd} gives the quantile function and
#'  \code{rnbd} generates random variables from the Neutrosophic Beta Distribution (NBD).
#' @references
#'    Sherwani, Rehan A. Khan, et al. Neutrosophic beta distribution with properties and applications. Infinite Study, 2021.
#' @importFrom stats runif dbeta pbeta qbeta
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pnbd(x)
#' pnbd(x, alpha = 2, beta = 1)
#' pnbd(x2, alpha = c(1, 2), beta = c(2, 2))
#' @export

pnbd <- function(x, alpha = 1, beta = 2) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    F0 <- stats::pbeta(x, shape1 = alpha[1], shape2 = beta[1])
  } else {
    if (length(alpha) < 2 || length(beta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- stats::pbeta(x[, 1], shape1 = alpha[1], shape2 = beta[1])
      F0[, 2] <- stats::pbeta(x[, 2], shape1 = alpha[2], shape2 = beta[2])
    }
  }
  return(F0)
}

#' @name NBD
#' @examples
#' dnbd(x, alpha = 1, beta = 2)
#' curve(dnbd, .1, 3)
#' dnbd(x2, alpha = c(1, 2), beta = c(2, 2))
#' @export
dnbd <- function(x, alpha = 1, beta = 2) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- stats::dbeta(x, shape1 = alpha[1], shape2 = beta[1])
  } else {
    if (length(alpha) < 2 || length(beta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- stats::dbeta(x[, 1], shape1 = alpha[1], shape2 = beta[1])
      df[, 2] <- stats::dbeta(x[, 2], shape1 = alpha[2], shape2 = beta[2])
    }
  }
  return(df)
}

#' @name NBD
#' @examples
#' qnbd(x, alpha = 1, beta = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnbd(x2, alpha = c(1, 2), beta = c(2, 2))
#' @export
qnbd <- function(q, alpha = 1, beta = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(alpha) < 2 || length(beta) < 2) {
    qf <- stats::qbeta(q, shape1 = alpha[1], shape2 = beta[1])
  } else {
    if (length(alpha) < 2 || length(beta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- stats::qbeta(q[, 1], shape1 = alpha[1], shape2 = beta[1])
      qf[, 2] <- stats::qbeta(q[, 2], shape1 = alpha[2], shape2 = beta[2])
    }
  }
  return(qf)
}

#' @name NBD
#' @examples
#' n <- 10
#' rnbd(n, alpha = 2, beta = 1)
#' rnbd(n, alpha = c(1, 2), beta = c(1, 1))
#' @export
rnbd <- function(n, alpha = 1, beta = 2) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (length(alpha) < 2 || length(beta) < 2) {
    u <- runif(n)
    X <- qnbd(u, alpha, beta)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnbd(u, alpha, beta)
  }
  return(X)
}
