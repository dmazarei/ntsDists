#' Neutrosophic Poisson Distribution (NPD)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the
#' \deqn{f_N(x)=e^{-\lambda_N} \cdot \frac{\left(\lambda_N\right)^x}{x !}}
#' for   \eqn{\lambda_N > 0}, the first shape parameter.
#' @name NPD
#' @param x scaler or vector or matrix lower and upper of values (non-negative integer) at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param lambda the value or vector lower and upper of (non-negative)the first shape parameter, must be positive.
#'
#' @return  \code{pnpd} gives the distribution function,
#'  \code{dnpd} gives the density,
#'  \code{qnpd} gives the quantile function and
#'  \code{rnpd} generates random variables from the Neutrosophic Poisson Distribution (NPD).
#' @references
#'        Alhabib, R., et al., Some neutrosophic probability distributions. Neutrosophic Sets and Systems, 2018. 22: p.30-38.
#' @importFrom stats runif dpois ppois qpois
#' @examples
#' x <- 1:10
#' x2 <- matrix(1:20, ncol = 2)
#' pnpd(x)
#' pnpd(x, lambda = 1)
#' pnpd(x2, lambda = c(2, 3))
#' @export

pnpd <- function(x, lambda = 2) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (any(x < 0) && any(x - floor(x) == 0)) stop(message = "[Warning] 0 < x or non-integer ")
  if (is.vector(x)) {
    F0 <- stats::ppois(x, lambda = lambda[1])
  } else {
    if (length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- stats::ppois(x[, 1], lambda = lambda[1])
      F0[, 2] <- stats::ppois(x[, 2], lambda = lambda[2])
    }
  }
  return(F0)
}

#' @name NPD
#' @examples
#' dnpd(x, lambda = 2)
#' curve(dnpd, .1, 3)
#' dnpd(x2, lambda = c(2, 2))
#' @export
dnpd <- function(x, lambda = 2) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (any(x < 0) && any(x - floor(x) == 0)) stop(message = "[Warning] 0 < x or non-integer ")
  if (is.vector(x)) {
    df <- stats::ppois(x, lambda = lambda[1])
  } else {
    if (length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- stats::dpois(x[, 1], lambda = lambda[1])
      df[, 2] <- stats::dpois(x[, 2], lambda = lambda[2])
    }
  }
  return(df)
}

#' @name NPD
#' @examples
#' x3 <- seq(0.1, 1, length.out = 40)
#' qnpd(x3, lambda = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnpd(x2, lambda = c(2, 2))
#' @export
qnpd <- function(q, lambda = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(lambda) < 2) {
    qf <- stats::qpois(q, lambda = lambda[1])
  } else {
    if (length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- stats::qpois(q[, 1], lambda = lambda[1])
      qf[, 2] <- stats::qpois(q[, 2], lambda = lambda[2])
    }
  }
  return(qf)
}

#' @name NPD
#' @examples
#' n <- 10
#' rnpd(n, lambda = 1)
#' rnpd(n, lambda = c(1, 2))
#' @export
rnpd <- function(n, lambda = 2) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (length(lambda) < 2) {
    u <- runif(n)
    X <- qnpd(u, lambda)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnpd(u, lambda)
  }
  return(X)
}
