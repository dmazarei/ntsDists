#' Neutrosophic Poisson Distribution (NPD)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic Poisson distribution with parameter \eqn{\lambda_N}.
#'
#' The neutrosophic Poisson distribution with parameter \eqn{\theta} has the density
#' \deqn{f_N(x)=e^{-\lambda_N}  \frac{\left(\lambda_N\right)^x}{x !}}
#' for   \eqn{\lambda_N > 0}.
#'
#' @name NPD
#' @param x,q a vector or matrix of (non-negative integer) quantiles at which the pdf or cdf needs to be computed.
#' @param p a vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param lambda the mean that must be positive.
#'
#' @return
#'  \code{pnpd} gives the distribution function,
#'  \code{dnpd} gives the density,
#'  \code{qnpd} gives the quantile function and
#'  \code{rnpd} generates random variables from the Neutrosophic Poisson Distribution (NPD).
#' @references
#'        Alhabib, R., Ranna, M. M., Farah, H., Salama, A. A. (2018).
#'        Some neutrosophic probability distributions.
#'        \emph{Neutrosophic Sets and Systems},  22, 30-38.
#' @importFrom stats runif dpois ppois qpois
#' @examples
#' x <- 1:10
#' x2 <- matrix(1:20, ncol = 2)
#' pnpd(x, lambda = 1)
#' pnpd(x2, lambda = c(2, 3))
#' @export

pnpd <- function(q, lambda) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (any(q < 0) && any(q - floor(q) == 0)) stop(message = "[Warning] 0 < q or non-integer ")
  if (is.vector(q)) {
    F0 <- stats::ppois(q, lambda = lambda[1])
  } else {
    if (length(lambda) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- stats::ppois(q[, i], lambda = lambda[i])
      }
    }
  }
  return(F0)
}

#' @name NPD
#' @examples
#' dnpd(x, lambda = 2)
#' dnpd(x2, lambda = c(2, 2))
#' @export
dnpd <- function(x, lambda) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (any(x < 0) && any(x - floor(x) == 0)) stop(message = "[Warning] 0 < x or non-integer ")
  if (is.vector(x)) {
    df <- stats::ppois(x, lambda = lambda[1])
  } else {
    if (length(lambda) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- stats::dpois(x[, i], lambda = lambda[i])
      }
    }
  }
  return(df)
}

#' @name NPD
#' @examples
#' q1 <- seq(0.1, 1, length.out = 40)
#' qnpd(q1, lambda = 2)
#' q2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnpd(q2, lambda = c(2, 2))
#' @export
qnpd <- function(p, lambda) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(lambda) < 2) {
    qf <- stats::qpois(p, lambda = lambda[1])
  } else {
    if (length(lambda) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == 2) {
        p <- matrix(p, nrow = 1, ncol = 2)
      }
      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- stats::qpois(p[, i], lambda = lambda[i])
      }
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
rnpd <- function(n, lambda) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (length(lambda) < 2) {
    u <- runif(n)
    X <- qnpd(u, lambda)
  } else {
    u <- matrix(runif(n * length(lambda)), nrow = n, ncol = length(lambda))
    X <- qnpd(u, lambda)
  }
  return(X)
}
