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
#' @param x,q vector or matrix lower and upper of (non-negative integer) quantiles at which the pdf or cdf needs to be computed.
#' @param p vector or matrix lower and upper of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param lambda the value or vector lower and upper of (non-negative)the first shape parameter, must be positive.
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
#' p <- 1:10
#' p2 <- matrix(1:20, ncol = 2)
#' pnpd(p, lambda = 1)
#' pnpd(p2, lambda = c(2, 3))
#' @export

pnpd <- function(q, lambda) {
  if (any(lambda < 0)) stop(message = "incompatible arguments.")
  if (any(q < 0) && any(q - floor(q) == 0)) stop(message = "[Warning] 0 < q or non-integer ")
  if (is.vector(q)) {
    F0 <- stats::ppois(q, lambda = lambda[1])
  } else {
    if (length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      F0[, 1] <- stats::ppois(q[, 1], lambda = lambda[1])
      F0[, 2] <- stats::ppois(q[, 2], lambda = lambda[2])
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
    if (length(lambda) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == 2) {
        p <- matrix(p, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(p), ncol = 2)
      qf[, 1] <- stats::qpois(p[, 1], lambda = lambda[1])
      qf[, 2] <- stats::qpois(p[, 2], lambda = lambda[2])
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
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnpd(u, lambda)
  }
  return(X)
}
