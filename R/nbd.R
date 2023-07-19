#' Neutrosophic Beta Distribution (NBD)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic Beta distribution with parameters \code{shape1}=\eqn{\alpha_N} and
#' \code{shape2}=\eqn{\beta_N}.
#'
#' The neutrosophic beta distribution with parameters \code{shape1}=\eqn{\alpha_N} and
#' \code{shape2}=\eqn{\beta_N} has density
#' \deqn{f_N(X)=\frac{1}{B\left(\alpha_N, \beta_N\right)}X^{\alpha_N-1}(1-X)^{\beta_N-1}}
#' for  \eqn{\alpha_N > 0}, the first shape parameter, and \eqn{\beta_N > 0},
#' the second shape parameter and \eqn{0 \le x \le 1}. The function \eqn{B(a,b)}
#' return the beta function that can be calculated by \code{\link{beta}}.
#'
#' @name NBD
#'
#' @param x,q vector or matrix  of quantiles at which the pdf or cdf needs to be computed.
#' @param p vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the positive value or vector lower and upper of the first shape parameter.
#' @param beta the positive value or vector lower and upper of the second shape parameter.
#'
#' @return
#' \code{pnbd} gives the distribution function,
#' \code{dnbd} gives the density,
#' \code{qnbd} gives the quantile function and
#' \code{rnbd} generates random variables from the neutrosophic Beta distribution.
#'
#' @references
#'    Sherwani, R. Ah. K., Naeem, M., Aslam, M., Reza, M. A., Abid, M., Abbas, S. (2021).
#'     Neutrosophic beta distribution with properties and applications. Neutrosophic Sets and Systems, 41, 209-214.
#' @importFrom stats runif dbeta pbeta qbeta
#' @examples
#' x <- seq(0.1, 1, length.out = 21) #'
#' pnbd(x, alpha = 2, beta = 1)
#'
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnbd(x2, alpha = c(1, 2), beta = c(2, 2))
#'
#' @export

pnbd <- function(q, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(q < 0)) stop(message = "[Warning] 0 < q ")
  if (is.vector(q)) {
    F0 <- stats::pbeta(q, shape1 = alpha[1], shape2 = beta[1])
  } else {
    if (length(alpha) < ncol(q) || length(beta) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- stats::pbeta(q[, i], shape1 = alpha[i], shape2 = beta[i])
      }
    }
  }
  return(F0)
}

#' @name NBD
#' @examples
#' dnbd(x, alpha = 1, beta = 2)
#'
#' dnbd(x2, alpha = c(1, 2), beta = c(2, 2))
#'
#' @export
dnbd <- function(x, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- stats::dbeta(x, shape1 = alpha[1], shape2 = beta[1])
  } else {
    if (length(alpha) < ncol(x) || length(beta) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- stats::dbeta(x[, i], shape1 = alpha[i], shape2 = beta[i])
      }
    }
  }
  return(df)
}

#' @name NBD
#' @examples
#'
#' qnbd(x, alpha = 1, beta = 2)
#' qnbd(x2, alpha = c(1, 2), beta = c(2, 2))
#'
#' @export
qnbd <- function(p, alpha, beta) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(alpha) < 2 || length(beta) < 2) {
    qf <- stats::qbeta(p, shape1 = alpha[1], shape2 = beta[1])
  } else {
    if (length(alpha) < ncol(p) || length(beta) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == 2) {
        p <- matrix(p, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- stats::qbeta(p[, i], shape1 = alpha[i], shape2 = beta[i])
      }
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
rnbd <- function(n, alpha, beta) {
  if (any(alpha <= 0) || any(beta <= 0)) stop(message = "incompatible arguments.")
  if (length(alpha) < 2 || length(beta) < 2) {
    u <- runif(n)
    X <- qnbd(u, alpha, beta)
  } else {
    a <- min(length(alpha), length(beta))
    u <- matrix(runif(n * a), nrow = n, ncol = a)
    X <- qnbd(u, alpha[1:a], beta[1:a])
  }
  return(X)
}
