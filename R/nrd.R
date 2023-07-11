#' Neutrosophic Rayleigh Distribution (NRD)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic Rayleigh distribution with parameter \eqn{\theta_N}.
#'
#' The neutrosophic Rayleigh distribution with parameter \eqn{\theta} has the density
#' \deqn{f_N(x)=\frac{x}{\theta_N^2} e^{-\frac{1}{2}\left(\frac{x}{\theta_N}\right)^2}}
#' for  \eqn{\theta_N > 0}.
#'
#' @name NRD
#' @param x,q vector or matrix lower and upper of quantiles at which the pdf or cdf needs to be computed.
#' @param p vector or matrix lower and upper of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param theta the value or vector lower and upper of the first shape parameter, must be positive.
#'
#' @return
#'  \code{pnrd} gives the distribution function,
#'  \code{dnrd} gives the density,
#'  \code{qnrd} gives the quantile function and
#'  \code{rnrd} generates random variables from the Neutrosophic Rayleigh Distribution (NRD).
#' @references
#' Khan, Z., Gulistan, M., Kausar, N. and Park, C. (2021). Neutrosophic Rayleigh Model With Some Basic Characteristics and Engineering Applications, in \emph{IEEE Access}, 9, 71277-71283.
#'
#' @importFrom stats runif
#' @examples
#' x <- seq(0.01, 1, length.out = 21)
#' pnrd(x, theta = 1)
#'
#' x2 <- matrix(seq(0.01, 1, length.out = 40), ncol = 2)#'
#' pnrd(x2, theta = c(2, 3))
#' @export

pnrd <- function(q, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (any(q < 0)) stop(message = "[Warning] 0 < q ")
  if (is.vector(q)) {
    F0 <- 1 - exp((-1 / 2) * (q / theta[1])^2)
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      F0[, 1] <- 1 - exp((-1 / 2) * (q[, 1] / theta[1])^2)
      F0[, 2] <- 1 - exp((-1 / 2) * (q[, 2] / theta[2])^2)
    }
  }
  return(F0)
}

#' @name NRD
#' @examples
#' dnrd(x, theta = 2)
#' dnrd(x2, theta = c(2, 2))
#' @export
dnrd <- function(x, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- (x / theta[1]^2) * exp((-1 / 2) * (x / theta[1])^2)
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- (x[, 1] / theta[1]^2) * exp((-1 / 2) * (x[, 1] / theta[1])^2)
      df[, 2] <- (x[, 2] / theta[2]^2) * exp((-1 / 2) * (x[, 2] / theta[2])^2)
    }
  }
  return(df)
}

#' @name NRD
#' @examples
#' qnrd(x, theta = 2)
#'
#' qnrd(x2, theta = c(2, 2))
#' @export
qnrd <- function(p, theta) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(theta) < 2) {
    qf <- theta[1] * sqrt(-2 * log(1 - p))
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == 2) {
        p <- matrix(p, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(p), ncol = 2)
      qf[, 1] <- theta[1] * sqrt(-2 * log(1 - p[, 1]))
      qf[, 2] <- theta[2] * sqrt(-2 * log(1 - p[, 2]))
    }
  }
  return(qf)
}

#' @name NRD
#' @examples
#' n <- 10
#' rnrd(n, theta = 1)
#' rnrd(n, theta = c(1, 2))
#' @export
rnrd <- function(n, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (length(theta) < 2) {
    u <- runif(n)
    X <- qnrd(u, theta)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnrd(u, theta)
  }
  return(X)
}
