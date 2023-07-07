#' Neutrosophic Rayleigh Distribution (NRD)
#'
#' Computes the pdf, cdf, quantile and random numbers of the
#' \deqn{f_N(z)=\frac{Z}{\theta_N^2} e^{-\frac{1}{2}\left(\frac{Z}{\theta_N}\right)^2}}
#' for   \eqn{\theta_N > 0}, the first shape parameter.
#' @name NRD
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param theta the value or vector lower and upper of the first shape parameter, must be positive.
#'
#' @return  \code{pnrd} gives the distribution function,
#'  \code{dnrd} gives the density,
#'  \code{qnrd} gives the quantile function and
#'  \code{rnrd} generates random variables from the Neutrosophic Rayleigh Distribution (NRD).
#' @references
#'        Khan, Zahid, et al. "Neutrosophic Rayleigh model with some basic characteristics and engineering applications." IEEE Access 9 (2021): 71277-71283.
#' @importFrom stats runif
#' @examples
#' x <- seq(0.01, 1, length.out = 21)
#' x2 <- matrix(seq(0.01, 2, length.out = 40), ncol = 2)
#' pnrd(x)
#' pnrd(x, theta = 1)
#' pnrd(x2, theta = c(2, 3))
#' @export

pnrd <- function(x, theta = 2) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    F0 <- 1 - exp((-1 / 2) * (x / theta[1])^2)
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- 1 - exp((-1 / 2) * (x[, 1] / theta[1])^2)
      F0[, 2] <- 1 - exp((-1 / 2) * (x[, 2] / theta[2])^2)
    }
  }
  return(F0)
}

#' @name NRD
#' @examples
#' dnrd(x, theta = 2)
#' curve(dnrd, .1, 3)
#' dnrd(x2, theta = c(2, 2))
#' @export
dnrd <- function(x, theta = 2) {
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
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnrd(x2, theta = c(2, 2))
#' @export
qnrd <- function(q, theta = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(theta) < 2) {
    qf <- theta[1] * sqrt(-2 * log(1 - q))
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- theta[1] * sqrt(-2 * log(1 - q[, 1]))
      qf[, 2] <- theta[2] * sqrt(-2 * log(1 - q[, 2]))
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
rnrd <- function(n, theta = 2) {
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
