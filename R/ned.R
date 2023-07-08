#' Neutrosophic Exponential Distribution (NED)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic exponential distribution with parameter \code{rate}=\eqn{\theta_N}.
#'
#' The neutrosophic exponential distribution with parameters \code{rate}=\eqn{\theta_N}
#' has density
#' \deqn{f_N(x)=\theta_N \exp \left(-x \theta_N\right)}
#' for \eqn{x \ge 0} and \eqn{\theta_N > 0}, the rate parameter.
#' @name NED
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param p scaler or vector of probabilities at which the quantile needs to be computed.
#' @param q scaler or vector of quantiles.
#' @param n number of random numbers to be generated.
#' @param theta the positive value or vector lower and upper of the first shape parameter.
#'
#' @return
#'  \code{pned} gives the distribution function,
#'  \code{dned} gives the density,
#'  \code{qned} gives the quantile function and
#'  \code{rned} generates random variables from the neutrosophic exponential distribution.
#' @references
#' Duan, W., Q., Khan, Z., Gulistan, M., Khurshid, A. (2021). Neutrosophic
#' Exponential Distribution: Modeling and Applications for Complex Data Analysis,
#' \emph{Complexity}, 2021, 1-8.
#'
#' @importFrom stats runif
#' @examples
#' x <- seq(0.01, 1, length.out = 21)
#' pned(x, theta = 1)
#'
#' x2 <- matrix(seq(0.01, 2, length.out = 40), ncol = 2)
#' pned(x2, theta = c(2, 3))
#' @export

pned <- function(p, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    F0 <- 1 - exp(-p * theta[1])
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- 1 - exp(-p[, 1] * theta[1])
      F0[, 2] <- 1 - exp(-p[, 2] * theta[2])
    }
  }
  return(F0)
}

#' @name NED
#' @examples
#'
#' dned(x, theta = 2)
#'
#' dned(x2, theta = c(2, 2))
#' @export
dned <- function(x, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- theta[1] * exp(-x * theta[1])
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- theta[1] * exp(-x[, 1] * theta[1])
      df[, 2] <- theta[2] * exp(-x[, 2] * theta[2])
    }
  }
  return(df)
}

#' @name NED
#' @examples
#'
#' qned(x, theta = 2)
#'
#' qned(x2, theta = c(2, 2))
#'
#' @export
qned <- function(q, theta) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(theta) < 2) {
    qf <- -log(1 - q) / theta[1]
  } else {
    if (length(theta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- -log(1 - q[, 1]) / theta[1]
      qf[, 2] <- -log(1 - q[, 2]) / theta[2]
    }
  }
  return(qf)
}

#' @name NED
#' @examples
#' n <- 10
#' rned(n, theta = 1)
#'
#' rned(n, theta = c(1, 2))
#' @export
rned <- function(n, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (length(theta) < 2) {
    u <- runif(n)
    X <- qned(u, theta)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qned(u, theta)
  }
  return(X)
}
