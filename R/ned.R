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
#' @param x a vector or matrix of quantiles at which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles at which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param theta the shape parameter, must be positive.
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
#' x2 <- matrix(seq(0.01, 1, length.out = 40), ncol = 2)
#' pned(x2, theta = c(2, 3))
#' @export

pned <- function(q, theta) {
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (any(q < 0)) stop(message = "[Warning] 0 < q ")
  if (is.vector(q)) {
    F0 <- 1 - exp(-q * theta[1])
  } else {
    if (length(theta) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- 1 - exp(-q[, i] * theta[i])
      }
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
    if (length(theta) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- theta[i] * exp(-x[, i] * theta[i])
      }
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
qned <- function(p, theta) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(theta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(theta) < 2) {
    qf <- -log(1 - p) / theta[1]
  } else {
    if (length(theta) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == ncol(p)) {
        p <- matrix(p, nrow = 1, ncol = ncol(p))
      }

      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- -log(1 - p[, i]) / theta[i]
      }
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
    u <- matrix(runif(n * length(theta)), nrow = n, ncol = length(theta))
    X <- qned(u, theta)
  }
  return(X)
}
