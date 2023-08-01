#' Neutrosophic Generalized Exponential Distribution (NGED)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic generalized exponential distribution with shape parameter
#' \eqn{\delta} and scale parameter \eqn{\nu}.
#'
#' The neutrosophic generalized exponential distribution with parameters
#' \code{shape}=\eqn{\delta} and \code{scale}=\eqn{\nu} has density
#' \deqn{f_n(x)=\frac{\delta_N}{\nu_N}\left(1-\exp \left\{-\frac{x_N}{\nu_N}\right\}\right)^{\delta_N-1} e^{\left\{-\frac{x_N}{\nu_N}\right\}}}
#' for\eqn{x \ge 0}, \eqn{\delta_N > 0}, the shape parameter,
#' and \eqn{\nu_N > 0}, the scale parameter.
#'
#' @name NGED
#' @param x,q quantiles at which the pdf or cdf needs to be computed.
#' @param p probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param nu scale parameter, must be positive.
#' @param delta  shape parameter, must be positive.
#'
#' @return
#'  \code{pnged} gives the distribution function,
#'  \code{dnged} gives the density,
#'  \code{qnged} gives the quantile function and
#'  \code{rnged} generates random variables
#'  from the neutrosophic generalized exponential distribution.
#'
#' @references
#'    Rao, G. S., Norouzirad, M., and Mazarei . D. (2023). Neutrosophic
#'    Generalized Exponential Distribution with Application.
#'    \emph{Neutrosophic Sets and Systems}, 55, 471-485.
#'
#' @importFrom stats runif
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' pnged(x, nu = 2, delta = 1)
#'
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnged(x2, nu = c(1, 2), delta = c(2, 2))
#' @export

pnged <- function(q, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(q < 0)) stop(message = "[Warning] 0 < q ")
  if (is.vector(q)) {
    F0 <- (1 - exp(-q / nu[1]))^delta[1]
  } else {
    if (length(nu) < ncol(q) || length(delta) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- (1 - exp(-q[, i] / nu[i]))^delta[i]
      }
    }
  }
  return(F0)
}

#' @name NGED
#' @examples
#' dnged(x, nu = 1, delta = 2)
#'
#' dnged(x2, nu = c(1, 2), delta = c(2, 2))
#' @export
dnged <- function(x, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- (delta[1] / nu[1]) * (1 - exp(-x / nu[1]))^(delta[1] - 1) * exp(-x / nu[1])
  } else {
    if (length(nu) < ncol(x) || length(delta) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- (delta[i] / nu[i]) * (1 - exp(-x[, i] / nu[i]))^(delta[i] - 1) * exp(-x[, i] / nu[i])
      }
    }
  }
  return(df)
}

#' @name NGED
#' @examples
#' qnged(x, nu = 1, delta = 2)
#'
#' qnged(x2, nu = c(1, 2), delta = c(2, 2))
#' @export
qnged <- function(p, nu, delta) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(nu <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(nu) < 2 || length(delta) < 2) {
    qf <- log(1 - (p^(1 / delta[1]))) * (-nu[1])
  } else {
    if (length(nu) < ncol(p) || length(delta) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == ncol(p)) {
        p <- matrix(p, nrow = 1, ncol = ncol(p))
      }

      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- log(-p[, i]^(1 / delta[i]) + 1) * (-nu[i])
      }
    }
  }
  return(qf)
}


#' @name NGED
#' @examples
#' n <- 10
#' rnged(n, nu = 2, delta = 1)
#' rnged(n, nu = c(1, 2), delta = c(1, 1))
#' @export
rnged <- function(n, nu = 1, delta = 2) {
  if (any(nu <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (length(nu) < 2 || length(delta) < 2) {
    u <- runif(n)
    X <- qnged(u, nu, delta)
  } else {
    a <- min(length(nu), length(delta))
    u <- matrix(runif(n * a), nrow = n, ncol = a)
    X <- qnged(u, nu[1:a], delta[1:a])
  }
  return(X)
}
