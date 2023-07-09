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
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param p scaler or vector of probabilities at which the quantile needs to be computed.
#' @param q  scaler or vector of quantiles.
#' @param n number of random numbers to be generated.
#' @param nu the positive value or vector lower and upper of the scale parameter.
#' @param delta the positive value or vector lower and upper of the shape parameter.
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
#' p <- seq(0.1, 1, length.out = 21)
#' pnged(p, nu = 2, delta = 1)
#'
#' p2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pnged(p2, nu = c(1, 2), delta = c(2, 2))
#' @export

pnged <- function(p, nu, delta) {
  if (any(nu <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(p < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(p)) {
    F0 <- (1 - exp(-p / nu[1]))^delta[1]
  } else {
    if (length(nu) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(p), ncol = 2)
      F0[, 1] <- (1 - exp(-p[, 1] / nu[1]))^delta[1]
      F0[, 2] <- (1 - exp(-p[, 2] / nu[2]))^delta[2]
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
    if (length(nu) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- (delta[1] / nu[1]) * (1 - exp(-x[, 1] / nu[1]))^(delta[1] - 1) * exp(-x[, 1] / nu[1])
      df[, 2] <- (delta[2] / nu[2]) * (1 - exp(-x[, 2] / nu[2]))^(delta[2] - 1) * exp(-x[, 2] / nu[2])
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
qnged <- function(q, nu, delta) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(nu <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(nu) < 2 || length(delta) < 2) {
    qf <- log(1 - (q^(1 / delta[1]))) * (-nu[1])
  } else {
    if (length(nu) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- log(-q[, 1]^(1 / delta[1]) + 1) * (-nu[1])
      qf[, 2] <- log(-q[, 2]^(1 / delta[2]) + 1) * (-nu[2])
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
   u <- matrix(runif(n * 2), nrow = n, ncol = 2)
   X <- qnged(u, nu, delta)
 }
 return(X)
}
