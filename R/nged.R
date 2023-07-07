#' Neutrosophic Generalized Exponential Distribution (NGED)
#'
#' Computes the pdf, cdf, hdf, shd, quantile and random numbers of the nuetrosophic generalized exponential distribution.
#' \deqn{f=\frac{\delta_N}{v_N}\left(1-\exp \left\{-\frac{x_N}{v_N}\right\}\right)^{\delta_N-1} \exp \left\{-\frac{x_N}{v_N}\right\}}
#' for   \eqn{v_N > 0}, the first shape parameter, and \eqn{\delta_N > 0}, the second shape parameter.
#' @name NGED
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param v the value or vector lower and upper of the first shape parameter, must be positive.
#' @param delta the value or vector lower and upper of the second shape parameter, must be positive.
#'
#' @return  \code{pnged} gives the distribution function,
#'  \code{dnged} gives the density,
#'  \code{qnged} gives the quantile function,
#'  \code{hnged} gives the hazard function,
#'  \code{snged} gives the survival function and
#'  \code{rnged} generates random variables from the Neutrosophic Generalized Exponential Distribution (NGED).
#' @references
#'    Rao, Gadde Srinivasa, Mina Norouzirad, and Danial Mazarei. "Neutrosophic Generalized Exponential Distribution with Application." Neutrosophic Sets and Systems 55.1 (2023): 28.
#' @importFrom stats runif
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pnged(x)
#' pnged(x, v = 2, delta = 1)
#' pnged(x2, v = c(1, 2), delta = c(2, 2))
#' @export

pnged <- function(x, v = 1, delta = 2) {
  if (any(v <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    F0 <- (1 - exp(-x / v[1]))^delta[1]
  } else {
    if (length(v) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- (1 - exp(-x[, 1] / v[1]))^delta[1]
      F0[, 2] <- (1 - exp(-x[, 2] / v[2]))^delta[2]
    }
  }
  return(F0)
}

#' @name NGED
#' @examples
#' dnged(x, v = 1, delta = 2)
#' curve(dnged, .1, 3)
#' dnged(x2, v = c(1, 2), delta = c(2, 2))
#' @export
dnged <- function(x, v = 1, delta = 2) {
  if (any(v <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    df <- (delta[1] / v[1]) * (1 - exp(-x / v[1]))^(delta[1] - 1) * exp(-x / v[1])
  } else {
    if (length(v) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- (delta[1] / v[1]) * (1 - exp(-x[, 1] / v[1]))^(delta[1] - 1) * exp(-x[, 1] / v[1])
      df[, 2] <- (delta[2] / v[2]) * (1 - exp(-x[, 2] / v[2]))^(delta[2] - 1) * exp(-x[, 2] / v[2])
    }
  }
  return(df)
}

#' @name NGED
#' @examples
#' qnged(x, v = 1, delta = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnged(x2, v = c(1, 2), delta = c(2, 2))
#' @export
qnged <- function(q, v = 1, delta = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(v <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(v) < 2 || length(delta) < 2) {
    qf <- log(1 - (q^(1 / delta[1]))) * (-v[1])
  } else {
    if (length(v) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- log(-q[, 1]^(1 / delta[1]) + 1) * (-v[1])
      qf[, 2] <- log(-q[, 2]^(1 / delta[2]) + 1) * (-v[2])
    }
  }
  return(qf)
}

#' @name NGED
#' @examples
#' n <- 10
#' rnged(n, v = 2, delta = 1)
#' rnged(n, v = c(1, 2), delta = c(1, 1))
#' @export
rnged <- function(n, v = 1, delta = 2) {
  if (any(v <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (length(v) < 2 || length(delta) < 2) {
    u <- runif(n)
    X <- qnged(u, v, delta)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnged(u, v, delta)
  }
  return(X)
}


#' @name NGED
#' @examples
#' hnged(x, v = 1, delta = 3)
#' curve(hnged, 0.1, 3)
#' hnged(x2, v = c(1, 2), delta = c(2, 2))
#' @export

hnged <- function(x, v = 1, delta = 2) {
  if (any(v <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")
  if (is.vector(x)) {
    h <- ((delta[1] / v[1]) * (1 - exp(-x / v[1]))^(delta[1] - 1) * exp(-x / v[1])) / (1 - (1 - exp(-x / v[1]))^delta[1])
  } else {
    if (length(v) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      h <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      h[, 1] <- ((delta[1] / v[1]) * (1 - exp(-x[, 1] / v[1]))^(delta[1] - 1) * exp(-x[, 1] / v[1])) / (1 - (1 - exp(-x[, 1] / v[1]))^delta[1])
      h[, 2] <- ((delta[2] / v[2]) * (1 - exp(-x[, 2] / v[2]))^(delta[2] - 1) * exp(-x[, 2] / v[2])) / (1 - (1 - exp(-x[, 2] / v[2]))^delta[2])
    }
  }
  return(h)
}

#' @name NGED
#' @examples
#' snged(x, v = 1, delta = 3)
#' curve(snged, 0.1, 3)
#' snged(x2, v = c(1, 2), delta = c(2, 2))
#' @export
snged <- function(x, v = 1, delta = 2) {
  if (any(v <= 0) || any(delta <= 0)) stop(message = "incompatible arguments.")
  if (any(x < 0)) stop(message = "[Warning] 0 < x ")

  if (is.vector(x)) {
    s <- 1 - (1 - exp(-x / v[1]))^delta[1]
  } else {
    if (length(v) < 2 || length(delta) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      s <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      s[, 1] <- 1 - (1 - exp(-x[, 1] / v[1]))^delta[1]
      s[, 2] <- 1 - (1 - exp(-x[, 2] / v[2]))^delta[2]
    }
  }
  return(s)
}
