#' Neutrosophic Normal Distribution (NND)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic generalized exponential distribution with parameters mean,
#' \eqn{\mu}, and standard deviation \eqn{\sigma}.
#'
#' The neutrosophic normal distribution with parameters with \code{mean}=\eqn{\mu_N}
#' and \code{std}=\eqn{\sigma_N} has density
#' \deqn{f_N(x) = \frac{1}{\sigma_N \sqrt{2 \pi}} e^{\left(\frac{\left(X-\mu_N\right)^2}{2 \sigma_N^2}\right)}}
#' for \eqn{-\infty < x < \infty}, \eqn{\mu_N}, the mean, and \eqn{\sigma_N > 0},
#' the standard deviation.
#'
#' @name NND
#' @param x,q vector or matrix of quantiles at which the pdf or cdf needs to be computed.
#' @param p vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param mu the value or vector of the first parameter.
#' @param sigma the positive value or vector of the second parameter.
#'
#' @return
#' \code{pnnd} gives the distribution function,
#'  \code{dnnd} gives the density,
#'  \code{qnnd} gives the quantile function and
#'  \code{rnnd} generates random variables from the neutrosophic normal distribution.
#'
#' @references
#'    Patro, S. and F. Smarandache, The Neutrosophic Statistical Distribution, More Problems, More Solutions. 2016: Infinite Study.
#' @importFrom stats runif dnorm pnorm qnorm
#'
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' pnnd(x, mu = c(2, 2), sigma = c(1, 1))
#'
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' pnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export

pnnd <- function(q, mu, sigma) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q)) {
    F0 <- stats::pnorm(q, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < ncol(q) || length(sigma) < ncol(q)) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = ncol(q))
      for (i in 1:ncol(q)) {
        F0[, i] <- stats::pnorm(q[, i], mean = mu[i], sd = sigma[i])
      }
    }
  }
  return(F0)
}

#' @name NND
#' @examples
#' dnnd(x = 0.5, mu = 1, sigma = 2)
#'
#' x1 <- c(-0.8, 0.2, 1.6, 3.9)
#' dnnd(x = 0.5, mu = c(1, 1), sigma = c(2, 2))
#'
#' x2 <- matrix(seq(-3, 3, length.out = 10), nrow = 5, ncol = 2)
#' dnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#'
#' x3 <- matrix(seq(-3, 3, length.out = 15), nrow = 5, ncol = 3)
#' dnnd(x3, mu = c(1, 2, 1), sigma = c(2, 2, 2))
#' @export
dnnd <- function(x, mu, sigma) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(x)) {
    df <- stats::dnorm(x, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < ncol(x) || length(sigma) < ncol(x)) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        df[, i] <- stats::dnorm(x[, i], mean = mu[i], sd = sigma[i])
      }
    }
  }
  return(df)
}

#' @name NND
#' @examples
#' q1 <- c(0.01)
#' qnnd(q1, mu = 1, sigma = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export
qnnd <- function(p, mu, sigma) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(mu) < 2 || length(sigma) < 2) {
    qf <- stats::qnorm(p, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < ncol(p) || length(sigma) < ncol(p)) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == ncol(p)) {
        p <- matrix(p, nrow = 1, ncol = ncol(p))
      }
      qf <- matrix(data = NA, nrow = nrow(p), ncol = ncol(p))
      for (i in 1:ncol(p)) {
        qf[, i] <- stats::dnorm(p[, i], mean = mu[i], sd = sigma[i])
      }
    }
  }
  return(qf)
}

#' @name NND
#' @examples
#' n <- 100
#' rnnd(n, mu = c(1, 1), sigma = c(2, 2))
#' rnnd(n, mu = c(1, 2), sigma = c(1, 1))
#' @export
rnnd <- function(n, mu, sigma) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (length(mu) < 2 || length(sigma) < 2) {
    u <- runif(n)
    X <- qnnd(u, mu, sigma)
  } else {
    a <- min(length(mu), length(sigma))
    u <- matrix(runif(n * a), nrow = n, ncol = a)
    X <- qnnd(u, mu[1:a], sigma[1:a])
  }
  return(X)
}
