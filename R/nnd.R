#' Neutrosophic Normal Distribution (NND)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic generalized exponential distribution with parameters mean,
#' \eqn{\mu}, and standard deviation \eqn{\sigma}.
#'
#' The neutrosophic normal distribution with parameters with \code{mean}=\eqn{\mu_N}
#' and \code{std}=\exp{\sigma_N} has density
#' \deqn{f_N(x) = \frac{1}{\sigma_N \sqrt{2 \pi}} e^{\left(\frac{\left(X-\mu_N\right)^2}{2 \sigma_N^2}\right)}}
#' for \eqn{-\infty < x < \infty}, \eqn{\mu_N}, the mean, and \eqn{\sigma_N > 0},
#' the standard deviation.
#'
#' @name NND
#' @param x,q vector or matrix lower and upper of quantiles at which the pdf or cdf needs to be computed.
#' @param p vector or matrix lower and upper of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param mu the value or vector lower and upper of the first parameter.
#' @param sigma the positive value or vector lower and upper of the second parameter.
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
#' p <- seq(0.1, 1, length.out = 21)
#' pnnd(x, mu = c(2,2), sigma = c(1,1))
#'
#' p2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export

pnnd <- function(q, mu, sigma) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q)) {
    F0 <- stats::pnorm(q, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < 2 || length(sigma) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      F0[, 1] <- stats::pnorm(q[, 1], mean = mu[1], sd = sigma[1])
      F0[, 2] <- stats::pnorm(q[, 2], mean = mu[2], sd = sigma[2])
    }
  }
  return(F0)
}

#' @name NND
#' @examples
#' dnnd(x = 0.5, mu = 1, sigma = 2)
#'
#' x1 = c(-0.8,0.2,1.6,3.9)
#' dnnd(x = 0.5, mu = c(1,1), sigma = c(2,2))
#'
#' x2 = matrix(seq(-3,3,length.out = 10), nrow = 2, ncol = 5)
#' dnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export
dnnd <- function(x, mu, sigma) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(x)) {
    df <- stats::dnorm(x, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < 2 || length(sigma) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      df <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      df[, 1] <- stats::dnorm(x[, 1], mean = mu[1], sd = sigma[1])
      df[, 2] <- stats::dnorm(x[, 2], mean = mu[2], sd = sigma[2])
    }
  }
  return(df)
}

#' @name NND
#' @examples
#' q1 <- c(0.01)
#' qnnd(q1, mu = 1, sigma = 2)
#'
#' qnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export
qnnd <- function(p, mu, sigma) {
  if (any(p < 0) || any(p > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(p) && length(mu) < 2 || length(sigma) < 2) {
    qf <- stats::qnorm(p, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < 2 || length(sigma) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(p) && length(p) == 2) {
        p <- matrix(p, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(p), ncol = 2)
      qf[, 1] <- stats::qnorm(p[, 1], mean = mu[1], sd = sigma[1])
      qf[, 2] <- stats::qnorm(p[, 2], mean = mu[2], sd = sigma[2])
    }
  }
  return(qf)
}

#' @name NND
#' @examples
#' n <- 100
#' rnnd(n, mu = c(1,1), sigma = c(2,2))
#' rnnd(n, mu = c(1, 2), sigma = c(1, 1))
#' @export
rnnd <- function(n, mu, sigma) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (length(mu) < 2 || length(sigma) < 2) {
    u <- runif(n)
    X <- qnnd(u, mu, sigma)
  } else {
    u <- matrix(runif(n * 2), nrow = n, ncol = 2)
    X <- qnnd(u, mu, sigma)
  }
  return(X)
}
