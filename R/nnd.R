#' Neutrosophic Normal Distribution (NND)
#'
#' Computes the pdf, cdf, quantile and random numbers of the nuetrosophic generalized exponential distribution.
#' \deqn{X_N \sim N_N\left(\mu_N, \sigma_N\right)=\frac{1}{\sigma_N \sqrt{2 \pi}} e^{\left(\frac{\left(X-\mu_N\right)^2}{2 \sigma_N^2}\right)}}
#' for   \eqn{\mu_N}, the first parameter, and \eqn{\sigma_N > 0}, the second parameter.
#' @name NND
#' @param x scaler or vector or matrix lower and upper of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param mu the value or vector lower and upper of the first parameter.
#' @param sigma the value or vector lower and upper of the second parameter, must be positive.
#'
#' @return  \code{pnnd} gives the distribution function,
#'  \code{dnnd} gives the density,
#'  \code{qnnd} gives the quantile function and
#'  \code{rnnd} generates random variables from the Neutrosophic Normal Distribution (NND).
#' @references
#'    Patro, S. and F. Smarandache, The Neutrosophic Statistical Distribution, More Problems, More Solutions. 2016: Infinite Study.
#' @importFrom stats runif dnorm pnorm qnorm
#' @examples
#' x <- seq(0.1, 1, length.out = 21)
#' x2 <- matrix(seq(0.1, 2, length.out = 40), ncol = 2)
#' pnnd(x)
#' pnnd(x, mu = 2, sigma = 1)
#' pnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export

pnnd <- function(x, mu = 1, sigma = 2) {
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(x)) {
    F0 <- stats::pnorm(x, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < 2 || length(sigma) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      F0 <- matrix(data = NA, nrow = nrow(x), ncol = 2)
      F0[, 1] <- stats::pnorm(x[, 1], mean = mu[1], sd = sigma[1])
      F0[, 2] <- stats::pnorm(x[, 2], mean = mu[2], sd = sigma[2])
    }
  }
  return(F0)
}

#' @name NND
#' @examples
#' dnnd(x, mu = 1, sigma = 2)
#' curve(dnnd, .1, 3)
#' dnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export
dnnd <- function(x, mu = 1, sigma = 2) {
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
#' qnnd(x, mu = 1, sigma = 2)
#' x2 <- matrix(seq(0.1, 1, length.out = 40), ncol = 2)
#' qnnd(x2, mu = c(1, 2), sigma = c(2, 2))
#' @export
qnnd <- function(q, mu = 1, sigma = 2) {
  if (any(q < 0) || any(q > 1)) stop(message = "[Warning] 0 < x < 1.")
  if (any(sigma <= 0)) stop(message = "incompatible arguments.")
  if (is.vector(q) && length(mu) < 2 || length(sigma) < 2) {
    qf <- stats::qnorm(q, mean = mu[1], sd = sigma[1])
  } else {
    if (length(mu) < 2 || length(sigma) < 2) {
      stop(message = "incompatible arguments.")
    } else {
      if (is.vector(q) && length(q) == 2) {
        q <- matrix(q, nrow = 1, ncol = 2)
      }

      qf <- matrix(data = NA, nrow = nrow(q), ncol = 2)
      qf[, 1] <- stats::qnorm(q[, 1], mean = mu[1], sd = sigma[1])
      qf[, 2] <- stats::qnorm(q[, 2], mean = mu[2], sd = sigma[2])
    }
  }
  return(qf)
}

#' @name NND
#' @examples
#' n <- 100
#' rnnd(n, mu = 2, sigma = 1)
#' rnnd(n, mu = c(1, 2), sigma = c(1, 1))
#' @export
rnnd <- function(n, mu = 1, sigma = 2) {
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
