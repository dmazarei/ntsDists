#' Neutrosophic Exponential Distribution (NED)
#'
#' Density, distribution function, quantile function and random generation for
#' the nuetrosophic exponential distribution with the parameter \code{rate}=\eqn{\theta_N}.
#'
#' The neutrosophic exponential distribution with parameters \code{rate}=\eqn{\theta_N}
#' has density
#' \deqn{f_N(x)=\theta_N \exp \left(-x \theta_N\right)}
#' for \eqn{x \ge 0} and \eqn{\theta_N > 0}, the rate parameter.
#' @name NED
#' @param x a vector or matrix of observations at which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles at which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities at which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param theta the shape parameter must be a positive interval.
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X \ge x)}; otherwise, \eqn{P(X >x)}.
#'
#' @return
#'  \code{pned} gives the distribution function,
#'  \code{dned} gives the density,
#'  \code{qned} gives the quantile function and
#'  \code{rned} generates random variables from the neutrosophic exponential distribution.
#'
#' @references
#' Duan, W., Q., Khan, Z., Gulistan, M., Khurshid, A. (2021). Neutrosophic
#' Exponential Distribution: Modeling and Applications for Complex Data Analysis,
#' \emph{Complexity}, 2021, 1-8.
#'
#' @importFrom stats runif
#'
#' @examples
#' # Example 4 of Duan et al. (2021)
#' data <- matrix(c(4,4,3.5,3.5,3.9,4.1,4.2,4.2,4.3,4.6, 4.7, 4.7), nrow = 6, ncol = 2, byrow = TRUE)
#' # Note that we have a matrix of nuetrosophic observations
#'
#' # The density function of data with the estimated value for parameters based
#' # on Duan et al. (2021)
#' dned(data, theta = c(0.23, 0.24))
#' # The density function for the nuetrosophic observation (4,4.1)
#' dned(x = c(4,4.1), theta = c(0.23, 0.24))
#'
#' # The density function for the nuetrosophic observation 4
#' # Here, 4 is equivalent to c(4,4).
#' dned(4, theta = c(0.23, 0.24))
#' @export
dned <- function(x, theta) {
  if (any(theta <= 0) || any(x < 0))
    stop("Arguments are incompatible.")

  theta <- rep(theta, length.out = 2)
  if(is.vector(x)){
    x <- matrix(rep(x, length.out = 2), ncol = 2)
    pdf <- theta * exp(-x * theta)
   }

  pdf <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  for (i in 1:ncol(x)) {
    pdf[, i] <- theta[i] * exp(-x[, i] * theta[i])
  }

  # Identify rows where col1 > col2
  swap_rows <- pdf[, 1] > pdf[, 2]
  # Swap values using logical indexing
  pdf[swap_rows, c(1, 2)] <- pdf[swap_rows, c(2, 1)]

  return(pdf)
}
#' @name NED
#' @examples
#'
#' # The cumulative distribution function for the nuetrosophic observation (4,4.1)
#' pned(c(4,4.1), theta = c(0.23, 0.24), lower.tail = TRUE)
#' pned2(c(4,4.1), theta = c(0.23, 0.24), lower.tail = TRUE)
#'
#'
#' pned(4, theta = c(0.23, 0.24))
#' pned1(4, theta = c(0.23, 0.24))
#' pned2(4, theta = c(0.23, 0.24))
#' @export
pned <- function(q, theta, lower.tail = TRUE) {
  if (any(theta <= 0) || any(q < 0))
    stop("Arguments are incompatible.")

  theta <- rep(theta, length.out = 2)
  if (is.vector(q)){
    q <- rep(q, length.out = 2)
  }
  cdf <- 1 - exp(-q * theta)

  if (!lower.tail)
    cdf <- 1 - cdf

  return(cdf)
}
#' @name NED
#' @examples
#' The first percentile
#' qned(p = 0.1, theta = 0.25)
#'
#' The quantiles
#' qned(p = c(0.25, 0.5, 0.75), theta = c(0.24, 0.25))
#'
#'
#'
#' @export
qned <- function(p, theta) {
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p should be in the interval [0,1].")
  }

  if (any(theta <= 0)){
    stop(message = "Arguments are incompatible.")
  }

  theta <- rep(theta, length.out = 2)
  p <- matrix(rep(p, each = 2), ncol = 2, byrow = TRUE)

  quantiles <- matrix(data = NA, nrow = nrow(p), ncol = 2)
  for (i in 1:ncol(p)) {
    quantiles[, i] <- -log(1 - p[, i]) / theta[i]
  }

  # Identify rows where col1 > col2
  swap_rows <- quantiles[, 1] > quantiles[, 2]
  # Swap values using logical indexing
  quantiles[swap_rows, c(1, 2)] <- quantiles[swap_rows, c(2, 1)]

return(quantiles)
}

#' @name NED
#' @examples
#' Simulate 10 numbers
#' rned(n, theta = c(0.23, 0.24))
#' @export
#'
rned <- function(n, theta) {
  if (any(theta <= 0))
    stop(message = "Arguments are incompatible.")
  theta <- rep(theta, length.out = 2)
  u <- matrix(runif(n * length(theta)), nrow = n, ncol = length(theta))
  X <- qned(u, theta)

  return(X)
}
