#' Neutrosophic Generalized Exponential (NGE)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the
#'
#' @name NGE
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param v the value of the first shape parameter.
#' @param delta the value of the second shape parameter, must be positive, the default is 2.
#'
#' @return  \code{pnge} gives the distribution function,
#'  \code{dnge} gives the density,
#'  \code{qnge} gives the quantile function,
#'  \code{hnge} gives the hazard function,
#'  \code{snge} gives the survival function and
#'  \code{rnge} generates random variables from the Neutrosophic Generalized Exponential (NGE).
#' @references
#'    Rao, Gadde Srinivasa, Mina Norouzirad, and Danial Mazarei. "Neutrosophic Generalized Exponential Distribution with Application." Neutrosophic Sets and Systems 55.1 (2023): 28.
#' @importFrom stats runif uniroot
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pnge(x)
#' pnge(x, v = 2, delta = 1)
#' @export

pnge <- function(x, v=1, delta=2) {
  F0 <- (1 - exp(-x / v))^delta
  return(F0)
}

#' @name NGE
#' @examples
#' dnge(x, v = 1, delta = 2)
#' curve(dnge, -3, 3)
#' @export

dnge <- function(x, v = 1, delta = 2) {
  df <- (delta / v) * (1 - exp(-x / v))^(delta - 1) * exp(-x / v)
  return(df)
}

#' @name NGE
#' @examples
#' qnge(x, v = 1, delta = 3)
#' @export

qnge <- function(q, v = 1, delta = 2) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pnge(t, v, delta)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(0, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#' @name NGE
#' @examples
#' n <- 10
#' rnge(n, v = 2, delta=1)
#' @export
rnge <- function(n, v = 1, delta = 2) {
  u <- runif(n)
  X <- log(-u^(1 / delta) + 1) * (-v)
  return(X)
}


#' @name NGE
#' @examples
#' hnge(x, v = 1, delta = 3)
#' curve(hnge, 0, 3)
#' @export

hnge <- function(x, v = 1, delta = 2) {
  h <- ((delta / v) * (1 - exp(-x / v))^(delta - 1) * exp(-x / v)) / (1 - (1 - exp(-x / v))^delta)
  return(h)
}

#' @name NGE
#' @examples
#' snge(x, v = 1, delta = 3)
#' curve(snge, 0, 3)
#' @export
snge <- function(x, v = 1, delta = 2) {
  s <- 1 - (1 - exp(-x / v))^delta
  return(s)
}
