#' Neutrosophic Generalized Exponential Distribution (NGED)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the
#'
#' @name NGED
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param v the value of the first shape parameter.
#' @param delta the value of the second shape parameter, must be positive, the default is 2.
#'
#' @return  \code{pnged} gives the distribution function,
#'  \code{dnged} gives the density,
#'  \code{qnged} gives the quantile function,
#'  \code{hnged} gives the hazard function,
#'  \code{snged} gives the survival function and
#'  \code{rnged} generates random variables from the Neutrosophic Generalized Exponential (NGE).
#' @references
#'    Rao, Gadde Srinivasa, Mina Norouzirad, and Danial Mazarei. "Neutrosophic Generalized Exponential Distribution with Application." Neutrosophic Sets and Systems 55.1 (2023): 28.
#' @importFrom stats runif uniroot
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' x2 <- matrix(seq(0, 2, length.out = 40),ncol=2)
#' pnged(x)
#' pnged(x, v = 2, delta = 1)
#' @export

pnged <- function(x, v=1, delta=2) {
  F0 <- (1 - exp(-x / v))^delta
  return(F0)
}

#' @name NGED
#' @examples
#' dnged(x, v = 1, delta = 2)
#' curve(dnged, -3, 3)
#' @export

dnged <- function(x, v = 1, delta = 2) {
  df <- (delta / v) * (1 - exp(-x / v))^(delta - 1) * exp(-x / v)
  return(df)
}

#' @name NGED
#' @examples
#' qnged(x, v = 1, delta = 3)
#' @export

qnged <- function(q, v = 1, delta = 2) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pnged(t, v, delta)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(0, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#' @name NGED
#' @examples
#' n <- 10
#' rnged(n, v = 2, delta=1)
#' @export
rnged <- function(n, v = 1, delta = 2) {
  u <- runif(n)
  X <- log(-u^(1 / delta) + 1) * (-v)
  return(X)
}


#' @name NGED
#' @examples
#' hnged(x, v = 1, delta = 3)
#' curve(hnged, 0, 3)
#' @export

hnged <- function(x, v = 1, delta = 2) {
  h <- ((delta / v) * (1 - exp(-x / v))^(delta - 1) * exp(-x / v)) / (1 - (1 - exp(-x / v))^delta)
  return(h)
}

#' @name NGED
#' @examples
#' snged(x, v = 1, delta = 3)
#' curve(snged, 0, 3)
#' snged(x2, v = c(1,2), delta = c(2,2))
#' @export
snged <- function(x, v = 1, delta = 2) {
  if (is.vector(x)){
    s <- 1 - (1 - exp(-x / v[1]))^delta[1]
  }else{
    if(length(v)<2||length(delta)<2){
      stop(message = "incompatible arguments.")
    }else{
      s <- matrix(data=NA,nrow=nrow(x),ncol=2)
      s[,1] <- 1 - (1 - exp(-x[,1] / v[1]))^delta[1]
      s[,2] <- 1 - (1 - exp(-x[,2] / v[2]))^delta[2]
    }
  }

  return(s)
}

