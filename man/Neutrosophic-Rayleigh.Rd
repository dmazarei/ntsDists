% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsRayleigh.R
\name{Neutrosophic Rayleigh}
\alias{Neutrosophic Rayleigh}
\alias{dnsRayleigh}
\alias{pnsRayleigh}
\alias{qnsRayleigh}
\alias{rnsRayleigh}
\title{Neutrosophic Rayleigh Distribution}
\usage{
dnsRayleigh(x, theta)

pnsRayleigh(q, theta, lower.tail = TRUE)

qnsRayleigh(p, theta)

rnsRayleigh(n, theta)
}
\arguments{
\item{x}{a vector or matrix of observations for which the pdf needs to be computed.}

\item{theta}{the shape parameter, which must be a positive interval.}

\item{q}{a vector or matrix of quantiles for which the cdf needs to be computed.}

\item{lower.tail}{logical; if TRUE (default), probabilities are
\eqn{P(X \leq x)}; otherwise, \eqn{P(X >x)}.}

\item{p}{a vector or matrix of probabilities for which the quantile needs to be computed.}

\item{n}{number of random values to be generated.}
}
\value{
\code{dnsRayleigh} gives the density function

 \code{pnsRayleigh} gives the distribution function

 \code{qnsRayleigh} gives the quantile function

 \code{rnsRayleigh} generates random variables from the Neutrosophic Rayleigh Distribution.
}
\description{
Density, distribution function, quantile function and random
generation for the neutrosophic Rayleigh distribution with
parameter \eqn{\theta_N}.
}
\details{
The neutrosophic Rayleigh distribution with parameter \eqn{\theta_N}
has the density
\deqn{f_N(x)=\frac{x}{\theta_N^2} \exp\{-\frac{1}{2}\left(\frac{x}{\theta_N}\right)^2\}}
for  \eqn{\theta_N \in (\theta_L, \theta_U)}, which must be a positive
interval and \eqn{x \ge 0}.
}
\examples{
data(remission)
dnsRayleigh(x = remission, theta = c(9.6432, 9.8702))

pnsRayleigh(q = 20, theta = c(9.6432, 9.8702))

# Calculate quantiles
qnsRayleigh(p = c(0.25, 0.5, 0.75), theta = c(9.6432, 9.8702))

# Simulate 10 values
rnsRayleigh(n = 10, theta = c(9.6432, 9.8702))

}
\references{
Khan, Z., Gulistan, M., Kausar, N. and Park, C. (2021).
Neutrosophic Rayleigh Model With Some Basic Characteristics and
Engineering Applications, in \emph{IEEE Access}, 9, 71277-71283.
}
