### error function
erf <- function(x) {
  2 * stats::pnorm(x * sqrt(2)) - 1
}

### inverses error function
inv_erf <- function (x) {
  stats::qnorm((1 + x)/2)/sqrt(2)
}
