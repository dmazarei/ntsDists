### error function
erf <- function(x) 2 * stats::pnorm(x * sqrt(2)) - 1

### inverses error function
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
