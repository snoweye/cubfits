# Generate random phi values from lognormal distribution and matched by the
# rank of SCUO indices.

scuo.random <- function(SCUO, phi.Obs = NULL,
    meanlog = .CF.PARAM$meanlog, sdlog = .CF.PARAM$sdlog){
#    meanlog = -1.125, sdlog = 1.5){
#    meanlog = -0.441473, sdlog = 1.393285){
  if(!is.null(phi.Obs)){
    meanlog <- mean(log(phi.Obs))
    sdlog <- sqrt(var(log(phi.Obs)))
  }

  n <- length(SCUO)
  x <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  ret <- sort(x)[rank(SCUO)]

  ret
} # End of scuo.random().
