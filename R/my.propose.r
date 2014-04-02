# Proposes value from multivariate normal with mean mu, and variance R'R.
# Returns proposal, log importance ratio, and z-values

my.propose.Norm <- function(prev, mu, R){
  # Draw from proposal
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R, zProp)

  # Calculate importance ratio
  zPrev <- R %*% (prev - mu)
  lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))
  ret <- list(prop = as.numeric(prop), lir = lir,
              zPrev = zPrev, zProp = as.numeric(zProp))
  ret
} # End of my.propose.Norm().


my.propose.RWNorm <- function(prev, R = NULL, A = NULL, sFactor = 1){
  # Draw from proposal. prev is PREVious value.
  zProp <- rnorm(length(prev))
  if(!is.null(A)){
    prop <- prev + (A %*% zProp) * sFactor
  } else if(!is.null(R)){
    prop <- prev + backsolve(R, zProp) * sFactor
  } else{
    prop <- prev + zProp * sFactor
  }
  ret <- list(prop = as.numeric(prop), lir = 0, zProp = as.numeric(zProp))
  ret
} # End of my.propose.RWNorm().

