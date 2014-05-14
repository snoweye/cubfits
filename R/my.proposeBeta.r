### Proposes value from multivariate normal with mean mu, and variance R'R.
### Returns proposal, log importance ratio, and z-values

my.propose.ID_Norm <- function(mu.prev, mu, R){
  ### Draw from proposal.
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R, zProp)

  ### Calculate importance ratio.
  zPrev <- R %*% (mu.prev - mu)
  lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))

  ### Return.
  ret <- list(prop = as.numeric(prop), lir = lir)
  ret
} # End of my.propose.ID_Norm().


my.propose.RW_Norm <- function(mu.prev, mu, R,
    b.DrawScale.aa, b.DrawScale.prev.aa){
  ### Draw from proposal.
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R * b.DrawScale.aa, zProp)

  ### Check if drawing from the same scale.
  if(b.DrawScale.aa == b.DrawScale.prev.aa){
    lir <- 0    # since symmetric random walk.
  } else{
    ### Calculate importance ratio since scaling was changed.
    zPrev <- (R * b.DrawScale.prev.aa) %*% (mu.prev - mu)
    lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))
  }

  ### Return.
  ret <- list(prop = as.numeric(prop), lir = lir)
  ret
} # End of my.propose.RW_Norm().
