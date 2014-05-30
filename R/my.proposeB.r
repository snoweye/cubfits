### Proposes value from multivariate normal with mean mu, and variance R'R.
### Returns proposal, log importance ratio, and z-values

my.proposeB.ID_Norm <- function(mu.prev, mu, R){
  ### Draw from proposal.
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R, zProp)

  ### Calculate importance ratio.
  zPrev <- R %*% (mu.prev - mu)
  lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))

  ### Return.
  ret <- list(prop = as.numeric(prop), lir = lir)
  ret
} # End of my.proposeB.ID_Norm().


### In this case, mu.prev == mu since my.drawBConditionalAll.RW_Norm() set this.
my.proposeB.RW_Norm <- function(mu.prev, mu, R,
    b.DrawScale.aa, b.DrawScale.prev.aa){
  ### Draw from proposal.
  zProp <- rnorm(length(mu))
  prop <- mu + backsolve(R / b.DrawScale.aa, zProp)

  ### Check if drawing from the same scale.
  lir <- 0    # no jacobin since no transformation
  if(b.DrawScale.aa != b.DrawScale.prev.aa){
    ### Calculate importance ratio since random walk scale was changed.
    zPrev <- (R / b.DrawScale.prev.aa) %*% (mu.prev - prop)
    lir <- -1/2 * (sum(zProp * zProp) - sum(zPrev * zPrev))
  }

  ### Return.
  ret <- list(prop = as.numeric(prop), lir = lir)
  ret
} # End of my.proposeB.RW_Norm().
