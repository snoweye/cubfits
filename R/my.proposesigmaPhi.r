### Propose sigma.Phi in log scale via random walk.
my.proposesigmaPhi.RW_Norm <- function(sigma.Phi.Curr,
    sigma.Phi.DrawScale = .CF.CONF$sigma.Phi.DrawScale,
    sigma.Phi.DrawScale.prev = .CF.CONF$sigma.Phi.DrawScale){
  ### Dispatch.
  log.sigma.Phi.Curr <- log(sigma.Phi.Curr)

  ### Draw from proposal.
  log.sigma.Phi.New <- rnorm(1, mean = log.sigma.Phi.Curr,
                                sd = sigma.Phi.DrawScale)
  sigma.Phi.New <- exp(log.sigma.Phi.New)
  nu.Phi.New <- -sigma.Phi.New^2 / 2

  ### Compute log ratio of prior since lognormal is not symmetric.
  ### This is too slow.
  # lir <- dlnorm(sigma.Phi.New, meanlog = log.sigma.Phi.Curr,
  #               sdlog = sigma.Phi.DrawScale, log = TRUE) -
  #        dlnorm(sigma.Phi.Curr, meanlog = log.sigma.Phi.New,
  #              sdlog = sigma.Phi.DrawScale.prev, log = TRUE)

  ### Faster since the next relations of normal and log normal
  ### x <- 1.5; m <- 2; s <- 3
  ### dnorm(log(phi), m, s, log = TRUE) - log(phi) ==
  ###   dlnorm(phi, m, s, log = TRUE)
  lir <- -log.sigma.Phi.New + log.sigma.Phi.Curr    # Jacobin
  if(sigma.Phi.DrawScale.prev != sigma.Phi.DrawScale){
    lir <- lir + dnorm(log.sigma.Phi.New, log.sigma.Phi.Curr,
                       sigma.Phi.DrawScale, log = TRUE) -
                 dnorm(log.sigma.Phi.Curr, log.sigma.Phi.New,
                       sigma.Phi.DrawScale.prev, log = TRUE)
  }

  ### Return.
  ret <- list(nu.Phi = as.numeric(nu.Phi.New),
              sigma.Phi = as.numeric(sigma.Phi.New),
              lir = lir)
  ret
} # End of my.proposesigmaPhi.RW_Norm().