### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw Gibbs Sampler given current status for measure error (sigmaW) and
### drew inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropType.lognormal_RW <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param, p.DrawScale = 0.1, p.DrawScale.prev = 0.1){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  log.phi.Curr <- log(phi.Curr)

  ### Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  ###                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.phi.Obs - log.phi.Curr)^2) / 2))

  ### Propose sigma.Phi.Curr.
  proplist <- my.propose.sigma.Phi.RW(sigma.Phi.Curr, p.DrawScale[1],
                                      p.DrawScale.prev[1])

  ### M-H step.
  list.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.draw.lognormal.hp.MH(proplist, list.Curr, phi.Curr)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", ret$accept)
  my.update.adaptive("p", ret$accept)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(sigmaW.Curr, ret$nu.Phi, ret$sigma.Phi)
  ret
} # my.pPropType.lognormal_RW().


### Propose sigma.Phi in log scale via random walk.
my.propose.sigma.Phi.RW <- function(sigma.Phi.Curr,
    sigma.Phi.DrawScale = 1, sigma.Phi.DrawScale.prev = 1){
  ### Draw from proposal.
  sigma.Phi.New <- exp(rnorm(1, mean = log(sigma.Phi.Curr),
                                sd = sigma.Phi.DrawScale))
  nu.Phi.New <- -sigma.Phi.New^2 / 2

  # Compute log ratio of prior since lognormal is not symmetric.
  lir <- dlnorm(sigma.Phi.New, meanlog = log(sigma.Phi.Curr),
                sdlog = sigma.Phi.DrawScale, log = TRUE) -
         dlnorm(sigma.Phi.Curr, meanlog = log(sigma.Phi.New),
                sdlog = sigma.Phi.DrawScale.prev, log = TRUE)

  # Return.
  ret <- list(nu.Phi = as.numeric(nu.Phi.New),
              sigma.Phi = as.numeric(sigma.Phi.New),
              lir = lir)
  ret
} # End of my.propose.sigma.Phi.RW().
