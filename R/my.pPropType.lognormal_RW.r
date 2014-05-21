### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw Gibbs Sampler given current status for measure error (sigmaW) and
### drew inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropType.lognormal_RW <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  log.phi.Curr <- log(phi.Curr)
  p.DrawScale <- .cubfitsEnv$all.DrawScale$p[1]
  p.DrawScale.prev <- .cubfitsEnv$all.DrawScale$p.prev[1]

  ### Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  ###                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  sigmaW.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.phi.Obs - log.phi.Curr)^2) / 2))

  ### Propose sigma.Phi.Curr.
  proplist <- my.propose.sigma.Phi.RW(
                sigma.Phi.Curr,
                sigma.Phi.DrawScale = p.DrawScale,
                sigma.Phi.DrawScale.prev = p.DrawScale.prev)

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
    sigma.Phi.DrawScale = .CF.CONF$sigma.Phi.DrawScale,
    sigma.Phi.DrawScale.prev = .CF.CONF$sigma.Phi.DrawScale){
  ### Dispatch.
  log.sigma.Phi.Curr <- log(sigma.Phi.Curr)

  ### Draw from proposal.
  log.sigma.Phi.New <- rnorm(1, mean = log.sigma.Phi.Curr,
                                sd = sigma.Phi.DrawScale)
  sigma.Phi.New <- exp(log.sigma.Phi.New)
  nu.Phi.New <- -sigma.Phi.New^2 / 2

  ### Check if drawing from the same scale.
  lir <- 0    # since symmetric random walk.
  if(sigma.Phi.DrawScale != sigma.Phi.DrawScale.prev){
    ### Calculate importance ratio since random walk scale was changed.
    lir <- dnorm(log.sigma.Phi.New, log.sigma.Phi.Curr,
                 sd = sigma.Phi.DrawScale) -
           dnorm(log.sigma.Phi.Curr, log.sigma.Phi.New,
                 sd = sigma.Phi.DrawScale.prev)
  }

  ### Compute log ratio of prior since lognormal is not symmetric.
  ### Is this wrong? Or, do I assume sigma.Phi ~ LN(?, ?).
  # lir <- dlnorm(sigma.Phi.New, meanlog = log(sigma.Phi.Curr),
  #               sdlog = sigma.Phi.DrawScale, log = TRUE) -
  #        dlnorm(sigma.Phi.Curr, meanlog = log(sigma.Phi.New),
  #               sdlog = sigma.Phi.DrawScale.prev, log = TRUE)

  ### Return.
  ret <- list(nu.Phi = as.numeric(nu.Phi.New),
              sigma.Phi = as.numeric(sigma.Phi.New),
              lir = lir)
  ret
} # End of my.propose.sigma.Phi.RW().
