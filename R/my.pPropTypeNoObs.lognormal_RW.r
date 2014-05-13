### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropTypeNoObs.lognormal_RW <- function(n.G, phi.Curr,
    p.Curr, hp.param, p.DrawScale = 0.1, p.DrawScale.prev = 0.1){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[1]
  sigma.Phi.Curr <- p.Curr[2]

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
  ret <- c(ret$nu.Phi, ret$sigma.Phi)
  ret
} # my.pPropTypeNoObs.lognormal_RW().
